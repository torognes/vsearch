# Design doc — unify the per-worker self-scheduling loop

**Status:** IMPLEMENTED (branch `tmp_20260710155615`).
**Date:** 2026-07-10
**Origin:** the note *"extend usage of thread pool"*.

## Implementation summary

`run_worker_loop()` was added in `src/utils/worker_loop.hpp` and all 8
pull-loop sites were converted, one command per commit: mask, allpairs,
search (CLI + batch), chimera (CLI + batch), search_exact, sintax. Per
review decisions: the two lambdas are **named** (`has_work_to_claim` +
a command-specific `process_sequence` / `process_query`), the batch API
sites were included, and the scaffolding comments (the RAII-unlock notes
and the `/* let other threads read input */` markers) were removed since
the helper's header documents the lock discipline. The `cluster` static
slicing and the `fastq_mergepairs` pipeline were left out of scope as
planned.

**Verification:** full debug build + release build clean (no new
warnings); `cppcheck` on all six files shows only pre-existing findings;
the 6 CLI-path sites are byte-identical to the `dev` baseline at
`--threads 1` (sintax with a fixed `--randseed`); the 2 batch sites pass
the api_examples "batch X matches sequential" checks; the full
vsearch-tests suite is green (9648 PASS, 0 FAIL).

---

## 1. Context: what is already done

The literal note is **already implemented**. `ThreadRunner`
(`utils/threads.hpp`) is a RAII fork-join pool that replaced every
hand-rolled `pthread` worker pool. Eight migration commits
(`9b47d923`, `ec0a1e85`, `d156b660`, `bfa6a05f`, `f1ab1eea`,
`a912e866`, `9d4212bd`, `4a0c4dbc`) plus the removal of the `xpthread`
layer (`83c98ac6`) finished the job. There is **no cruder
thread-lifecycle model left to migrate**: no raw `std::thread` outside
`threads.hpp`, no `pthread` (only a comment in `chimera.cc`), no
OpenMP.

So this doc is *not* about `ThreadRunner`. It is about the layer
**on top of** it: the work-distribution loop each command hand-rolls
inside its worker body. That code is duplicated and is where the
"less flexible / cruder" feeling actually lives.

## 2. The duplicated pattern

`ThreadRunner` deliberately owns no task queue — it runs one function
once per thread. Every command therefore embeds the same
self-scheduling loop in its worker body:

```
while (true) {
  lock(input_mutex);
  if (claim_next())          // fastx_next()  OR  index = counter++
    { copy/read work; unlock(input_mutex); do_work(); [lock output; commit] }
  else
    { unlock(input_mutex); break; }
}
```

### Inventory of the 8 sites

| # | Site | Input source | Output sync | Per-thread scratch |
|---|------|--------------|-------------|--------------------|
| 1 | `search.cc` `search_thread_run` | shared fastx reader | `mutex_output` (queries, qmatches, progress) | `si_plus[t]`, `si_minus[t]` |
| 2 | `search_exact.cc` `search_exact_thread_run` | shared fastx reader | `mutex_output` | `si_plus[t]`, `si_minus[t]` |
| 3 | `sintax.cc` `sintax_thread_run` | shared fastx reader | `mutex_output` | `si_plus[t]`, `si_minus[t]` |
| 4 | `chimera.cc` `chimera_thread_core` | fasta reader **or** DB counter (ref vs denovo) | `mutex_output` | `chimera_info_s * ci` + local aligner |
| 5 | `allpairs.cc` `allpairs_thread_run` | DB counter (triangular) | `mutex_output`, **non-linear progress** | per-thread `searchinfo` + local vectors |
| 6 | `mask.cc` `dust_all_worker` | DB counter | `state.mutex` (progress only) | none (works on DB in place) |
| 7 | `search.cc` `search_batch_worker_fn` | counter `next_query` | **none** — disjoint `results[qi]` slice | `batch_si_plus[tid]` |
| 8 | `chimera.cc` `chimera_batch_worker_fn` | counter `next_query` | **none** — disjoint `results[qi]` slice | `ci_array[tid]` |

Two sub-families fall out:

- **A. Index counter** (5, 6, 7, 8, and chimera-denovo): claim an
  integer index under the lock; the payload is read afterward from the
  immutable loaded DB / caller arrays.
- **B. Stateful reader** (1, 2, 3, chimera-ref): `fastx_next` /
  `fasta_next` mutates a shared handle, so the record **must be copied
  into per-thread scratch while still holding the lock** — another
  worker overwrites the reader's buffers on its next pull.

The batch sites (7, 8) are a clean variant: counter input **and no
output mutex at all**, because each query writes to its own disjoint
`results[qi]` slot.

## 3. Proposed helper

The only thing genuinely common is the *loop skeleton and the input
mutex discipline*. Everything else (what "claim" means, what "work"
means, whether output is locked) is command-specific. So the helper
should own **only** the skeleton and hand two callables back to the
caller.

Header-only, C++11, sits next to `ThreadRunner` (either in
`utils/threads.hpp` or a new `utils/worker_loop.hpp`):

```cpp
/*
  Drives one worker's self-scheduling loop. `claim` runs under
  `input_mutex` and returns true after it has claimed/copied the next
  unit of work into thread-local scratch, or false when the input is
  exhausted. `work` then runs WITHOUT the lock held; it performs the
  computation and does its own output synchronisation (if any).

  Centralising the loop guarantees the input lock is always released
  before `work()` runs — the property every hand-rolled copy currently
  re-establishes by careful placement of `input_lock.unlock()`.
*/
template <typename ClaimFn, typename WorkFn>
auto run_worker_loop(std::mutex & input_mutex,
                     ClaimFn claim,
                     WorkFn work) -> void {
  while (true) {
    {
      std::unique_lock<std::mutex> const lock(input_mutex);
      if (not claim()) { break; }   // claim under lock; false = done
    }
    work();                          // no lock held
  }
}
```

`claim` and `work` are lambdas that capture the enclosing worker's
`t` (thread id), `state`, and per-thread scratch — exactly what the
bodies already capture. The helper deliberately does **not** template
on a "work item" type or own an output mutex: sites differ too much
(reader vs counter, output-mutex vs lock-free), and forcing a common
item type would be more coupling, not less.

### Before / after (site 6, `mask.cc`, the simplest)

Before:
```cpp
static auto dust_all_worker(struct dust_state_s & state, uint64_t) -> void {
  while (true) {
    std::unique_lock<std::mutex> lock(state.mutex);
    const auto seqno = state.nextseq;
    if (seqno < state.seqcount) {
      ++state.nextseq;
      state.progress->update(seqno);
      lock.unlock();
      dust(db_getsequence(seqno), ..., *state.parameters);
    } else {
      break;
    }
  }
}
```

After:
```cpp
static auto dust_all_worker(struct dust_state_s & state, uint64_t) -> void {
  uint64_t seqno = 0;
  run_worker_loop(state.mutex,
    /*claim*/ [&]() -> bool {
      if (state.nextseq >= state.seqcount) { return false; }
      seqno = state.nextseq++;
      state.progress->update(seqno);
      return true;
    },
    /*work*/ [&]() {
      dust(db_getsequence(seqno), ..., *state.parameters);
    });
}
```

## 4. Which sites fit, and the edge cases

- **All 8 sites fit** the two-callable shape.
- **Chimera dual input source (4).** The `uchime_ref` (reader) vs
  denovo (counter) branch lives *inside* `claim()`. Fits cleanly; no
  special support needed.
- **Chimera progress-position race (comment at `chimera.cc:2176`).**
  The fix reads `fasta_get_position()` under the input lock. With the
  helper that read stays inside `claim()`, i.e. still under the lock —
  the fix is preserved by construction.
- **Allpairs non-linear progress (5).** `progress += seqcount -
  query_no - 1` stays in `work()` under the output lock. Unaffected.
- **Batch sites have no output mutex (7, 8).** `work()` simply doesn't
  lock; writes to the disjoint `results[qi]`. The helper never forces
  an output mutex, so this is natural.
- **Reader sites must copy under the lock (B).** The `strcpy` into
  `si_plus[t]` / `populate_si` / `realloc_arrays` all move into
  `claim()`, which runs under the lock — same as today.
- **No worker throws / no `std::exit` from a worker.** History shows
  these were deliberately removed (cooperative abort in mergepairs,
  `5355315c`/`95b70429`). The helper adds no exception path; `claim`
  and `work` behave as they do now.

## 5. Explicitly out of scope

- **`cluster.cc` (`cluster_work_pool_s`).** Uses *static* per-thread
  slicing per round, not a pull-loop. Converting it to dynamic
  scheduling is a **behaviour/perf change** (better load-balancing only
  if per-query cost varies within a round) and needs hyperfine
  benchmarks. Not part of this refactor.
- **`fastq_mergepairs.cc`.** Roles are assigned by thread index
  (read / process / write pipeline over a chunk ring buffer) to overlap
  I/O with compute. This is a deliberate pipeline, not a parallel-for;
  genericising it would be a risky rewrite with dubious benefit. Leave
  as is.

## 6. Rollout & verification

Incremental, one commit per site, lowest-risk first:

1. Sites 7, 8 (batch API, counter-only, no output lock) — simplest.
2. Site 6 (`mask`).
3. Site 5 (`allpairs`).
4. Sites 1–3 (reader family).
5. Site 4 (`chimera`, dual source) last.

Each step is a pure structural rewrite — **no observable behaviour
change**. Verify per step:

- `--threads 1` vs `--threads 8` output **byte-identical** (the
  standing cluster/search invariant).
- Relevant vsearch-tests suites green.
- Full suite + ASan/UBSan on the library examples before the final
  merge.
- `cppcheck` + `clang-tidy` clean on touched files.

## 7. Honest cost / benefit

- **Benefit:** removes 8 copies of the loop/lock scaffolding; makes the
  "release input lock before doing work" rule a property of one tested
  helper instead of something each site re-derives; matches the
  CLAUDE.md goals (complexity reduction, decoupling, RAII). Low risk.
- **Limit:** the scaffolding is only ~6 lines per site; the large,
  heterogeneous bodies stay in the `claim`/`work` lambdas. This is a
  **uniformity + correctness-guard** win, not a big LOC reduction. It
  does not make the pool "more flexible" in a functional sense — vsearch
  has no need for heterogeneous task queues, so a full task-pool would
  be over-engineering.

## 8. Open questions for review

1. Worth doing at all, given the modest LOC payoff? (I lean **yes** for
   the correctness-guard + uniformity, but it is a judgement call.)
2. Helper location: extend `utils/threads.hpp`, or a new
   `utils/worker_loop.hpp`?
3. Naming: `run_worker_loop`? `dispatch_worker`? `claim_loop`?
4. Do you want the batch sites (7, 8) included, or keep the
   library-API path untouched for ABI-stability caution?
