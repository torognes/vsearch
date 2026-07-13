# Remaining global variables: inventory and elimination plan

Status: in progress (Phase 1 + 2 landing on branch `tmp_20260713082154`).
Date: 2026-07-13.

Goal (from `CLAUDE.md`): "avoid non const global variables". This document
inventories every mutable global with static storage duration that remains in
`src/`, and lays out a staged plan to eliminate them, ordered easy/low-risk to
hard/high-risk.


## Method

The inventory was produced by cross-referencing `nm -C` over all built object
files (which lists every static-storage-duration symbol landing in `.data` /
`.bss`, including anonymous-namespace and file-`static` state that a source
grep alone misses) against the source declarations. Read-only symbols (`r`/`R`
sections) and function-local `const` statics were filtered out: they are not
mutable globals and are left alone.


## Inventory: the 34 mutable globals that remain

| Group | File | Variables | Linkage | Notes |
|-------|------|-----------|---------|-------|
| **A** | `core/db.cpp` | `datap`, `seqindex` (both `extern` in `db.hpp`), `seqindex_alloc`, `h`, `is_fastq`, `sequences`, `nucleotides`, `longest`, `shortest`, `longestheader`, `dataalloc`, `datalen` — **12** | 2 external, 10 file-static | The in-memory sequence database. Already wrapped by free-function + inline accessors (`db_get*`); **31 TUs** consume it, but only `udb.cpp` touches `datap`/`seqindex` directly. |
| **B** | `core/dbhash.cpp` | `dbhash_table` (external, no header decl), `dbhash_bitmap`, `dbhash_size`, `dbhash_shift`, `dbhash_mask` — **5** | 1 external-leak, 4 static | DB dedup hash index; tightly coupled to Group A. |
| **C** | `core/mergepairs.cpp` | tables `merge_qual_same`, `merge_qual_diff`, `match_score`, `mism_score`, `q2p` (all external-linkage, no header decl); abort state `merge_abort`, `merge_error_claimed`, `merge_error_reason`, `merge_error_value` — **9** | 5 external-leak, 4 static | Tables filled once from `Parameters` at init; abort state coordinates worker threads (exposed via `request_merge_abort`/`merge_aborted`). |
| **D** | `core/otutable.cpp` | `otutable` (pointer to heap `otutable_s`) — **1** | static | Singleton; the struct already exists. |
| **E** | `core/getseq.cpp` | `labels_data` (external, no header decl) — **1** | external-leak | `std::vector<std::vector<char>>`. |
| **F** | `core/derep.cpp`, `core/kmerhash.cpp`, `commands/derep_smallmem.cpp` | `hash_function` ×3; plus `hashtable`, `hashtablesize` (derep_smallmem) — **5** | static | The three `hash_function` pointers are never reassigned -> const-able immediately. |
| **G** | `utils/random.cpp` | `base_seed` — **1** | static | RNG seed, set from `--randseed` at init, read via `random_base_seed()`. |

### Latent defect worth flagging

The Group C tables, plus `dbhash_table` (B) and `labels_data` (E), have
**external linkage with no `extern` declaration** anywhere. They pollute the
global symbol namespace and are one accidental duplicate name away from an
ODR/link collision — `q2p` in particular already coexists with unrelated
external `q2p` *functions* in `commands/fastq_stats.cpp` and
`commands/fastq_eestats2.cpp`. Giving them internal linkage is a correctness
fix, not just hygiene.


## Not targets (documented, left alone)

- `vsearch.cc:session_mutex` — intentional process-wide `std::mutex`
  serialising the C-API `vsearch_session_begin`/`vsearch_session_end`. A
  legitimate singleton.
- `default_quality_padding` / `alternative_quality_padding` /
  `default_sequence_padding` (`vsearch.h`) — `std::string const`, so already
  **const** (not mutable). Being `std::string const` in a header gives them
  internal linkage, i.e. one heap-allocating copy per TU. Converting to
  `constexpr char[]` is optional hygiene, deferred (touches a widely-included
  header; not a mutable-global concern).
- `bz2_libname`/`gz_libname` (const), `userfields_names`/`help_message`
  (const), `fastx.cpp:blanks` (function-local `const` statics) — all const.


## Plan

### Phase 1 — const-ify the set-once (no ABI impact)  [this branch]

- Group F: the three `hash_function` function pointers are never reassigned ->
  `static constexpr Hash hash_function = ...` (a function address is a valid
  constant-expression initializer). Removes them from the mutable-global set.

### Phase 2 — plug the external-linkage leaks (low risk)  [this branch]

Give internal linkage to the file-scope state that currently leaks external
symbols, matching the `static` convention already used by sibling declarations
in each file. Independently commit-worthy; also de-risks Phase 3/4.

- `core/mergepairs.cpp`: `static` on `merge_qual_same`, `merge_qual_diff`,
  `match_score`, `mism_score`, `q2p` (the abort state below them is already
  `static`).
- `core/dbhash.cpp`: `static` on `dbhash_table` (siblings already `static`).
- `core/getseq.cpp`: `static` on `labels_data`.

### Phase 3 — encapsulate self-contained TU state into per-run structs (medium)

Each is independent and low-blast-radius; thread an owned instance by reference
through the relevant functions (the proven `Dbindex` RAII-struct playbook).

- **D** `otutable` -> own an `Otutable` in the command entry point; thread
  `Otutable&` through the `otutable_*` functions.
- **C** mergepairs -> bundle tables + abort state into a `Merge_context` owned
  by `pair_all()`; replace the abort free-functions with methods / parameter
  passing. Must preserve the atomics and the join happens-before that currently
  avoid the `std::exit()`-during-worker data race documented in the file.
- **E** `labels_data` -> thread through the `getseq` functions.
- **F** `hashtable`/`hashtablesize` -> fold into the derep_smallmem state struct
  (aligns with the pending derep refactor, `TBD_20260712_derep_refactoring.md`).
- **G** `base_seed` -> a small seed holder threaded to `random_*`.

### Phase 4 — the DB core (large, ABI-breaking; do last)

Groups A + B -> a `Db` type holding `datap`/`seqindex`/counts with `db_get*`
as members (or free functions taking `Db const &`), and `Dbhash` alongside it.
Mirrors the completed `Dbindex` refactor exactly: instantiate in the command
layer, thread a mutable `Db &`, update the library API, and break the
`libvsearch_core` ABI per the standing E1 decision. Blast radius: 31 TUs +
`udb.cpp` + public API.

### Cross-cutting risks

- **Runtime-mutated-globals trap**: DB state is mutated after load by
  `udb_read` and the `db_sortby*` functions — Phase 4 must keep a *mutable*
  `Db &` reaching those mutation points, not a `const` snapshot.
- **Threading**: Group C's abort state is deliberately atomic file-scope state
  to avoid a `std::exit()`-during-worker data race; any `Merge_context` must
  preserve the atomics and the happens-before on join.
- Do Phase 4 on its own `tmp_` branch with the full test suite plus a
  release-build library-API run (debug libs segfault the `api_examples` net).
