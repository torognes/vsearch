# Remaining global variables: inventory and elimination plan

Status: Phases 1, 2, and 4 **DONE and merged into `dev`**. Phase 3 is now
**mostly done**: Group D (`otutable` → the owned `OtuTable` class), the Group C
quality-score tables (→ the `QualityTables` struct), and all of Group F
(`hash_function` const-ified + `hashtable`/`hashtablesize` localised) have all
landed on `dev`. What is left is a small residue of genuinely-mutable file-scope
state: **only Group C-abort** (the 4 merge worker-coordination atomics). Groups E
and G from the original set and all three re-sweep additions (showalign buffers,
`derep_sort_db`, `unique.cpp`'s `hash_function`) have since been eliminated.
Original inventory: 2026-07-13. Revised 2026-07-18 (re-swept on `dev` @
`37c2b5d7`). Revised 2026-07-21: `derep_sort_db` eliminated (commit `090751fe`),
Group G (`base_seed`) eliminated via the owned `RandomSeed` class, then Group E
(`labels_data`) and the showalign row buffers eliminated (threaded as a local /
bundled into an `AlignmentRows` struct). An `nm` re-sweep over the current build
confirms the mutable-global residue is now **4** with no new mutable global
having appeared.

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

The 2026-07-18 revision re-ran the same `nm` sweep over the 107 linked object
files of the current `dev` build. Two useful results: (1) it confirmed that
every eliminated group is genuinely gone from `.data`/`.bss`, and (2) it caught
three mutable file-scope symbols the 2026-07-13 pass had not listed (the
showalign row buffers, `derep_sort_db`, and `unique.cpp`'s `hash_function` — see
"Surfaced by the re-sweep" below). All remaining mutable symbols now have
**internal** linkage (lowercase `b`/`d`): the Phase 2 external-linkage leaks
stay plugged and nothing new leaked an external symbol.


## Progress since the original inventory (34 → 4 of the original set)

| Group | File | Original count | Now | Landed |
|-------|------|---------------:|-----|--------|
| **A** ✅DONE | `core/db.cpp` | 12 | 0 | → owned `Database` class (`02bf3975` … `2ac3f27f`, `c16d80ef`, `782a16fd`); private `std::vector` via `FatalAllocator`, const read API. ABI break → api 0.10.0. |
| **B** ✅DONE | `core/dbhash.cpp` | 5 | 0 | → owned `Dbhash` class (`9d8d498a`, `2625497b`); threaded into `--search_exact`. Internal-only. |
| **C** ⚠ partial | `core/mergepairs.cpp` | 9 | 4 | 5 lookup tables → the `QualityTables` struct returned by `mergepairs_init`/`precompute_qual` and passed by `const &` (`d94aa9a7`, api 0.11.0). The 4 abort-state atomics remain — **deliberate** (see below). |
| **D** ✅DONE | `core/otutable.cpp` | 1 | 0 | The `otutable` singleton pointer → the owned `OtuTable` class (`0dc01686`; name buffers to `std::string` `a4d6b5b3`; stale include dropped `302ac812`). |
| **E** ✅DONE | `core/getseq.cpp` | 1 | 0 | `labels_data` → a local `std::vector<std::vector<char>>` in `getseq()`, passed by ref to `read_labels_file()` (fill) and `test_label_match()` (`const &`, match). File-local functions only; no header/ABI change. |
| **F** ✅DONE | `core/derep.cpp`, `core/kmerhash.cpp`, `commands/derep_smallmem.cpp` | 5 | 0 | The 3 `hash_function` pointers → `static constexpr` (Phase 1). `hashtable`/`hashtablesize` are now function-locals of `derep_smallmem` (`92e3dc7f` per-invocation struct, `2414ea0e` `std::vector`). |
| **G** ✅DONE | `utils/random.cpp` | 1 | 0 | `base_seed` + `random_init()` + `random_base_seed()` → the owned `RandomSeed` class (ctor resolves the seed from `Parameters`, `value()`/`substream()` read it; both `noexcept`, ctor not — `std::random_device` may throw). Constructed per consuming command (shuffle/fastx_subsample locals, sintax a `sintax_state_s` member read lock-free by workers); `random_init` call dropped from `vsearch.cc`. `random_substream_seed` moved to an anonymous namespace (internal linkage). Internal-only, no ABI change. |

Of the original 34, **30 are eliminated** and **4 remain** (only Group C's 4
abort-state atomics; Groups E and G were both eliminated 2026-07-21).


## Surfaced by the 2026-07-18 re-sweep (not in the original 34)

The fresh `nm` pass listed three mutable file-scope symbols the original sweep
had not captured. All have internal (anonymous-namespace / `static`) linkage, so
none is an ODR/link hazard; they are encapsulation candidates, not correctness
bugs.

- **`core/showalign.cpp`** — `q_line`, `a_line`, `d_line`, three
  anonymous-namespace `std::vector<char>` row buffers
  (showalign.cpp:78-80) resized and overwritten per alignment. They pre-date the
  original inventory (blame: 2025-08-21) but were simply not listed. Internal
  linkage; single reader (output is single-threaded). **[DONE, 2026-07-21]** —
  bundled into an `AlignmentRows` struct (members `query`/`symbols`/`target`)
  owned as a local by `align_show()` and threaded by reference through
  `putop()`/`putop_final()`/`print_alignment_block()`; the manual `clear()` calls
  are gone (RAII). **3 symbols → 0.**
- **`commands/derep_prefix.cpp`** — `derep_sort_db`, a file-scope
  `Database const *` set to `&db` purely so the C-style `std::sort` comparator
  could reach the database. Introduced by the Group A work itself (`9cef6395`,
  "give derep_prefix its own Database"): that migration removed the global `db`
  but reintroduced this narrower global bridge. It was the **only** such
  comparator-bridge global left (searchcore threads its `Database &` through
  `hit_compare_bysize_typed` as a parameter instead). **[DONE]** — eliminated in
  commit `090751fe` exactly as planned: the free-function comparator + global were
  replaced by a stateful `compare_prefix` lambda capturing `&db`
  (derep_prefix.cpp:313), passed to `std::sort` (line 342), done alongside the
  `std::qsort` → `std::sort` conversion. Verified absent from the source and from
  `.data`/`.bss` in the 2026-07-21 `nm` re-sweep. **1 symbol → 0.**
- **`core/unique.cpp`** — `hash_function` (unique.cpp:86), a non-`const`
  anonymous-namespace function pointer, never reassigned. The original Group F
  counted "×3" and missed this fourth one.  **[DONE]** — const-ified in place to
  `constexpr Hash hash_function = CityHash64;` (a function address is a valid
  constant-expression initializer; `static` omitted as redundant inside the
  anonymous namespace). Matches the Group F siblings' `.data.rel.ro`/RELRO
  placement. cppcheck clean; binary relinks. **1 symbol → 0.**

Adding these, the current mutable-global residue is **4**: Group C-abort (4) only.
(`derep_sort_db`, Group G's `base_seed`, Group E's `labels_data`, and the three
showalign buffers were all eliminated on 2026-07-21 — see their bullets/rows above.)


## Why Group C keeps its 4 abort-state symbols

`merge_abort`, `merge_error_claimed`, `merge_error_reason`, `merge_error_value`
(mergepairs.cpp:111-114) coordinate the merge worker threads and are exposed via
`request_merge_abort` / `merge_aborted` / `report_merge_abort`. `merge_abort` and
`merge_error_claimed` are `std::atomic`; the state is deliberately file-scope so
that a worker hitting an out-of-range quality value can signal the others without
racing a `std::exit()` mid-run. Any future `Merge_context` bundling these must
preserve the atomics and the join happens-before that currently guards that race;
this is lower-priority and higher-risk than E/G, so it stays last.


## Not targets (documented, left alone)

- `vsearch_api.cpp:session_mutex` — intentional process-wide `std::mutex`
  serialising the C-API `vsearch_session_begin`/`vsearch_session_end`. A
  legitimate singleton.
- `utils/logfile.cpp:log_handle` — the encapsulated log-file `std::FILE *`
  (anonymous namespace, reached via `handle()`/`set_handle()`). This is the
  deliberate end-state of the `fp_log`-global elimination: `fatal()`/`warn()`
  cannot take a `Parameters`, so the log sink is a single encapsulated
  process-wide handle behind accessors, analogous to `session_mutex`. Not a
  target.
- `fatal_detail::throw_on_fatal()::enabled`, `OtuTable::print_biomout(...)::date`,
  `fastx.cpp:blanks` — function-local statics (the last two are `const`); local
  scope, left alone.
- `default_quality_padding` / `alternative_quality_padding` /
  `default_sequence_padding` (`vsearch.hpp`) — `std::string const`, so already
  **const** (not mutable). Being `std::string const` in a header gives them
  internal linkage, i.e. one heap-allocating copy per TU (the many
  `default_*_padding` rows in the `nm` dump are these per-TU copies, not distinct
  globals). Converting to `constexpr char[]` is optional hygiene, deferred.
- `bz2_libname`/`gz_libname`, `userfields_names`/`help_message`,
  `cli.cc:option_specs` (`constexpr`), `attributes.cpp:attributes` (`constexpr`),
  `fastx_syncpairs.cpp:default_read_separators` (`const char * const`), the
  fasta/fastq char maps (`char_actions`, `chrmap_identity`,
  `char_fq_action_seq`, `char_fq_action_qual`, all `const std::vector`),
  `compare_strings_nocase.cpp:compare_chars` (a stateless captureless lambda) —
  all const or stateless.


## Resolved: the external-linkage leak (Phase 2)

The original inventory flagged that the Group C tables, `dbhash_table` (B) and
`labels_data` (E) had **external linkage with no `extern` declaration**, an
ODR/link hazard (`q2p` in particular coexisted with unrelated external `q2p`
*functions*). Phase 2 (`2ce96b81`) gave them all internal linkage, and the
2026-07-18 re-sweep confirms **no mutable global carries external linkage any
more** — every remaining mutable symbol is lowercase `b`/`d` in `nm`.


## Plan

### Phase 1 — const-ify the set-once (no ABI impact)  [DONE, commit dc22fb25]

- Group F: the three `hash_function` function pointers →
  `static constexpr Hash hash_function = ...` in `core/derep.cpp`,
  `core/kmerhash.cpp`, `commands/derep_smallmem.cpp`.
- `core/unique.cpp:hash_function` → `constexpr Hash hash_function = CityHash64;`
  (2026-07-18 follow-up, closing the "×3 vs ×4" gap).  [DONE]

### Phase 2 — plug the external-linkage leaks (low risk)  [DONE, commit 2ce96b81]

Gave internal linkage to the file-scope state that leaked external symbols. `nm`
confirmed each symbol moved from `B`/`D` (external) to `b`/`d` (local), and the
2026-07-18 re-sweep confirms it held.

### Phase 3 — encapsulate self-contained TU state into per-run structs (medium)

The proven `Dbindex`/`Database` RAII-struct playbook: thread an owned instance by
reference through the relevant functions.

- **D** `otutable` → owned `OtuTable` class.  [DONE, `0dc01686`]
- **C-tables** the 5 lookup tables → `QualityTables`, built by `mergepairs_init`
  and passed by `const &`.  [DONE, `d94aa9a7`]
- **F** `hashtable`/`hashtablesize` → folded into the derep_smallmem
  per-invocation struct.  [DONE, `92e3dc7f`, `2414ea0e`]
- **E** `labels_data` → thread through the `getseq` functions.  [DONE, 2026-07-21
  — local in `getseq()`, passed by ref to `read_labels_file`/`test_label_match`]
- **G** `base_seed` → the owned `RandomSeed` class (ctor from `Parameters`,
  `value()`/`substream()`), constructed per consuming command; `random_init()`/
  `random_base_seed()` removed, `random_substream_seed` made file-local.
  [DONE, 2026-07-21]
- **C-abort** bundle the 4 atomics into a `Merge_context` owned by `pair_all()`,
  preserving the atomics + join happens-before.  [pending, last — highest risk]
- **showalign** bundle `q_line`/`a_line`/`d_line` into a `ShowAlign` state or
  pass them as parameters.  [DONE, 2026-07-21 — `AlignmentRows` struct owned by
  `align_show()`, threaded by ref through `putop`/`putop_final`/`print_alignment_block`]
- **derep_prefix** replace the `derep_sort_db` global + free-function comparator
  with a stateful comparator capturing `db`.  [DONE, `090751fe` — `compare_prefix`
  lambda captures `&db`, passed to `std::sort`; landed with the `std::qsort` →
  `std::sort` conversion]

### Phase 4 — the DB core (large, ABI-breaking)  [DONE, merged to `dev`]

Done as two owned classes rather than one combined `Db`:
- **Group A** → an owned `Database` class (`datap`/`seqindex`/counts +
  `getX()` members), threaded through all 31 TUs + `udb.cpp` + the public API;
  later polished to private `std::vector` (via `FatalAllocator`) storage,
  `std::sort`, and const-correctness. ABI broken → api 0.10.0.
- **Group B** → an owned `Dbhash` class (`core/dbhash.cpp`), threaded into
  `--search_exact` (the sole consumer). Internal only, so no ABI change.

Both mirror the earlier `Dbindex` refactor. Full CLI suite + release
`api_examples` + Windows/POWER/mips64el cross-compiles all passed.

### Cross-cutting risks

- **Runtime-mutated-globals trap**: DB state is mutated after load by
  `udb_read` and the `Database::sortby*` members — Phase 4 kept a *mutable*
  `Database &` reaching those mutation points, not a `const` snapshot. This also
  applied to the (now-completed) `derep_sort_db` removal: `compare_prefix` reads
  `db` during the sort, so it captures `db` by reference and the sort runs before
  `db.clear()`, keeping the reference valid for the sort's lifetime.
- **Threading**: Group C's abort state is deliberately atomic file-scope state
  to avoid a `std::exit()`-during-worker data race; any `Merge_context` must
  preserve the atomics and the happens-before on join.
- Do the remaining Phase 3 items on a `tmp_` branch with the full test suite
  plus a release-build library-API run (debug libs segfault the `api_examples`
  net).
