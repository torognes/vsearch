# TBD — Eliminate the globals in `dbindex.h` (2026-07-09)

Refactoring proposal, written on branch `tmp_20260709102359`. Part of the
Band 7 (architecture) "globals → owned state" effort, following the same
playbook as `userfields: replace globals with a Parameters member`
(`a8391422`) and, most closely, `dynlibs: replace globals with a RAII
DynamicLibraries class` (`4c6a8231`).

**Status:** proposal, awaiting human review before implementation.

---

## 1. Goal

`dbindex.h` exports **ten mutable globals** that together form a single
process-wide k-mer index singleton:

```
extern unsigned int  * kmercount;      /* matches per kmer                 */
extern uint64_t      * kmerhash;       /* per-kmer offset into kmerindex   */
extern unsigned int  * kmerindex;      /* the seqno lists                  */
extern struct bitmap_s * * kmerbitmap; /* per-kmer bitmaps (dense kmers)   */
extern unsigned int  * dbindex_map;    /* index-slot -> seqno              */
extern unsigned int    dbindex_count;  /* number of indexed sequences      */
extern unsigned int    kmerhashsize;   /* 4^wordlength                     */
extern uint64_t        kmerindexsize;  /* total entries in kmerindex       */
extern uhandle_s     * dbindex_uh;     /* unique-kmer finder (build only)  */
extern unsigned int    dbindex_wordlength; /* effective index width        */
```

(`bitmap_mincount` is an eleventh piece of shared state, a file-`static` in
`dbindex.cc`.)

The goal is to fold these into a single owned object so that no command relies
on process-global index state — matching the reentrancy/thread-safety aim of
E4 and the "no non-const global variables" rule in `CLAUDE.md`.

---

## 2. Why these globals are awkward today

- **They are one object pretending to be ten.** Every allocation in
  `dbindex_prepare` and every free in `dbindex_free` moves them as a group;
  their invariants (`kmerhash[i] + j < kmerindexsize`, `dbindex_count <=
  seqcount`, `kmerhashsize == 4^dbindex_wordlength`) are cross-field and
  entirely implicit.
- **`dbindex_wordlength` is derived index state, not config.** It is set by
  `dbindex_prepare` (from `parameters.opt_wordlength`) **and** overwritten by
  `udb_read` when a UDB file declares a different width (`udb.cc:407`). It must
  **not** be migrated to `Parameters` — that is exactly the trap recorded in
  the `e1_f3_runtime_mutated_globals` note (mutating `opt_wordlength` caused a
  live orient SIGSEGV). It belongs *inside* the index object.
- **The read API is a thread-safety contract in disguise.** The search engines
  read the index concurrently from worker threads (`searchcore.cc`
  `search_topscores`). That is safe only because the reads are const. Nothing
  in the current types says so.

---

## 3. Consumer map

Grouped by how they touch the index. This is the surface the refactor must
thread an object reference through.

| File | Role | Touches |
|------|------|---------|
| `dbindex.cc` | owner/impl | all globals + `bitmap_mincount` |
| `udb.cc` | owner + **direct builder** | fills the raw arrays from a UDB file (`udb_read`), reads them back when writing/reporting (`udb_make`, `udb_stats`, `udb_info`, `udb_fasta`) |
| `searchcore.cc` | **hot-path reader** | getters in `search_topscores`; `dbindex_wordlength` at `:849` |
| `sintax.cc` | owner + hot-path reader | lifecycle + getters (`:308-356`) + `dbindex_wordlength` (`:430`) |
| `orient.cc` | owner + reader | lifecycle + `dbindex_getmatchcount` (`:229-230`) + `dbindex_wordlength` (`:101,215`) |
| `cluster.cc` | owner (incremental) | `dbindex_prepare` (no `addall`), repeated `dbindex_addsequence`, `dbindex_wordlength` (`:597`) |
| `chimera.cc` | owner | `dbindex_prepare`/`addsequence`/`addallsequences`/`free` |
| `search.cc` | owner | `dbindex_prepare`+`addallsequences`+`free`, or `udb_read` |

**Key structural facts:**

- **Each command owns its own index** for the duration of one run
  (`prepare … use … free`). There is **no cross-command sharing**, so each
  command can own a local object and thread a reference down its own call
  tree — the refactor is per-command and independently committable.
- **`searchinfo_s` already carries `Parameters const *`** (`searchcore.h:162`,
  added in the E1 shared-infra phase). Adding a `Dbindex const *` beside it is
  the same, already-established pattern.
- **`udb_read` is the odd one out.** It builds the index by filling the raw
  arrays directly (bypassing `dbindex_prepare`) and is called by four commands
  (`orient`, `sintax`, `search`, `chimera`) plus `udb.cc` internals. It must
  gain a `Dbindex &` out-parameter to populate the caller's object.
- **The cluster *session* holds the index across calls** (`cluster.h:81-118`):
  `cluster_session_init` requires `dbindex_prepare` to have been called, and
  centroids are added incrementally via the session. The session must store a
  `Dbindex &` (like it already stores a `Parameters const &`).

---

## 4. The library-API constraint (important)

`dbindex.h` is included by the **public** `vsearch_api.h` (`:156`), and the
lifecycle functions are called **directly** by the shipped examples:

```
api_examples/example_search.cc, example_cluster.cc, example_chimera.cc,
example_lifecycle.cc, example_reinit.cc, example_dbinfo.cc
```

So changing these signatures is an **ABI/API break**. Per the
`e1_abi_decision` note, Band 7 E1 already chose to **break the
`libvsearch_core` ABI** (Option B), so this is acceptable — **but** the change
must be done in lockstep:

1. update `api_examples/*.cc`,
2. update `LIBRARY_API.md` (bump the documented API version),
3. re-run the library test net **against a release build** (the
   `library_tests_need_release_build` note: a debug `libvsearch.a` segfaults
   the examples via the `_GLIBCXX_DEBUG` ABI mismatch).

---

## 5. Target design

Introduce one type in `dbindex.h`:

```cpp
struct Dbindex
{
  /* owned buffers (public: this is a data-carrying struct, and udb_read fills
     them in place — see below). East-const, RAII-freed in the destructor. */
  unsigned int      * kmercount   = nullptr;
  uint64_t          * kmerhash    = nullptr;
  unsigned int      * kmerindex   = nullptr;
  bitmap_s * *        kmerbitmap  = nullptr;
  unsigned int      * dbindex_map = nullptr;
  uhandle_s         * uh          = nullptr;  /* build-time only */

  unsigned int   count        = 0;
  unsigned int   hashsize     = 0;
  uint64_t       indexsize    = 0;
  unsigned int   wordlength   = 0;   /* effective index width (derived state) */

  Dbindex() = default;
  ~Dbindex();                        /* RAII: frees everything (was dbindex_free) */
  Dbindex(Dbindex const &) = delete; /* owns raw buffers: non-copyable */
  auto operator=(Dbindex const &) -> Dbindex & = delete;

  /* lifecycle (was the free functions) */
  auto prepare(bool use_bitmap, int seqmask, Parameters const & parameters) -> void;
  auto add_sequence(unsigned int seqno, int seqmask) -> void;
  auto add_all_sequences(int seqmask, Parameters const & parameters) -> void;
  auto clear() -> void;              /* explicit free for reuse (was dbindex_free) */

  /* read API (const == thread-safe for concurrent search workers) */
  auto getbitmap(unsigned int kmer) const -> unsigned char *;
  auto getmatchcount(unsigned int kmer) const -> unsigned int;
  auto getmatchlist(unsigned int kmer) const -> unsigned int *;
  auto getmapping(unsigned int index) const -> unsigned int;
  auto getcount() const -> unsigned int;
};
```

Notes on the design choices:

- **RAII class, mirroring `dynlibs`'s `DynamicLibraries`.** The destructor
  subsumes `dbindex_free`; `clear()` stays for the in-place reuse case
  (`dbindex_prepare` currently calls `dbindex_free` first for idempotency, and
  `example_reinit` relies on prepare→free→prepare).
- **Public data members.** `udb_read` builds the index by writing the arrays
  directly; keeping them public lets that code become `dbindex.kmercount[i] =
  …` with no behavioural change and no `friend`. It also honours the project's
  "`struct` so all members are public" leaning.
- **`bitmap_mincount` becomes a local** in `prepare` (it is only read inside the
  same function's loop) — a small bonus cleanup; the file-`static` disappears.
- **`bitmap_threshold`** stays a `constexpr` in `dbindex.cc`.
- **The getters stay out-of-line** in `dbindex.cc` (defined as before, now as
  members) to keep codegen — and therefore the hot path — identical.

**API-shape decision to confirm (see §8, Q1).** The alternative to member
functions is to keep free functions that take an explicit `Dbindex &`
(`dbindex_prepare(Dbindex &, …)`, `dbindex_addsequence(Dbindex &, …)`). That is
a smaller diff for the library examples and closer to the current C-ish style,
but the RAII-class form is what `CLAUDE.md` ("prefer RAII") and the `dynlibs`
precedent point to. **Recommendation: the RAII class (Option B above).**

---

## 6. Incremental commit plan

Each step compiles and passes tests on its own. Commit messages follow the
`dbindex: …` convention; add the `Co-Authored-By: Florian FILLOUX` trailer.

1. **`dbindex: introduce the Dbindex struct alongside the globals`**
   Define `struct Dbindex` and its methods in `dbindex.{h,cc}`, implemented by
   *delegating to the existing globals* (thin wrappers). Nothing else changes
   yet. This lets every later step migrate one consumer at a time while the
   tree stays green. (Optional: skip if the big-bang per-file migration below is
   preferred — but the wrapper keeps commits small.)

2. **`dbindex: make bitmap_mincount a local in prepare`**
   Independent micro-cleanup; removes one file-`static`.

3. **`dbindex: thread Dbindex through udb_read`**
   Add `Dbindex & dbindex` out-param to `udb_read` and rewrite its raw-array
   fills as `dbindex.<member>`. Update all four external callers + the two
   internal ones to pass a local `Dbindex`. This is the biggest single step
   because `udb_read` is the direct builder — do it early so the UDB path is
   settled.

4. **`dbindex: migrate the search command`** (`search.cc` + `searchcore.cc`)
   Add `Dbindex const * dbindex` to `searchinfo_s`; owner in `search.cc`
   creates a local `Dbindex`, sets it on each per-thread `searchinfo_s`; the
   getters in `search_topscores` become `si->dbindex->getX(…)`; `:849` reads
   `si->dbindex->wordlength`.

5. **`dbindex: migrate the sintax command`** (`sintax.cc`)
   Same pattern (owner + hot-path reader in one file).

6. **`dbindex: migrate the orient command`** (`orient.cc`)

7. **`dbindex: migrate the chimera command`** (`chimera.cc`)

8. **`dbindex: migrate the cluster command + session`** (`cluster.cc`,
   `cluster.h`) Store `Dbindex &` in `cluster_session_s`; thread it to the
   incremental `add_sequence` sites and the `:597` `wordlength` read. Update the
   `cluster.h` doc comments (ask before touching comments — see `CLAUDE.md`).

9. **`dbindex: migrate the udb command family`** (`udb_make`, `udb_stats`,
   `udb_info`, `udb_fasta`) to read/write through a local `Dbindex`.

10. **`dbindex: delete the globals`**
    Remove the `extern` declarations and definitions; drop the thin wrappers
    from step 1 so the methods do the real work. Grep confirms zero remaining
    references. `dbindex_free` and the free-function getters are gone.

11. **`dbindex: update the library examples and LIBRARY_API.md`**
    Migrate `api_examples/*.cc` to the new API, bump the documented API
    version, refresh the numbered lifecycle comments in `vsearch_api.h`
    (`:74-84`, `:186`).

(If Option A / free-functions is chosen in §8, steps 4-10 pass `Dbindex &`
explicitly instead of via `searchinfo_s`/`cluster_session_s`, and step 11 is a
much smaller diff.)

---

## 7. Verification

- **Build matrix** (per `CLAUDE.md`): debug `--enable-debug`, then the three
  cross-compiles (mingw/Windows, POWER, RISC-V). `kmerhashsize`/`kmerindexsize`
  width and the `bitmap_s` pointer array are the portability-sensitive bits.
- **`cppcheck`** each modified `.cc`; **`clang-tidy`** for the new class
  (rule-of-five, const-correctness).
- **Test net:** `vsearch-tests/run_all_tests.sh`; run `scripts/orient.sh`
  manually (orient tests are disabled by default — `orient_tests_disabled`
  note) since orient is a migrated consumer.
- **Library net:** `api_examples` `make test` **against a release build**
  (`library_tests_need_release_build`).
- **Performance:** `search`/`cluster`/`sintax` share `search_topscores` — the
  hot loop. The getters must stay out-of-line and const so codegen is
  unchanged; confirm with `hyperfine` (`usearch_global.sh`, a `cluster` test)
  if the disassembly is not obviously identical.

---

## 8. Open questions for human review

- **Q1 — API shape.** RAII class with member functions (Option B,
  recommended) vs. free functions taking an explicit `Dbindex &` (Option A)?
  This decides the shape of the public library API and how steps 4-11 read.
- **Q2 — scope of this branch.** Land the whole sequence here, or split the
  library-example migration (step 11) into its own follow-up, as the `fp_log`
  work was split into per-area commits?
- **Q3 — comment edits.** Steps 8/11 touch explanatory comments in
  `cluster.h`, `dbindex.h`, and `vsearch_api.h` that describe the *global*
  lifecycle. `CLAUDE.md` requires asking before modifying comments — flagging
  now so the review can pre-approve the wording.
