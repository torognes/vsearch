# `Database` polish: deferred improvements

Status: the migration is DONE; the items below are follow-ups. Branch context: `tmp_20260713082154`.
Date: 2026-07-13.

The DB-core refactor (Group A / Phase 4 of `TBD_20260713_globals_inventory.md`)
is complete: a `Database` struct in `src/core/db.hpp`/`db.cpp` now wraps the
former `datap`/`seqindex` externs and the `db.cpp` file-statics, and **every
command and library caller owns its own instance** and threads a reference
through the code. The transitional `db_global` instance and all the `db_*`
free functions have been **deleted**; the original `datap`/`seqindex` externs
no longer exist. The library API is the `Database`-owned model (api version
0.9.0). Main build + full CLI test suite are green.

The design is intentionally conservative where it had to be: `Database` is a
`struct` with **public** data members and **raw `xrealloc` buffers**. This
document lists the improvements deliberately left as follow-ups, and explains
why the two most valuable ones are **blocked by `udb.cpp`**.

## Verification (all done)

- **CLI test suite** (`vsearch-tests/run_all_tests.sh`): 0 failures (debug build).
- **Release-build `api_examples`**: `configure CXXFLAGS=-O2` (no `_GLIBCXX_DEBUG`)
  + `cd api_examples && make test` → "All tests passed" (search/cluster/chimera
  batch-vs-sequential, dbinfo db.read/db.add/accessors/sort/indexing, reinit,
  merge, dust, derep — all match the `vsearch` binary).
- **Cross-compiles**: Windows (`x86_64-w64-mingw32` → `vsearch.exe`), POWER
  (`powerpc64le-linux-gnu`), and `mips64el-linux-gnuabi64` all build clean. (The
  change is pure reference-parameter threading with no CPU-intrinsic edits.)
- Tree restored to the standard `--enable-debug -O0 -ggdb3` config afterward.


## The blocker: `udb.cpp` writes the buffers directly

`udb_read()` (in `src/core/udb.cpp`) is a second database loader: instead of
building the DB record-by-record through `add()`, it reads a binary `.udb` file
and populates the buffers directly. In step 1 it does so through two local
reference-aliases bound to the global instance:

```cpp
char * & datap = db_global.datap;
seqinfo_t * & seqindex = db_global.seqindex;
```

Within `udb_read()` it then:

- `xmalloc`s `seqindex` as `seqcount * sizeof(seqinfo_t)` and writes individual
  fields (`seqindex[i].header_p`, `.headerlen`, `.seq_p`, `.seqlen`, `.qual_p`,
  `.size`) in indexed loops;
- `xmalloc`s `datap`, `largeread`s the header block straight into `datap` and
  the sequence block into `datap + udb_headerchars`;
- performs an **in-place `std::memmove` shuffle** over `datap`, rewriting each
  `seqindex[i].seq_p` offset and writing a terminating `\0` at computed
  positions.

This direct, low-level access to both buffers is what blocks the two
improvements below. Any of them requires first giving `udb_read()` a
**controlled population path** on `Database` (see "Prerequisite" below).


## Blocked improvement 1 — encapsulation (private members)

**Goal:** make `datap`, `seqindex`, the statistics, and the allocation
bookkeeping (`dataalloc`, `datalen`, `seqindex_alloc`) **private**, exposing
only the const read API and the build/sort mutators. This is the natural end
state (and what swarm's `Data` class already does).

**Why blocked:** `udb_read()` writes `datap`/`seqindex` (and per-record
`seqinfo_t` fields) directly. Private members would break it. The current
public-member `struct` mirrors the accepted `Dbindex` design, which has the
same "second loader" (`udb_read`) constraint and also kept public members.

**Unblock:** give `Database` a UDB population API and route `udb_read` through
it — e.g. a `friend`, or a builder such as
`Database::reserve_from_udb(seqcount, header_bytes, nt_bytes)` returning
writable spans plus a `finalize_udb()` that runs the memmove/terminator pass.
Once `udb_read` no longer touches raw members, they can become private.


## Blocked improvement 2 — `std::vector` storage (RAII containers)

**Goal:** replace the raw `xrealloc`/`memchunk` buffers with
`std::vector<char> data_` and `std::vector<seqinfo_t> seqindex_`, per the
project's "prefer RAII / prefer std" guidance and matching swarm's `Data`.

**Why blocked (two reasons):**

1. `udb_read()`'s raw buffer manipulation — `xmalloc`, `largeread` into
   `datap`/`datap + offset`, in-place `memmove`, manual `\0` terminators —
   would all have to be reworked against `std::vector::data()` and its size
   invariants. This is feasible but tightly coupled to the UDB binary format
   and must land together with the loader rework above.
2. Allocation model mismatch: the current buffers grow by fixed `memchunk`
   (2^24) via `xrealloc`, which calls `fatal()` on OOM. `std::vector` grows
   geometrically and throws `std::bad_alloc` — but vsearch is compiled
   `-fno-exceptions`. Moving to `std::vector` means accepting its growth policy
   and terminate-on-OOM behaviour (or a custom allocator that calls `fatal()`).
   Worth a short design note before committing.


## DONE: `db_global` eliminated (the migration)

The transitional `db_global` and the `db_*` forwarders have been removed; every
command owns a `Database` (a local, or a RAII member of its per-run state
struct), and the shared core threads a `Database const &` (read) / `Database &`
(load/sort) through the call graph. The encapsulation/`std::vector` work below
can now proceed on a clean `Database` with a single owner per command.

**How it was done (topology):** the read API is *not* consumed per-command in
isolation — the hot getters live in the **shared** search core (`chimera.cpp`,
`cluster.cpp`, `results.cpp`, `msa.cpp`, `searchcore.cpp`, `fasta.cpp`,
`dbhash.cpp`, `dbindex.cpp`); even `--sortbylength` reaches them through
`fasta_print_db()`. So the migration threaded `Database` through the shared core
first (adding a `Database const * db` to `searchinfo_s` beside `dbindex`, and a
`Database const & db` param elsewhere — mirroring the `Dbindex` refactor), then
flipped each owner to a local instance, then deleted `db_global`. The library
session-init functions (`search_session_init`, `cluster_session_init`,
`chimera_detect_init`/`_thread_init`/`_batch`, `search_batch`) each gained a
`Database const & db` parameter (the 0.9.0 ABI break).


## Other deferred items (not udb-blocked)

- **`std::qsort` → `std::sort`.** `std::qsort` was kept per instruction. Because
  a C comparator (a plain function pointer) cannot capture the `Database` and
  `qsort_r` is a non-portable linuxism, the comparators now reach the data
  through small **transient file-scope pointers** set immediately before each
  sort: `sort_datap` in `db.cpp` (the three `compare_*`) and `derep_sort_db` in
  `commands/derep_prefix.cpp`. These are single-threaded and internal, but they
  are non-const file-scope state — mild irony in a globals-elimination refactor,
  and the reason to revisit. Switching to `std::sort` with a comparator
  capturing the `Database` removes them and satisfies "prefer std algorithms" —
  **but** the comparators' final tie-break is on element *addresses*
  (`lhs < rhs`), whose result differs once elements move under `std::sort`. That
  only reorders true duplicates (equal length, size, and header), so output is
  equivalent in practice, but it is a behaviour change on duplicate-heavy input
  and must be verified against the full test suite before landing.

- **Stale comments (DONE).** The comments that named the deleted `db_*` free
  functions (in `db.cpp` `init()`/`add()`/`read()` + the `seqinfo_s` note,
  `fastx_subsample.cpp`, `fastx.cpp`, `search.hpp`, `dbindex.cpp`) were reworded
  to the `Database` method names, and the `chimera.cpp` `#if 0` debug block was
  updated to `ci->db->getheader(...)`.

- **`db_add` grow-if-needed false positive.** `cppcheck` flags
  `dataalloc > dataalloc_old` and `seqindex_alloc > seqindex_alloc_old` as
  "always false" (`knownConditionTrueFalse`) — a false positive: the preceding
  `while` loop can raise the value. The logic is correct and was copied verbatim
  from the original `db_add`; the warning only appeared once the variables
  became members. A behaviour-preserving cleanup that also silences it (and
  drops the `_old` locals):

  ```cpp
  if (needed > dataalloc)
    {
      while (dataalloc < needed) { dataalloc += memchunk; }
      datap = static_cast<char *>(xrealloc(datap, dataalloc));
    }
  ```

  Left out of step 1 to keep it a pure structural wrap; apply once reviewed.

- **`db_is_fastq` asymmetry.** All queries became `Database` const accessors
  (`getheader()`, `getsequencecount()`, …) except the FASTQ flag: an accessor
  named `is_fastq()` would clash with the `is_fastq` member, and renaming the
  member would falsify a preserved comment in `add()`. So `db_is_fastq()`
  forwards by reading the public member directly. When members go private,
  rename the member (e.g. `fastq_format`) and add an `is_fastq()` accessor.

- **Comment audit.** A comment in `Database::add()` still refers to "the global
  `is_fastq` flag". It stays accurate while `db_global` exists; once `db_global`
  is gone it should be reworded (requires human sign-off per `CLAUDE.md`).
