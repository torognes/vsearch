# `Database` polish: COMPLETE

Status: **all items done and verified.** Branch context: `tmp_20260713082154`.
Date: 2026-07-13.

The DB-core refactor (Group A / Phase 4 of `TBD_20260713_globals_inventory.md`)
is complete: a `Database` struct in `src/core/db.hpp`/`db.cpp` wraps the former
`datap`/`seqindex` externs and the `db.cpp` file-statics; **every command and
library caller owns its own instance** and threads a reference through the code.
The transitional `db_global` and all the `db_*` free functions are gone.

The follow-up polish that this document originally deferred has now also landed
(commits `0b8646f7`..`b9f1150b`) and been verified:

- **Encapsulation** — `Database`'s data members are now **private**; access is
  through const getters + an `is_fastq()` accessor (member renamed
  `fastq_format`), with non-const `mutatesequence()`/`mutateheader()` for the
  in-place masking passes. `udb_read` is a `friend`.
- **`std::vector` storage** — `data_`/`seqindex_` are `std::vector` parameterised
  with a `FatalAllocator` (`utils/fatal_allocator.hpp`) that routes through
  `xmalloc`/`xfree`, preserving the raw buffers' fatal-on-OOM behaviour instead
  of throwing `std::bad_alloc` (which `-fno-exceptions` would turn into
  `std::terminate`). `add()` grows in `memchunk` steps via `reserve_in_chunks`.
- **`std::qsort` → `std::sort`** — the sort members use `std::sort`/`std::stable_sort`
  with capturing comparators; the transient `sort_datap`/`derep_sort_db` file
  globals are gone.
- **const-correctness** — getters return `char const *`; a read-only `View`
  (`utils/view.hpp`) carries the sequence spans into the aligner without a
  `const_cast`.
- Library API reconciled and bumped to **0.10.0**.

The `udb.cpp` blocker (below) was resolved by the `udb_reserve`/`udb_finalize`
seam + `friend` access, so the loader fills the private `std::vector` buffers in
place. Re-verification of this final state:

- **CLI test suite**: 0 failures (debug build).
- **cppcheck** (db.cpp/hpp, udb.cpp, fatal_allocator.hpp, view.hpp): only
  pre-existing findings; the old `db_add` `knownConditionTrueFalse` false
  positive is gone (the `std::vector` rewrite removed the flagged comparison).
- **Release `api_examples`** `make test`: 30/30 pass at api 0.10.0.
- **Cross-compiles**: Windows / POWER / mips64el all build clean.
- Tree left in the standard `--enable-debug -O0 -ggdb3` config.

The rest of this document is retained as the historical rationale for how the
`udb.cpp` coupling was broken.

## (Historical) The blocker: `udb.cpp` wrote the buffers directly

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


## (DONE) Encapsulation (private members) — how the udb blocker was broken

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


## (DONE) `std::vector` storage (RAII containers)

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


## Other items — all resolved

Every follow-up originally listed here has since been resolved:

- **`std::qsort` → `std::sort`** — done; the comparators capture the `Database`,
  and the transient `sort_datap` (db.cpp) / `derep_sort_db` (derep_prefix.cpp)
  file globals were removed. `std::stable_sort` is used where the old
  address-tie-break mattered, so duplicate ordering is preserved.
- **`db_add` grow-if-needed false positive** — gone: the `std::vector` rewrite of
  `add()` (insert + `reserve_in_chunks`) no longer has the `_old` comparison
  `cppcheck` flagged.
- **`db_is_fastq` asymmetry** — resolved: the member was renamed `fastq_format`
  and a proper `is_fastq()` accessor added.
- **Stale comments** — reworded to the `Database` method names; the
  `chimera.cpp` `#if 0` debug block uses `ci->db->getheader(...)`.

## Truly remaining (optional, future)

- **Group B `Dbhash`** (`core/dbhash.cpp`) is still process-global — the sibling
  of this refactor, not part of it. See `TBD_20260713_globals_inventory.md`.
- A full per-session context object (so more than one library session could be
  active at once) is still future work; some lower-level compute state remains
  process-global.
