# Design doc — finish the RAII file-handling layer

**Status:** PROPOSED (branch `tmp_20260710155615`).
**Date:** 2026-07-10
**Origin:** the in-progress refactoring visible in `src/utils/open_file.{hpp,cpp}`.

---

## 1. Context: what exists today

Two new helpers in `src/utils/` are meant to replace the codebase's
discipline-based `FILE *` handling with RAII:

* `open_file.{hpp,cpp}`
  * `open_input_file(char const *) -> FileHandle`, where
    `FileHandle = unique_ptr<FILE, CloseFileHandle>` (plain `fclose`).
  * `open_output_file(char const *) -> OutputFileHandle`, where
    `OutputFileHandle = unique_ptr<FILE, CheckedCloseOutputHandle>`; the
    deleter runs `fflush` + `ferror` + `fclose` and `fatal()`s on a deferred
    write error rather than leaving a silently truncated file.
  * A filename of `nullptr` yields an empty handle; a filename of `"-"` is
    served by a `dup()` of stdin/stdout, with `EBADF`/`EMFILE` diagnostics
    centralized in `check_file_descriptor()`.
* `check_output_filehandle.{hpp,cpp}` — post-open validators
  `check_mandatory_output_handle` / `check_optional_output_handle` /
  `check_mandatory_fastq_output_handle`, taking `(filename, bool is_empty)`.

These replace **three** legacy mechanisms:

| Legacy | Location | Style |
|--------|----------|-------|
| `fopen_input` | `fastx.cc` | raw `FILE *`, manual `fclose` |
| `fopen_output` / `open_optional_output` / `fclose_output` | `util.cc`, `util.h` | raw `FILE *`, manual close |
| `xopen_read` / `xopen_write` | `system.h`, `os/*/system.cc` | fd-based (udb `mmap` path) |

### Migration state (measured 2026-07-10)

**Done or partial** (call `open_output_file`/`open_input_file`): `sortbysize`,
`sortbylength`, `rereplicate`, `shuffle`, `fasta2fastq`, `getseq` (input only),
`mask` (partial), `cut` (partial), `subsample` (partial), `fastx_syncpairs`
(partial — see §3.3).

**Still on the legacy output trio** (open-count = `fopen_output` +
`open_optional_output`; close-count = `fclose_output`):

| File | open | close | File | open | close |
|------|-----:|------:|------|-----:|------:|
| `cluster.cc` | 18 | 20 | `orient.cc` | 4 | 4 |
| `search.cc` | 16 | 17 | `fastqops.cc` | 3 | 3 |
| `search_exact.cc` | 15 | 16 | `eestats.cc` | 2 | 2 |
| `allpairs.cc` | 10 | 11 | `derep_prefix.cc` | 2 | 2 |
| `filter.cc` | 8 | 8 | `sintax.cc` | 1 | 1 |
| `fastq_mergepairs.cc` | 7 | 8 | `derep_smallmem.cc` | 1 | 1 |
| `chimera.cc` | 7 | 5 | `sff_convert.cc` | 1 | 1 |
| `derep.cc` | 5 | 4 | `dbindex.cc` | 1 | 1 |
| `subsample.cc` | 4 | 1 | `vsearch.cc` | 0 | 1 |
| `getseq.cc` | 4 | 4 | | | |
| `cut.cc` | 4 | 1 | | | |
| `mask.cc` | 2 | 2 | | | |
| `fastx_syncpairs.cc` | 2 | 1 | | | |
| `fastq_join.cc` | 2 | 1 | | | |

`udb.cc` (`fopen_output` ×1, `xopen_write` ×1) is **out of scope** — see §3.1.

## 2. Assessment: why finish

The approach is sound and worth completing. The wins:

* **RAII removes the manual-`fclose` bug class** — leaks and double-closes on
  early returns / `fatal()` paths / new branches. Today `fclose`/`fclose_output`
  is scattered across ~26 files by discipline.
* **Close semantics are encoded in the type.** The legacy rule "use
  `fclose_output` on outputs, plain `fclose` on inputs" (warned about in
  `util.h:186`) becomes the `OutputFileHandle` vs `FileHandle` distinction —
  impossible to get wrong.
* **Centralized `"-"` and errno handling** — legacy `fopen_input` dropped the
  `EBADF`/`EMFILE` distinction on `dup()` failure.

The one caveat: **today's half-migrated state is the worst outcome** — three
parallel mechanisms, no canonical one. The value is only banked when the legacy
trio is deleted (§7). This doc drives to that endpoint.

## 3. Scope decisions (agreed 2026-07-10)

### 3.1 udb is out of scope
`udb.cc` keeps `xopen_read`/`xopen_write` + `mmap` (needs the raw fd) and its
`fopen_output`. `xopen_read`/`xopen_write` stay in `system.h` with a one-line
comment on why they remain.

### 3.2 Test wording will be updated in sync
The new validators emit a generic `unable to open output file <name>`, dropping
the legacy per-option descriptions. A few tests grep the exact text:

* `cluster_fast.sh` → `Unable to open msaout file for writing`, `…consout…`,
  `…profile…`
* `fastx_getseqs.sh` → `Unable to open labels file`

These tests **will be edited in the same commit** as the code that changes the
message. This removes the need for a description-carrying validator variant —
the helpers stay as they are (a real simplification the decision buys us).

### 3.3 fastx_syncpairs local wrapper is replaced
`fastx_syncpairs.cc:141` defines a file-local `open_output_file(char *, char
const *)` that shadows the global name and calls `fopen_output` internally. It
will be removed and its call sites routed through the global
`open_output_file` + `check_*` helpers.

## 4. Phase 1 — improve `open_file.{hpp,cpp}` (start here)

Four changes, no behavior change except the `noexcept` markers:

1. **Collapse the two open paths into one raw-`FILE *` core** `open_stream()`.
   Today `open_output_file`'s `"-"` branch builds a `FileHandle`, `.release()`s
   it, and rewraps into an `OutputFileHandle` — because the internal
   `open_file`/`open_file_descriptor` helpers are hardcoded to `FileHandle`.
   A core that returns a raw `FILE *` removes the rewrap dance and the
   input/output duplication (nullptr-check + dash-branch + fopen).
2. **`is_dash()` predicate** — single `strcmp` site; the existing C++17
   `string_view` note attaches to it.
3. **Mark both deleters `noexcept`.** `std::fclose` cannot throw; `fatal()`
   terminates rather than throws (vsearch is exception-free), so the checked
   deleter is honestly `noexcept` too. Aligns with the CLAUDE.md noexcept
   guidance and lets the compiler enforce the `unique_ptr` destructor contract.
4. **Unify the empty-handle return** (input used a named `empty` local, output
   used `{nullptr}`) and **add header doc-comments** documenting the
   nullptr/`"-"`/failure contract, matching the `util.h` doc style.

Proposed `open_file.cpp` body (license header unchanged):

```cpp
#include "fatal.hpp"
#include "open_file.hpp"
#include <unistd.h>  // dup, STDIN_FILENO, STDOUT_FILENO
#include <cassert>
#include <cerrno>  // errno
#include <cstdio>  // std::fopen, fdopen
#include <cstring>  // std::strcmp


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  // C++17 refactoring:
  // constexpr std::string_view a_dash = "-";
  // simpler string comparisons: if (filename == a_dash) {

  // type-safe string mode wrapper
  struct ModeString {
    explicit constexpr ModeString(char const * str) noexcept
      : mode(str) {
    }
    char const * mode;
  };


  auto is_dash(char const * filename) -> bool {
    assert(filename != nullptr);
    return std::strcmp(filename, "-") == 0;
  }


  auto check_file_descriptor(int const file_descriptor) -> void {
    assert(file_descriptor >= -1);
    if (file_descriptor != -1) {
      return;
    }
    if (errno == EBADF) {
      fatal("original fd is not an open file descriptor.");
    }
    if (errno == EMFILE) {
      fatal("too many open file descriptors.");
    }
    fatal("cannot duplicate input or output stream.");
  }


  // Open a stdio stream and return a raw FILE * for the caller to wrap in the
  // matching RAII handle. A filename of '-' is served by a duplicate of the
  // given standard stream (STDIN_FILENO for input, STDOUT_FILENO for output).
  auto open_stream(char const * filename,
                   ModeString const & mode,
                   int const standard_fileno) -> std::FILE * {
    assert(filename != nullptr);
    assert(mode.mode != nullptr);
    if (is_dash(filename)) {
      auto const file_descriptor = dup(standard_fileno);
      check_file_descriptor(file_descriptor);
      return fdopen(file_descriptor, mode.mode);
    }
    return std::fopen(filename, mode.mode);
  }

}  // end of anonymous namespace


// read_file, file to read, open_input_file, open_istream
auto open_input_file(char const * filename) -> FileHandle {
  if (filename == nullptr) {
    return FileHandle{nullptr};
  }
  auto const mode = ModeString{"rb"};  // r: reading, b: non-UNIX environments
  /* open the input stream given by filename, but if name is '-' then
     use a duplicate of stdin (fd = STDIN_FILENO = 0) */
  return FileHandle{open_stream(filename, mode, STDIN_FILENO)};
}


auto CheckedCloseOutputHandle::operator()(std::FILE * file_handle) noexcept -> void {
  if (file_handle == nullptr) {
    return;
  }
  /* A write error (full disk, quota, broken pipe) is often deferred by stdio
     until the buffer is flushed, so check fflush and the error flag before
     closing; fclose also flushes and can report the same error. Fail loudly
     rather than leave a silently truncated output file. */
  if ((std::fflush(file_handle) != 0) or (std::ferror(file_handle) != 0)) {
    fatal("Unable to write to output file (disk full, quota exceeded, or broken pipe?)");
  }
  if (std::fclose(file_handle) != 0) {
    fatal("Unable to close output file (disk full or quota exceeded?)");
  }
}


// write_file, file to write, open_output_file, open_ostream
auto open_output_file(char const * filename) -> OutputFileHandle {
  if (filename == nullptr) {
    return OutputFileHandle{nullptr};
  }
  auto const mode = ModeString{"wb"};  // w: writing, b: binary (no \n->\r\n on non-UNIX), matches input "rb"
  /* open the output stream given by filename, but if name is '-' then
     use a duplicate of stdout (fd = STDOUT_FILENO = 1) */
  return OutputFileHandle{open_stream(filename, mode, STDOUT_FILENO)};
}
```

Proposed `open_file.hpp` changes (deleter signatures + doc-comments):

```cpp
struct CloseFileHandle {
  auto operator()(std::FILE * file_handle) noexcept -> void {
    static_cast<void>(std::fclose(file_handle));
  }
};

// ...

struct CheckedCloseOutputHandle {
  auto operator()(std::FILE * file_handle) noexcept -> void;  // defined in open_file.cpp
};

// ...

/* Open a named input stream in binary read mode. A null filename yields an
   empty handle; "-" reads a duplicate of stdin. On failure the handle is empty
   (the caller checks it). */
auto open_input_file(char const * filename) -> FileHandle;

/* Open a named output stream in binary write mode. A null filename yields an
   empty handle; "-" writes a duplicate of stdout. On failure the handle is
   empty; validate with the check_*_output_handle helpers. */
auto open_output_file(char const * filename) -> OutputFileHandle;
```

Note on the retained `nullptr` guard in `CheckedCloseOutputHandle`: a
`unique_ptr` never invokes its deleter on a null pointer, so the guard is
defensive (unreachable via the handle) — kept because the deleter is a
named type a maintainer might call directly.

**Phase 1 verification:** debug + release build clean (no new warnings);
`cppcheck --force --enable=warning,style` on `open_file.cpp` clean; the seven
already-migrated commands (`sortbysize`, `sortbylength`, `rereplicate`,
`shuffle`, `fasta2fastq`, `mask`, `cut`) produce byte-identical output to the
pre-change binary, exercising both the named-file and `"-"` (stdin/stdout)
paths; mingw cross-compile clean (the `"wb"` binary mode is the reason this
layer exists).

## 5. Phase 2 — merge open + check into the openers (API pivot, 2026-07-10)

Superseding the original two-step design (`open_output_file` +
`check_*_output_handle`): an old maintainer note proposed folding filehandle
creation and validation into one call parameterized by the option name. Adopted
after review — it deletes the duplicated mandatory validators (which differed
only by the hardcoded `--output`/`--fastqout` string), removes the awkward
`(not handle)` bool argument, and fixes the real wart hit in the first Phase 3
commits: commands whose option is neither `--output` nor `--fastqout` (sintax
`--tabbedout`, derep_smallmem `--fastaout`) had no fitting mandatory validator.

**Decisions (all recommended options taken):**
* **Timing** — pivot now, before the bulk of Phase 3, so the ~20 remaining files
  are migrated once and the ~10 already-migrated files get one mechanical update.
* **Messages** — option-named and standardized. Null: `output file must be
  specified with <option>`. Open failure: `unable to open output file for
  <option> (<filename>)`. Replaces the bespoke null messages (sintax,
  derep_smallmem, sff_convert) and updates the few wording-asserting tests
  (`cluster_fast`, `fastx_getseqs`) in the same commits.
* **Option type** — wrapped in `struct OutputOption { char const * name; }` so it
  cannot be transposed with the `char const *` filename (mirrors `ModeString`,
  satisfies CLAUDE.md's "avoid easily swappable parameters").
* **`check_output_filehandle.{hpp,cpp}`** — deleted once absorbed.

**New API (in `open_file.{hpp,cpp}`, added 2026-07-10):**
```cpp
struct OutputOption { explicit constexpr OutputOption(char const *) noexcept; char const * name; };
auto open_mandatory_output_file(char const * filename, OutputOption option) -> OutputFileHandle;
auto open_optional_output_file (char const * filename, OutputOption option) -> OutputFileHandle;
```
`open_output_file` stays as the low-level primitive both wrappers call; the
single-arg `fatal` is used for the two-substitution failure message (built as a
`std::string`, since `fatal` is not variadic — only `fatal(msg)` and
`fatal(fmt, one-string)` exist), which is also `%`-safe.

**Sequencing:** (C1) add the openers ✅ → (C2) retrofit already-migrated
mandatory callers (sortbysize, sortbylength, rereplicate, shuffle, mask,
fasta2fastq) → (C3) retrofit Phase-3 mandatory callers (sintax, derep_smallmem,
sff_convert; drop their bespoke guards) → (C4) retrofit eestats (keeps its
explicit pre-input null-guard for fail-fast ordering, uses
`open_optional_output_file` for the open) → (C5) delete `check_output_filehandle`
+ update `Makefile.am`. Then resume Phase 3 on the new API.

## 6. Phase 3 — migrate output streams (one command per commit)

Per file: `fopen_output`/`open_optional_output` → `open_mandatory_output_file` /
`open_optional_output_file` (§5); store handles as `OutputFileHandle` (in the
command's file structs where they outlive a single scope); delete the matching
`fclose_output`; update any test asserting the old wording **in the same
commit** (§3.2).

Order — smallest/lowest-risk first, wording-sensitive `cluster` last:

1. **Single-stream:** `sintax`, `derep_smallmem`, `sff_convert`, `dbindex`
   (`kmercounts.txt`), `eestats`, `derep_prefix`.
2. **Finish the partials:** `mask` (`fastx_mask`), `cut`, `subsample`,
   `fastx_syncpairs` (+ remove the local wrapper, §3.3), `fastq_join`.
3. **Multi-stream:** `filter`, `fastqops`, `derep`, `orient`, `getseq`,
   `vsearch.cc`.
4. **Large:** `fastq_mergepairs`, `chimera`, `allpairs`, `search_exact`,
   `search`.
5. **Last (wording-sensitive tests):** `cluster` (`msaout`/`consout`/`profile`).

**Watch during Phase 3:**
* Structs that hold handles become move-only (`unique_ptr` member) — audit for
  accidental copies; add defaulted moves where a struct is returned/reseated.
* Handles handed to worker threads must outlive the `ThreadRunner` scope.
* `exit()`-via-`fatal()` does **not** run local destructors, so an output
  handle held as a function local is *not* flushed on an unrelated fatal-error
  path (OS closes the fd, stdio buffer is lost). This matches current behavior
  and the legacy path; note it but do not try to fix it here.

### 6a. Execution progress (branch `tmp_20260710155615`)

* **Single-stream — DONE.** sintax, eestats, derep_smallmem, sff_convert,
  derep_prefix, dbindex (the last is `#if 0` debug code, compile-checked by
  temporarily enabling it).
* **Partials — DONE.** cut, subsample, fastx_syncpairs (local `open_output_file`
  shadow renamed to `open_output`, §3.3 resolved), fastq_join. Struct-held
  handles became `OutputFileHandle` members; each command's `close_*` helper was
  kept but converted to reset the handles **in the original fixed order** (RAII
  scope-exit runs in reverse, which would flip flush order for streams sharing
  stdout); the structs are now move-only.
* **Multi-stream — DONE.** filter, fastqops, derep, orient, getseq. `vsearch.cc`
  needed **no change** — its only `fclose_output` mention is an explanatory
  comment in `main()`; it opens no output files. Raw-`std::FILE *`-local files
  used one of two shapes: (a) *direct* — make `fp_*` an `OutputFileHandle` and
  `.get()` at the write sites (filter, since its writes guard on `opt_* !=
  nullptr`, avoiding 16 alias variables); (b) *alias* — keep the raw `fp_*` and
  add owning handles beside them, used where opens are conditional/separated and
  write guards are best left untouched (orient, derep, fastqops). Both preserve
  close order via `handle.reset()` at each original `fclose_output` point.
* **Verification.** Every command A/B byte-identical vs the pre-pivot baseline
  (named + `-` stdout, and both-to-stdout ordering where two outputs can share
  stdout); full `vsearch-tests` green (9648/0) after the pivot and after the
  partials; `orient.sh` (disabled by default) 150/0.
* **Remaining:** the **large** group (fastq_mergepairs, chimera, allpairs,
  search_exact, search) and **cluster** (wording-sensitive tests), then Phase 4.

## 7. Phase 4 — input streams, then delete the legacy trio

* Route local-lifetime input `fopen` through `open_input_file` (`getseq` is the
  template). For `fastx_open`, the `fp` lives inside `fastx_s` alongside gz/bz
  handles closed in `fastx_close`; at minimum adopt `open_input_file` for the
  `"-"`/`dup`/errno handling, or leave it as a documented boundary if converting
  the ownership proves invasive.
* Once no callers remain, **delete** `fopen_output` / `open_optional_output` /
  `fclose_output` from `util.cc`/`util.h` and `fopen_input` from `fastx.cc`.
  This is the commit that banks the value. Also reword the `vsearch.cc` `main()`
  comment that currently says stdout "is not closed through `fclose_output`" —
  after deletion the checked close lives in `CheckedCloseOutputHandle`.

## 8. Phase 5 — final verification

`cppcheck` every modified file; debug build; **mingw cross-compile**; run
`vsearch-tests` and grep `FAIL` (expect the updated-wording tests to pass). No
perf runs — open/close is not a hot loop.

## 9. Risks & open questions

* **Q (input-side error messages):** input open failures are currently reported
  ad hoc by each caller (`fastx.cc` → `Unable to open file for reading (%s)`;
  `getseq` → its own labels message). Keep them caller-side, or add a symmetric
  `check_*_input_handle` for consistency? Proposed: keep caller-side for now
  (input surface is small); revisit if it turns out noisy.
* **Move-only structs (Phase 3):** the main mechanical risk; caught at compile
  time, but may cascade into signatures of the per-command file structs.
* **`fastx_open` ownership (Phase 4):** the only genuinely invasive site; the
  fallback (reuse `open_input_file` for `"-"`/errno only) is cheap and keeps the
  win without restructuring `fastx_s`.
