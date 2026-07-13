# TBD — Factor the FASTA/FASTQ line-reading loops into a line iterator

Date: 2026-07-13
Status: **IMPLEMENTED** on branch `tmp_20260713082154` (commits `9a5b9c5b`..`5e352abd`),
awaiting human review. See "Implementation status" below.

## Implementation status (2026-07-13)

Done, in these commits:

1. `9a5b9c5b` — add the line primitive to `core/fastx.hpp`.
2. `c7a2ecc3` — convert the FASTA sequence loop (loop B).
3. `7093d73c` — convert the FASTA header loop (loop A); `fasta_next` fully done.
4. `e5cb287d` — convert all four FASTQ loops (C/D/E/F) + the two downstream
   `line_end` uses.
5. `5e352abd` — perf fix: split the primitive into fill + scan (see below).

Two deviations from the plan below, both deliberate:

- **Phases 4 and 5 were done in one commit.** All four FASTQ loops and the two
  downstream error-line adjustments (the `+`-line validity check and the final
  sequence/quality length check) shared the single `line_end` variable;
  splitting header/`+` from sequence/quality would have left fragile
  intermediate states. The `+`-line check is always reached with an
  LF-terminated plus line, so its number simplifies to `lineno - 1` (documented
  in code); the final length check keeps the conditional via the quality loop's
  `last_line_complete`, which mirrors the old `line_end != nullptr` at every
  exit path.

- **The primitive is `scan_line_fragment` (fill + scan split), not
  `peek_line_fragment`.** The first cut used a single `peek_line_fragment()`
  that refilled and scanned together. Benchmarking (release, `hyperfine`) showed
  a **~4 % regression on a FASTQ read** (`fastq_chars` over 2M single-line
  reads: 1.045 s vs 1.004 s), because `peek` ran a `std::memchr` on the
  sentinel-break iteration where the old code tested the boundary byte and broke
  *before* scanning; FASTA was flat (it amortises over multi-line records). The
  Risks section below anticipated exactly this ("split the peek into a fill-only
  step + a separate scan"). After the split — reader loops call the existing
  `fastx_file_fill_buffer()` (EOF via its `0` return), test the sentinel byte,
  then `scan_line_fragment()` — the regression is gone (1.008 s vs 1.006 s,
  within noise) and FASTA is unchanged. `Line_fragment` therefore has no
  `end_of_input` flag.

Validation performed:

- Full test suite `vsearch-tests/run_all_tests.sh`: **9698 PASS, 0 FAIL**.
- Byte-for-byte output parity vs the pre-refactor binary on plain/gzip/bzip2/
  long-line FASTA, `fastq_chars`/`fastx_revcomp` (including an input with no
  trailing newline), and every error path (identical messages *and* reported
  line numbers): illegal FASTA char, unterminated header, illegal FASTQ
  sequence char, short quality, bad `+` line, sequence/quality length mismatch.
- `cppcheck --force --enable=warning,style --std=c++11` clean on all three files
  (the one remaining `fastx.cpp:178` `wrongPrintfScanfArgNum` is the pre-existing
  `warn()` bug, present in the baseline, untouched here).
- Cross-compiles all green (full build): Windows `x86_64-w64-mingw32`, POWER
  `powerpc64le-linux-gnu`, `mips64el-linux-gnuabi64` (the CLAUDE.md "RISC-V"
  target), and ARM `aarch64-linux-gnu`. No intrinsics are involved (pure scalar
  `std::memchr`), so there was no per-ISA SIMD work.

The invariant checklist and the "free functions vs a class" / Layer-2 open
questions below still stand for the reviewer. The implementation used free
functions over `fastx_handle` (open question 1, the recommended option) and did
**not** add Layer-2 convenience readers (open question 2).

## Origin

This started from the question "could swarm's `NwAligner` and `Line_buffer`
classes be used in vsearch?".

- `NwAligner` (Needleman–Wunsch, *not* Newick): **no** — vsearch already has a
  SIMD 8-lane aligner (`core/align_simd.cpp`) plus a linear-memory Hirschberg
  scalar aligner (`core/linmemalign.cpp`, already a proper class). Swarm's
  scalar `O(n²)`-memory full-matrix aligner would be a regression. Not pursued
  here.
- `Line_buffer` (swarm, `src/utils/line_buffer.{h,cc}`): the **class itself is
  not reusable** (it wraps POSIX `getline()` over a plain `FILE*` and knows
  nothing about gzip/bzip2), but the *problem it names* — a reusable line
  abstraction instead of hand-rolled scan loops — is real in vsearch. This note
  designs a vsearch-native version.

## Verdict

Do **not** port swarm's `Line_buffer`. Instead, extract vsearch's own repeated
line-scanning body into a small primitive that sits **on top of** the existing
`fastx_file_fill_buffer` refill path, so gzip/bzip2 and file-position tracking
keep working unchanged. This is a readability/DRY refactor, not a behaviour
change.

## Current state — six near-identical loops

Line/record parsing today is a growable byte buffer (`fastx_buffer_s`) plus one
shared refill primitive (`fastx_file_fill_buffer`, `core/fastx.cpp:604`).
On top of that, the *same line-scanning loop is hand-written six times*:

| # | Location | Reads | Copy fn | Terminator | Record-boundary sentinel | EOF policy | `++lineno`? |
|---|----------|-------|---------|------------|--------------------------|------------|-------------|
| A | `fasta.cpp:288` | header | `buffer_extend` (raw) | `line_end == nullptr` (one line) | — | fatal "header must be terminated with newline" | yes, on LF |
| B | `fasta.cpp:325` | sequence | `buffer_extend` (raw) | `while(true)` | `>` | normal `break` | **no** (counted later in `fasta_filter_sequence`) |
| C | `fastq.cpp:366` | header | `buffer_extend` (raw) | one line | — | `fastq_fatal` "Unexpected end of file" | yes, on LF |
| D | `fastq.cpp:398` | sequence | `buffer_filter_extend` (filtered) | `while(true)` | `+` | `fastq_fatal` | yes, on LF |
| E | `fastq.cpp:466` | `+` line | `buffer_extend` (raw) | one line | — | `fastq_fatal` | yes, on LF |
| F | `fastq.cpp:527` | quality | `buffer_filter_extend` (filtered) | `while(true)` | `@` **and** `qual.length == seq.length` | normal `break` | yes, on LF |

Every one of the six contains this identical, error-prone core:

```cpp
auto * const current_position =
  std::next(input_handle->file_buffer.data,
            static_cast<std::ptrdiff_t>(input_handle->file_buffer.position));
line_end = static_cast<char *>(std::memchr(current_position, '\n', rest));
uint64_t len = rest;
if (line_end != nullptr)
  {
    len = static_cast<uint64_t>(line_end - current_position + 1);   // <-- pointer arithmetic + narrowing cast
  }
/* ... copy len bytes ... */
input_handle->file_buffer.position += len;
rest -= len;
```

That block is exactly the collection of anti-patterns `CLAUDE.md` calls out —
raw pointer arithmetic (`line_end - current_position + 1`), a narrowing
`ptrdiff_t → uint64_t` conversion papered over with `static_cast`, and a magic
`'\n'`. Writing it once and getting it right once is the main win.

### What genuinely varies (and must be preserved)

The six loops are *not* trivially identical — the abstraction has to leave these
decisions with the caller:

1. **Copy policy**: raw `buffer_extend` vs char-action-filtering
   `buffer_filter_extend` (which needs `char_action` + `char_mapping` tables and
   an `ok`/`illegal_char` out-pair).
2. **FASTA asymmetry**: FASTA sequence (loop B) is copied **raw** in the loop and
   filtered afterwards in a *second pass* (`fasta_filter_sequence`,
   `fasta.cpp:195`), while FASTQ filters inline. So the primitive must **not**
   bake in filtering.
3. **`lineno` accounting**: loop B does **not** bump `lineno` (newlines are
   counted later via `Action::count`); all others do. Easy to break.
4. **EOF policy**: fatal (A, C, D, E) vs normal loop exit (B, F). And the fatal
   helper differs: `fatal()`/deferred-error in FASTA vs `fastq_fatal()` in FASTQ.
5. **Record-boundary sentinel**: none / `>` / `+` / (`@` + equal-length) — checked
   *before* consuming, only once the previous line completed.
6. **Extra guards**: loop F breaks early if the quality run already exceeds the
   sequence length (`fastq.cpp:569`); loops D/F run an `ok` error check after each
   filtered copy.

## Why swarm's `Line_buffer` cannot be dropped in

- It is built on POSIX `getline()` over a `std::FILE *`. vsearch never uses
  `getline`; it refills its own buffer via `fread` / `gzread` / `BZ2_bzRead`
  (`fastx.cpp:630-672`). Swapping in `Line_buffer` would **lose gzip and bzip2
  support** (loaded at runtime through the `DynamicLibraries` dlopen wrapper) and
  bypass the `file_position` tracking + old-zlib `dup()`/`xlseek` workaround.
- It hands back whole raw lines; vsearch interleaves line-scanning with the
  char-action filter and multi-line record assembly. A pure "next line" API only
  covers part of the job.

We borrow swarm's *spirit* (encapsulate the buffer contract, expose read-only
views, RAII) but not its code.

## Proposed design

Keep the record-shape logic (what a header/sequence/plus/quality line *is*) in
`fasta_next`/`fastq_next`. Extract only the shared "find and consume the next
line fragment" mechanics into one small primitive living beside the existing
buffer helpers in `core/fastx.{cpp,hpp}`.

### Layer 1 — the core primitive (recommended, sufficient on its own)

A view of the next line fragment sitting in the file buffer, plus a consume
step. Both operate on the existing `fastx_handle`, matching the surrounding
C-style-over-a-struct idiom:

```cpp
// A read-only view of the next line fragment in the input buffer.
//   view          : bytes available now, up to and including the LF if present
//   has_newline   : true when `view` ends at an LF (the line is complete)
//   end_of_input  : true when the buffer could not be refilled (no more data)
// Reuses the existing Span<> view type (see find_header_end, parse_cigar_string).
struct Line_fragment
{
  Span<char const> view;
  bool has_newline = false;
  bool end_of_input = false;
};

// Refill the file buffer if empty and locate the next LF. Does NOT consume:
// the caller inspects the fragment (e.g. peeks view[0] for a record-boundary
// sentinel), copies it into whichever destination buffer it wants, then calls
// consume_fragment(). Marked inline in the header so cross-TU callers in
// fasta.cpp / fastq.cpp keep the current inlined hot path (no LTO assumed).
inline auto peek_line_fragment(fastx_handle input_handle) -> Line_fragment;

// Advance the file-buffer read position past a fragment just copied out.
inline auto consume_fragment(fastx_handle input_handle,
                             Line_fragment const & fragment) -> void;
```

`peek_line_fragment` folds the `fastx_file_fill_buffer` + `std::memchr` + length
computation into one place; `consume_fragment` folds the `position += len`
advance. All of the pointer arithmetic and the narrowing cast disappear from the
six call sites and live once, tested once.

The caller keeps full control of policy. Loop B (FASTA sequence) becomes:

```cpp
bool previous_line_complete = false;
while (true)
  {
    auto const fragment = peek_line_fragment(input_handle);
    if (fragment.end_of_input) { break; }                 // EOF: normal exit
    if (previous_line_complete and fragment.view[0] == '>') { break; }  // next record
    buffer_extend(&input_handle->sequence_buffer,          // raw copy, no filter
                  fragment.view.data(), fragment.view.size());
    consume_fragment(input_handle, fragment);
    previous_line_complete = fragment.has_newline;         // loop B: no ++lineno
  }
```

and loop A (FASTA header, one line, fatal on EOF, bumps lineno) becomes:

```cpp
while (true)
  {
    auto const fragment = peek_line_fragment(input_handle);
    if (fragment.end_of_input)
      { /* fatal / deferred-error "header must be terminated with newline" */ }
    buffer_extend(&input_handle->header_buffer,
                  fragment.view.data(), fragment.view.size());
    consume_fragment(input_handle, fragment);
    if (fragment.has_newline) { ++input_handle->lineno; break; }
  }
```

Filtered loops (D, F) call `buffer_filter_extend` on `fragment.view` instead,
then run their existing `ok`/`illegal_char` checks — unchanged.

### Layer 2 — optional convenience readers (only if Layer 1 leaves too much boilerplate)

If, after Layer 1, the three raw single-line loops (A, C, E) still look
copy-pasted, they could collapse into one helper — but the FASTA-vs-FASTQ
error-helper split (`fatal`+defer vs `fastq_fatal`) has to be resolved first
(see open questions). I recommend **doing Layer 1 first, measuring, and only
then** deciding whether Layer 2 earns its keep. Do not build Layer 2
speculatively.

## Invariants the refactor must preserve (review checklist)

- [ ] Loop B still does **not** touch `lineno`; loops A/C/D/E/F still bump it on LF.
- [ ] FASTA sequence is still copied raw and filtered in the second pass; only
      FASTQ filters inline.
- [ ] EOF is still fatal for A/C/D/E and a clean `break` for B/F, via the *same*
      error helper each uses today (`fatal`/deferred vs `fastq_fatal`) and the
      *same* message text and reported line number.
- [ ] Record-boundary sentinels fire only after a completed previous line
      (`previous_line_complete` mirrors the old `line_end != nullptr` guard).
- [ ] Loop F's "quality already longer than sequence" early break stays.
- [ ] Deferred-error mode (`defer_errors`, `fastx.hpp:128`) behaviour is
      byte-for-byte unchanged for worker threads.
- [ ] `file_position` tracking and the gzip `dup()`/`xlseek` path are untouched
      (they live below this layer in `fastx_file_fill_buffer`).

## Risks and considerations

- **Hot path / performance.** This is the input parser — run on every record of
  every input file. `peek_line_fragment`/`consume_fragment` must inline. Because
  `fasta.cpp` and `fastq.cpp` are separate TUs from `fastx.cpp` and the build
  does not assume LTO, the two helpers should be `inline` **in the header**
  (`fastx.hpp`), not defined in `fastx.cpp`. Verify with the hyperfine reader
  benchmarks per `CLAUDE.md` (`fastq_chars.sh`, and a plain-FASTA read path) on
  a release build before/after — target: no measurable regression.
- **One extra `memchr` per record.** `peek_line_fragment` runs `memchr` even on
  the iteration where a sentinel then breaks the loop (today the sentinel breaks
  *before* the `memchr`). That is at most ~2 extra `memchr` calls per FASTQ
  record over a partial buffer — expected to be lost in the noise, but it is a
  real difference; confirm in the benchmark. If it shows up, split the peek into
  a `fill`-only step + a separate scan.
- **No SIMD / cross-arch concern.** This code is pure scalar C++ (`std::memchr`);
  there are no intrinsics here, so no Power/ARM/x86 variants to keep in sync.
- **Scope creep.** Resist folding `buffer_filter_extend`, `fasta_filter_sequence`
  and the char-action tables into this change. They are a separate cleanup.

## Phased plan (small commits, per `CLAUDE.md`)

1. Add `Line_fragment` + `peek_line_fragment` + `consume_fragment` to
   `core/fastx.{hpp,cpp}` with no callers; `cppcheck` + build.
2. Convert loop B (FASTA sequence) — the simplest (raw copy, clean EOF, no
   lineno). Run conformity + reader benchmarks.
3. Convert loop A (FASTA header). Run tests.
4. Convert loops C, E (FASTQ header, plus — raw single-line). Run tests.
5. Convert loops D, F (FASTQ sequence, quality — filtered, sentinels, guards).
   Run tests. This is the highest-risk step; review the checklist above.
6. Full test suite (`vsearch-tests/run_all_tests.sh`), cross-compile smoke
   (Windows/POWER/RISC-V per `CLAUDE.md`), final `cppcheck` on all three files.

## Testing

- Correctness: `vsearch-tests/run_all_tests.sh` (fastx read paths are covered by
  many commands) + the conformity harness in `perf-results/hyperfine/conformity/`.
- Performance: `perf-results/hyperfine/fastq_chars.sh` and a plain-FASTA read
  command, release build, before vs after.
- Edge cases to eyeball: file with no trailing newline; a single record larger
  than `file_buffer.alloc` (forces mid-line refill); empty input; gzip and
  bzip2 inputs (exercise the refill path under the new primitive); CRLF line
  endings (the existing `'\r'` handling in the `+`-line check must still pass).

## Open questions for human review

1. **Free functions vs a thin class?** The proposal uses free functions over
   `fastx_handle` to match the existing `fastx_*` / `buffer_*` style. A
   `Line_reader` class wrapping the handle would echo swarm's `Line_buffer` more
   closely but adds no state (everything already lives in `fastx_s`). Preference?
2. **Layer 2**: worth collapsing the three raw single-line loops, given the
   `fatal` vs `fastq_fatal` error-helper split? Or leave them as thin Layer-1
   call sites?
3. Naming: `peek_line_fragment` / `consume_fragment` / `Line_fragment` — open to
   better names.
4. Is a preparatory step to unify `fatal`/deferred-error handling between the
   FASTA and FASTQ readers in scope, or strictly out of scope for this note?
