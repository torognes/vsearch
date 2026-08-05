# Plan: replace the C-style output calls with typed writers (2026-08-04)

## Goal

Remove the run-time format string from vsearch's output paths, as swarm
did between 2026-08-01 and 2026-08-04. Every write becomes a call whose
argument types are checked by the compiler instead of matched against a
format string by convention:

```c++
fprint(handle, header)         // View<char>          was "%.*s" / fputs
fprint(handle, '\t')           // one character       was fputc
fprint(handle, "\t*\n")        // a literal           was fputs
fprint_integer(handle, count)  // one number, decimal was "%" PRIu64 / "%d"
```

The direct consequences, in order of value:

1. **`PRIu64` / `PRId64` and `<cinttypes>` disappear from the tree.**
   282 macro uses over 164 call sites in 32 files, plus 34 `#include
   <cinttypes>` lines. These macros exist only because a format string
   has to name the width of its argument; with no format string there is
   nothing to name. They are also the last reason vsearch has to care
   that MinGW's msvcrt `printf` does not understand `%llu`.
2. **`-Wformat=2` becomes possible.** `src/Makefile.am:28` records why it
   is off: `-Wformat-nonliteral` fires on `fatal()`, which passes a
   run-time format string to `std::fprintf` (`utils/fatal.cpp:86`,
   `:106`, 26 callers). That is the same class of defect the migration
   removes, so the blocker goes with it.
3. **A measured win in the record writers.** swarm's finding, which
   applies here: fortification stops GCC from folding `fprintf` into
   anything cheaper (while it *does* still fold `fputs` of a literal into
   an `fwrite`), so each call parses its format at run time, once per
   record. The local toolchain predefines `_FORTIFY_SOURCE 3` at `-O2`
   (`gcc -O2 -dM -E -x c /dev/null`), so release builds here are in that
   regime; confirm in the disassembly for vsearch's own flags before
   citing it as the reason for a measurement. The View/Span propagation
   pass measured 1.02–1.18× on this effect alone when it dropped `"%.*s"`
   (`TBD_20260729`, landed).
4. **No varargs, no narrowing to the `int` a `%.*s` precision needs, no
   silent `int64_t`→`uint64_t` casts inserted only to satisfy a macro.**

vsearch is *not* going to stop calling `std::fprintf`: 167 call sites
format a `double`, and formatting a `double` byte-identically is the one
thing nothing simpler does (see Decision 1). The target is the integer,
string, character and literal writes — 285 of the 476 sites — plus the
integer fields *inside* the 167 float-bearing lines, which is what
actually retires the `PRI` macros.


## The document this supersedes

`TBD_20260725_replace_C_functions.md:5` states: *"the **printf family is
excluded** throughout (`fprintf`/`snprintf`/… are the project's chosen
output mechanism and are not refactoring targets)"*. That exclusion is
**reversed** (Decision 0, 2026-08-04): the printf family is a target, in
full scope. Amend that line to point here in Phase 1 so the two documents
do not contradict each other while this runs.

Two inventory corrections to that document, while it is open: its
`<cstdio>` table (`:149`) lists 3 `std::fputc` and 1 `std::fputs` sites;
the measured counts are 115 and 124. Its statement that write errors are
already checked (`open_file.cpp:158` does `fflush`/`ferror` before
`fclose`) is correct and still holds — vsearch does not need swarm's
`065a2e0` fix.


## What swarm did, in order

| commit | what |
|---|---|
| `71f6ead` | `utils/decimal_digits.hpp`: the digit loop, sink-agnostic. No call site converted. |
| `e9f5595` | `fputc`/`fputs` wrapped as `fprint` overloads. No call site converted. |
| `e604dad` | first converted family: the five internal-structure columns |
| `e97ec88` | abundance annotations |
| `164e592` | the `-s`, `-u`, `-r` output files |
| `11fab73` | the log/report block; `<cinttypes>` leaves the tree |

Two primitives, then one family per commit, cheapest first, each verified
byte-identical before the next. That ordering is worth copying exactly:
the primitives land with *no* call site converted, so their shape can be
reviewed before ~300 mechanical edits depend on it.


## Inventory (measured 2026-08-04, on `dev`, comments stripped)

Over `src/**/*.{cc,cpp,hpp,h}`, C comments and `//` comments removed
before counting, so prose mentioning `%s` does not inflate the numbers.

| | sites |
|---|---|
| `std::fprintf` | 436 |
| `std::printf` | 40 |
| `std::fputs` | 124 (116 take a literal, 8 a run-time pointer) |
| `std::fputc` | 115 (**all** take a character literal) |
| `std::snprintf` | 30 |
| `fprint(handle, View)` already in use | 35, in 27 files |

The 476 `printf`-family sites classified by what they format:

| bucket | sites | disposition |
|---|---|---|
| plain convertible (int / string / char, no width) | 227 | convert fully |
| literal only | 35 | convert fully |
| constant or dynamic field width, no float | 47 | needs a width writer (Decision 2) |
| formats at least one `double` | 167 | convert every non-float field, keep one `fprintf` |

`PRI` macros: 282 uses, 164 sites, 32 files, 34 `<cinttypes>` includes.
41 of those 164 sites also format a float — they still lose their `PRI`
macros, because only the `double` field keeps a format string. Two of the
34 includes are already dead today (`cli.cc`, `core/searchcore.cpp`: the
include and its comment, no `PRI` use) and can go in the first commit of
any phase, or on their own.

`double` conversions: 237, dominated by `%.1f` (80), `%.2f` (40), `%.0f`
(17), then fixed-width table forms (`%5.1f` 12, `%10.2f` 9, `%10.1f` 6,
`%8.2f` 6, `%9.1f` 5, `%6.2f` 4) and the `fastq_eestats` precision ladder
(`%.5f` … `%.13f`, 3 each).

Everything else in the format-specifier census is one or two sites:
`%02x` (`core/fasta.cpp:295`) and `%lx` (`commands/version.cpp:92`, zlib
compile flags); no site formats a `%p`.

Reproduce with:

```sh
cd ~/Documents/src/vsearch/src
grep -rn 'PRIu64\|PRId64' --include=*.cpp --include=*.cc --include=*.hpp --include=*.h . | wc -l
grep -rlc 'include <cinttypes>' --include=*.cpp --include=*.cc --include=*.hpp . | wc -l
```

The bucket classification came from a throw-away Python pass that
bracket-matches each call, concatenates its string literals (rendering
`PRI*` as `d`) and reads the conversions; it is not worth keeping in the
tree, but the per-file table it produced is in the phase list below.


### Per-file, largest first (total / plain / padded / float / literal-only)

```
core/results.cpp                 49  32   2  15   0     hot: every hit, every format
commands/fastq_mergepairs.cpp    40   1  18  21   0     cold: the stats block
core/chimera.cpp                 39  11   7  19   2     warm: per query + alignments
commands/udbstats.cpp            30   8   9  13   0     cold: one table
core/align_simd.cpp              24   2   5   0  17     debug dumps, callers under #if 0
commands/orient.cpp              19  11   0   8   0     hot: per sequence
core/fasta.cpp                   19   7   1  11   0     hot: per sequence
commands/fastq_stats.cpp         14   6   0   8   0     cold
commands/udbinfo.cpp             14  12   0   2   0     cold
core/fastq.cpp                   14   4   0  10   0     hot: per sequence
core/otutable.cpp                13  12   0   0   1     hot: per cell (OTUs x samples)
commands/derep_smallmem.cpp      12   6   0   6   0     hot: per unique
core/cluster.cpp                 12   6   0   5   1     warm: per cluster
core/derep.cpp                   12   8   0   4   0     hot: per unique
commands/version.cpp             11   1   0   0  10     cold
```

then a tail of 34 files with 10 or fewer each (`search_exact`,
`sff_convert`, `usearch_global`, `db`, `fastq_eestats2`, `sintax`,
`logfile`, `derep_prefix`, `fastq_chars`, `fatal`, `progress.hpp`,
`allpairs_global`, `fastx_mask`, `fastx`, `cli.cc`, `udb`,
`fastx_subsample`, `dbindex`, `getseq`, `showalign`, `vsearch.cc`,
`fastq_convert`, `rereplicate`, `sortbylength`, `sortbysize`, `filter`,
`msa`, `cigar`, `cut`, `fastq_eestats`, `fastq_join`, `fastx_syncpairs`,
`help`).


## The primitives to add

### `src/utils/decimal_digits.hpp` (new)

Port swarm's file, which is already C++11 and already carries the
reasoning; adapt to vsearch conventions (`#pragma once` rather than an
include guard, vsearch's licence header, and the *global* namespace to
match today's `fprint(View)` — Decision 4; the `decimal` namespace around
`to_decimal`/`Buffer` stays, it is swarm's and it earns its keep). It
fills a caller-supplied `std::array<char, 20>` and returns a
`View<char>` over the part it used, knowing nothing about the sink:

```c++
namespace decimal {
  using Buffer = std::array<char, max_width>;     // 20: UINT64_MAX and INT64_MIN both fit
  template <typename Integer>
  auto to_decimal(Buffer & buffer, Integer value) noexcept -> View<char>;
}
```

Sink-agnostic on purpose: `utils/cigar.cpp` and `core/linmemalign.cpp`
append digits to a `std::string`/buffer and must not pull in `<cstdio>`.

This is a C++11 stand-in for `std::to_chars` (C++17). The three details
swarm's version already gets right, all of which apply verbatim here:
`INT64_MIN` (take the magnitude as `0 - static_cast<uint64_t>(value)`, so
the conversion precedes the negation), unsigned arguments (tag-dispatch
on signedness, because `if (value < 0)` is always false for an unsigned
type and `-Wtype-limits` is on in vsearch's debug build), and `char`
(integral, so `to_decimal('A')` would silently print 65 — a
`static_assert` rejects it).

Not `std::to_string`: before libstdc++ 11 its integral overloads are
`__gnu_cxx::__to_xstring(&std::vsnprintf, "%ld", …)`, i.e. the run-time
format string this exercise removes, on the oldest compiler vsearch
supports (GCC 4.9). That also makes the 43 existing `std::to_string`
calls a (separate) part of the same problem — see Out of scope.

### `src/utils/print_view.hpp` (extend)

Already holds `fprint(std::FILE *, View<char>)` and its rationale. Add:

```c++
inline auto fprint(std::FILE *, char) -> void;                        // was fputc
template <std::size_t Size>
auto fprint(std::FILE *, char const (&literal)[Size]) -> void;        // was fputs
template <typename Integer>
auto fprint_integer(std::FILE *, Integer) -> void;                    // was "%d" / "%" PRIu64
template <typename Integer>
auto fprint_integer(std::FILE *, Integer, std::size_t width) -> void; // was "%10u" / "%*d"
inline auto fprint_spaces(std::FILE *, std::size_t count) -> void;    // was "%*s" with ""
```

The literal overload writes with `std::fwrite` and the array's own bound
rather than calling `fputs`, so the length is a compile-time constant by
construction rather than by optimisation; its contract is therefore the
array bound, not a terminator (`static_assert(Size > 0)` plus an
`assert` on the terminator). Unlike swarm, **vsearch declares
partially-filled `char` arrays** — `core/cluster.hpp:70` `char
centroid_label[1024]`, `:72` `char cigar[4096]`, and
`core/chimera.hpp:74`–`:77`'s four `char [1024]` label buffers — so the
"the contract is the bound, not the terminator" gap is reachable here,
not theoretical. Passing `centroid_label` to this overload would emit all
1024 bytes. Those buffers go through `View` (with the length the writer
already knows) and never through the literal overload; the `assert`
catches the unterminated case but not the filled-then-truncated one,
which is why this has to be stated at the function and in the commit.

The last two are what swarm did not need (Decision 2). `fprint_integer`
with a width pads on the left with spaces to match `%10u`; it must not
try to reproduce `%-10u` (left-aligned; no site needs it) or `%02x`
(zero-padded hex; two sites, which keep `fprintf`).


## Decisions taken (maintainer, 2026-08-04)

All six are settled; the rationale each one was decided on is kept below,
because the phases and the verification protocol depend on it.

| | decision |
|---|---|
| 0. Scope | **Reverse `TBD_20260725:5`, full scope.** The printf family is a refactoring target; amend that line to point here. |
| 1. `double` fields | **One `fprintf` per line, for its `double` field(s); convert everything around it.** swarm's rule, applied in the hot writers too. |
| 2. Field width | **Add `fprint_integer(handle, value, width)` and `fprint_spaces()`.** Covers all 47 padded sites, run-time widths included. |
| 3. `align_simd.cpp` | **Out of scope**, all 24 debug-dump sites. |
| 4. Namespace | **Global now**, matching today's `fprint(View)`; all four overloads move together when `TBD_20260725_vsearch_namespace.md` runs. |
| 5. `Makefile.am:28` | **Approved**: replace the "evaluated and left disabled" note with `-Wformat=2` in the debug warning set, in Phase 8. |

**0. Scope — reverse the exclusion, full scope.**
`TBD_20260725_replace_C_functions.md:5` states that the printf family is
not a refactoring target. That line is superseded by this plan and should
be amended to point here in Phase 1, so the two documents do not
contradict each other for the length of the migration.

**1. `double` fields — keep `std::fprintf`, per field.**
167 sites format a `double`; swarm had 4 and kept `fprintf` at all of
them. Writing a fixed-precision `double` byte-identically is not a digit
loop: `printf` converts the *exact* binary value and breaks ties
half-to-even, so the obvious `value * 10^p + 0.5` differs on real
inputs, and C++11 offers neither `to_chars` nor anything else that avoids
`vsnprintf`. So: a converted line keeps exactly one `fprintf`, for its
`double` field(s), and everything around it becomes `fprint*` — which is
what retires the `PRI` macros at those 41 mixed sites. The cost is
accepted with the decision, not discovered later:
`core/results.cpp`'s blast6 line becomes one `fprintf` plus ~10
`fprint_integer`/`fprint` calls, and the same rule applies in the hot
writers, not only the cold ones. If a *specific* hot line reads badly
enough to be worth an exception, raise it in that file's commit rather
than bending the rule tree-wide.

**2. Field width — add the width parameter and the space run.** 47
non-float sites use one: 42 with a constant width (`%10d` mergepairs
stats ×18, `%10u` udbstats ×9, `%5d`/`%8d` chimera ×7, `%6d`
fastq_eestats2 ×2, the align_simd dumps ×5, plus `fasta.cpp`'s `%02x`,
which keeps its `fprintf`) and 5 with a run-time width
(`core/showalign.cpp:206`, `:216`, `:224` and `core/results.cpp:769`,
`:771`, all `%*s`/`%*" PRId64`). Both are exactly reproducible with
spaces, so both are in. Leaving them out would have kept `PRI` alive in
six files (`fastq_mergepairs`, `udbstats`, `fastq_eestats2`, `chimera`,
`showalign`, `results.cpp`), which is most of the point of the exercise.
Bonus: `results.cpp` currently computes those widths with
`std::to_string(qseqlen).size()`, which the digits view gives for free.

**3. Where to stop — `core/align_simd.cpp`'s 24 sites are out of scope.**
They are `dprofile_dump16()` / `dumpscorematrix()` debug dumps whose only
callers sit under `#if 0` (`:541`, `:709`); they carry no `PRI` macro, and
no test or differential run can exercise a converted version. The end
state records them as deliberately left, so a later reader does not read
them as a missed sweep.

**4. Namespace — global now, moved as a group later.** The new overloads
join today's `fprint(View)` in the global namespace, so this migration
introduces no namespace-qualification decisions at ~300 call sites and its
diffs stay mechanical. The trade-off is explicit: this widens the ADL
surface that `TBD_20260725_vsearch_namespace.md` is about, and that pass
now has four `fprint*` names to move instead of one — a bigger but purely
mechanical move, and it must not be dropped. Note it in that document's
inventory when this branch lands.

**5. The `-Wformat=2` comment in `src/Makefile.am:28` — approved.** In
Phase 8, once `fatal()`'s run-time format strings are gone, the
"evaluated and left disabled" note is replaced by `-Wformat=2` in the
debug warning set. Expect the flag to surface further sites the first
time it is enabled; they belong to that phase, not to a follow-up.


## Traps

**The `char const *` overload trap.** Do **not** add
`fprint(std::FILE *, char const *)`. Against the array template, a
string-literal argument gives both candidates an Exact Match sequence
(array-to-pointer is an lvalue transformation, excluded from the
subsequence comparison), and a non-template beats a template
specialization on the tiebreaker — so *every* literal would silently
take the `strlen` path instead of its compile-time bound. The 8
`fputs`-with-a-run-time-pointer sites therefore either go through
`View`, or keep `std::fputs`, or get a differently-named function. Verify
the claim with a two-line compile test before relying on either
direction; swarm records the same dispatch matrix as a test block in
`print_view.hpp` and vsearch should too.

**Thread interleaving.** This is the one hazard swarm did not have. A
single `std::fprintf` is one stdio-locked call; splitting a line into six
`fprint` calls is six, and another thread can now interleave between
them. vsearch's threaded commands already serialize whole output blocks
under their own `mutex_output` (`usearch_global.cpp:467`,
`search_exact.cpp:534`, `sintax.cpp:602`, `allpairs_global.cpp:508`,
`chimera.cpp:1130`/`:1726`/`:2253`, `mask.cpp:232`, and mergepairs' chunk
lock), so the split is safe — but **check the lock covers the whole line
at each converted site**, per file, and say so in the commit. `--threads
1` differential runs cannot catch a regression here; a `--threads 8`
run with a large query set is what would.

**Embedded NUL.** `fwrite` emits the full byte count; `"%.*s"` stops at
an embedded NUL as well as at the precision. `print_view.hpp` already
documents that vsearch's readers reject embedded NULs so the two agree
on every reachable input. Nothing new, but it applies to each newly
converted `%s`/`%.*s`.

**Empty views.** `fprint(View)` returns before `fwrite` when the view is
empty, because an empty view may carry a null pointer (the `--profile`
path in `msa.cpp` does). Converting a `%s` of a possibly-null `char *`
must go through that guard, not around it.

**`src/Makefile.in` is version-controlled.** Adding
`utils/decimal_digits.hpp` means editing `src/Makefile.am` *and*
regenerating and committing `src/Makefile.in`.

**`-fno-exceptions`.** The CLI build sets it (`AM_CXXFLAGS`). The new
headers allocate nothing and throw nothing, so they are also fine for the
library build; keep `to_decimal` `noexcept`.


## Phasing

One family per commit, cheapest verification first. Branch
`tmp_$(date +%Y%m%d%H%M%S)` off `dev`, `Co-Authored-By: Florian
FILLOUX`.

| # | scope | risk | notes |
|---|---|---|---|
| 1 | `decimal_digits.hpp` + the `print_view.hpp` overloads, **no call site converted**; `src/Makefile.am` + the tracked `src/Makefile.in`; amend `TBD_20260725:5`; drop the two dead `<cinttypes>` includes | none | object files must be byte-identical: the templates are uninstantiated |
| 2 | the 116 literal `fputs` and 115 literal `fputc` sites → `fprint` | very low | no format string involved; per-file commits; mechanical |
| 3 | cold, float-free report writers: `udbinfo`, `version`, `help`, `cli.cc`, `dbindex`, `fastx`, `db`, `otutable`, `logfile`, `msa`, `filter`, `cut`, `fastq_join`, `fastx_syncpairs`, `fastx_subsample`, `rereplicate`, `fastq_convert` | low | first `PRI` files retired; `otutable` is also the widest-table hot loop, so measure it |
| 4 | hot record writers: `results.cpp`, `fasta.cpp`, `fastq.cpp`, `derep.cpp`, `derep_prefix`, `derep_smallmem`, `orient`, `sintax`, `cluster`, `search_exact`, `usearch_global`, `allpairs_global`, `progress.hpp` | medium | one file per commit, `hyperfine` per file (see below); the thread-interleaving check lives here |
| 5 | padded tables: `fastq_mergepairs` stats block, `udbstats`, `fastq_eestats`, `fastq_eestats2`, `fastq_chars`, `fastq_stats`, `chimera` labels | medium | uses the width writer from Phase 1; widest diffs, easiest to verify (cold, one table each) |
| 6 | run-time widths: `showalign.cpp` ×3, `results.cpp` ×2 | medium | drops `to_string().size()` width computation |
| 7 | the 30 `snprintf` sites: `linmemalign.cpp` cigar flush (hot), `fasta.cpp` abundance annotation (hot), `getseq.cpp` field buffer, `cluster.cpp` centroid label / cigar, `vsearch.cc`, `cli.cc` | medium | same primitive, different sink; `linmemalign` is the clearest perf case |
| 8 | `fatal()`'s two format-taking overloads (26 callers) and `warn()` → message built by the caller; then `-Wformat=2` in `src/Makefile.am`, replacing its note (approved) | medium | touches the library error path (`fatal_throw.cpp`); the sites `-Wformat=2` surfaces on first enabling belong to this phase |
| 9 | sweep: delete all `<cinttypes>`, prune `<cstdio>` from files that no longer need it (79 include it today), `cppcheck`, `run_clang_tidy.sh`, all four cross-compiles, final counts | low | the commit that states the end state |

Phases 3–7 are independent of each other and can be reordered freely;
2 must follow 1, and 9 must be last.


## Verification (per phase)

1. **Build clean**, debug flags as in `CLAUDE.md`, plus the debug warning
   set (`-Wconversion`, `-Wsign-conversion`, `-Wtype-limits`,
   `-Wdouble-promotion` all matter for a digit loop).
2. **Byte-identical output**, which is the whole safety net here. Use the
   recipe from the View/Span pass: `git worktree add --detach` a baseline
   (vsearch-tests already holds `dev`), run both binaries with
   `--threads 1`, strip the command-line first line, `diff -r`. Port
   swarm's `byte_identity.sh` (it refuses to compare a DEBUG binary
   against a RELEASE one, which is the mistake that produces
   intermittent phantom differences) and drive it from the per-command
   option matrices. Every commit in every phase must be byte-identical;
   there is no ordering-change exception in this migration.
3. **A `--threads 8` run** for each converted threaded writer, diffed
   after sorting, for the interleaving hazard.
4. **Test suite**: `(cd ~/Documents/src/vsearch-tests/ ; bash
   run_all_tests.sh ../vsearch/bin/vsearch | grep FAIL)`. Two cases flake
   under concurrent build load (getseqs `--labels` log, mergepairs
   `--threads 1024`) — re-run those two in isolation before believing
   them.
5. **Perf, release build, `hyperfine`**, only for phases 4 and 7, and
   only for the command whose writer changed: `derep_id`,
   `usearch_global`, `cluster_size`, `fastq_chars`, `cut`. Expect neutral
   to slightly better; a regression means the split introduced a
   per-field call where the format parse was previously amortised over a
   long line, which is a signal to batch that line differently, not to
   revert the phase.
6. **`cppcheck` and `run_clang_tidy.sh`** on each touched file (both
   report *current* warnings, not a diff — compare against the baseline
   commit).
7. **Cross-compile all four targets** (MinGW, POWER, RISC-V/MIPS,
   plus native ARM if available) at the end of each phase that adds a
   template instantiation. MinGW is the interesting one: it is the reason
   the `PRI` macros were there.


## End state

**Projected, before the work** (kept as written; the measured column was
added on completion, and three of these were wrong):

| | projected | actual |
|---|---|---|
| `PRIu64` / `PRId64` | 0 uses | **0** |
| `<cinttypes>` | 0 includes | **0** |
| `std::fputs` | 0, or 8 if the run-time-pointer sites keep it | **107** |
| `std::fputc` | 0 | **2** (1 is `align_simd`, 1 is `fprint(char)`'s own wrapper) |
| `std::fprintf` / `std::printf` | ~167 | **257** |
| `std::snprintf` | not projected | **4** |
| `-Wformat=2` | on, note replaced | **on, and clean** |
| new/extended headers | two | **three** |

Why the three differ, all decided and recorded in the commits that caused
them:

- **`std::fputs` 107, not 0–8.** The projection assumed the run-time-pointer
  sites were the only ones. They were not: every `"%s"` of a `char const *`
  whose length the caller does not know became one, and the Traps section
  above is what rules out the alternative (an `fprint(std::FILE *,
  char const *)` overload would silently capture every string literal from
  the array template). They carry no format string, which is the point; the
  count rose because `"%s"` sites were converted *into* them.
- **`std::fprintf` 257, not ~167.** Arithmetic, not scope creep: the
  projection counted *sites* that format at least one `double`, and
  Decision 1 gives each `double` *field* its own call, so a line with three
  of them becomes three calls.
- **Three headers, not two.** `utils/decimal_digits.hpp` and the extended
  `utils/print_view.hpp` were planned; `utils/print_record.hpp` was added
  late, when phase 4's measurement showed the split is only a win where the
  compiler inlines the writer. It is the "batch that line differently"
  the verification section prescribes.

Deliberately left, so a later reader does not read them as a missed sweep:

- **`core/align_simd.cpp`, 27 sites** (24 `fprintf`/`printf`, 2 `snprintf`,
  1 `fputc`) — Decision 3. Their only callers sit under `#if 0`, they carry
  no `PRI` macro, and no differential run can exercise a converted version.
  The plan said 24; the extra 3 are the `snprintf`/`fputc` in the same dumps,
  which the original census counted separately.
- **the three hex sites** — `core/fasta.cpp`'s `"%02x"`,
  `commands/version.cpp`'s `"%lx"` and `core/fastx.cpp`'s `"0x%2x"`. The
  first two were named in Out of scope; the third is the same case, found
  during phase 8, and is noted at the call.
- **two `snprintf`** outside `align_simd` — `vsearch.cc`'s and
  `core/fastx.cpp`'s, each formatting one field a digit loop cannot produce
  (a `"%.1f"` and that `"0x%2x"`).

All three headers are still in the global namespace and are listed in
`TBD_20260725_vsearch_namespace.md` — twelve `fprint*` overloads across
`print_view.hpp` and `print_record.hpp`, which must move together, plus
`Record`/`OutputRecord` and the `decimal` namespace.


## Out of scope (worth their own plans)

- **`std::to_string`, 43 uses.** On libstdc++ ≤ 10 these *are* a
  `vsnprintf` call with a format string, so "it starts with `std::`" is
  not evidence the C library is not being called. Five in `results.cpp`
  exist only to measure a width and die with Phase 6; the rest build
  error messages and are cold. swarm converted its one hot case
  (`cigar.cpp`) in the same pass.
- **`std::FILE *` itself.** Deliberate, RAII-wrapped
  (`utils/open_file.hpp`), already checks `fflush`/`ferror` before
  `fclose`, and `std::fstream` was evaluated and rejected for the
  record-by-record readers. Not a target.
- **The two hex sites** (`core/fasta.cpp:295` `%02x`,
  `commands/version.cpp:92` `%lx`): one error path, one zlib version
  line, no `PRI` macro either side. Leave them, documented.
- **`assert()`**: `CLAUDE.md` asks for more of it, not less.


## Results, measured on the migration branch

Recorded as each phase landed, so a later reader does not have to re-derive
them. Release builds (`-O2`, this toolchain predefines `_FORTIFY_SOURCE 3`
at `-O2`), `hyperfine`, `--threads 1`.

### Codegen (phase 2)

Checked in the disassembly at vsearch's release flags, which settles which
of these calls are worth replacing:

| written as | compiles to |
|---|---|
| `std::fputs("\t*\n", h)` | `jmp fwrite` |
| `fprint(h, "\t*\n")` | `jmp fwrite` — identical |
| `std::fputc('\t', h)` | `jmp fputc` |
| `fprint(h, '\t')` | `jmp fputc` — identical |
| `std::fprintf(h, "\t*\n")` | `jmp __fprintf_chk` — **not** folded |
| `std::fprintf(h, "%" PRIu64, v)` | `jmp __fprintf_chk` |
| `fprint_integer(h, v)` | inlined digit loop + `call fwrite` |

So the `fputs`/`fputc` family is a readability change with no performance
component, and every `fprintf` — even one whose format has no conversion at
all — pays a run-time format parse that fortification prevents GCC from
folding away. That is the cost phases 3–7 remove.

### `core/otutable.cpp` (phase 3)

The widest table vsearch writes. Measured by differencing a run with the
three table outputs against the same run without them, so the clustering
that dominates the command is subtracted out rather than diluting the
result: 3468 sequences, each its own sample, `--id 1.0`, giving a 24 MB
classic table.

| | without tables | with tables | table writers |
|---|---|---|---|
| `dev` | 4.681 s ± 0.048 | 5.396 s ± 0.072 | 0.715 s |
| branch | 4.671 s ± 0.063 | 5.248 s ± 0.085 | 0.577 s |

**1.24× on the writers**, 1.028× on the whole command. The
without-tables pair is identical to within noise, which is the control: no
code on the clustering path changed.

### The hot record writers (phase 4)

End-to-end, on real inputs, the phase is neutral — which is the headline,
because these are the commands users run:

| command | `dev` | branch |
|---|---|---|
| `--usearch_global` 8652 q / 57605 db, `--alnout --blast6out --uc` | 15.141 s ± 0.250 | 14.803 s ± 0.063 |
| `--cluster_size` 17314 seqs, `--centroids --uc` | 7.636 s ± 0.093 | 7.621 s ± 0.105 |
| `--derep_fulllength` 286789 seqs, `--output --uc` | 597.2 ms ± 8.4 | 588.8 ms ± 10.9 |
| `--derep_fulllength`, no `--uc` (control) | 473.9 ms ± 8.5 | 476.1 ms ± 3.8 |

The `--uc` writer differenced out of the last two pair is 123.3 ms → 112.7 ms,
1.09× faster. `--usearch_global` is dominated by alignment: 137k hits are
written in a 15 s run, so the writer is under 2 % of it either way.

**Where it is not neutral, and why.** Isolated on a writer-dominated
synthetic (300 k `--search_exact` queries against a one-sequence database,
so the lookup is trivial), the picture splits:

| output | `dev` | branch | |
|---|---|---|---|
| `--userfields id` (one field; **code path unchanged**) | 808.9 ms ± 11.5 | 855.3 ms ± 7.3 | +5.7 % |
| `--blast6out` (12 fields) | 839.4 ms ± 16.5 | 923.2 ms ± 9.2 | +10.0 % |
| `--uc` (10 fields) | 831.3 ms ± 12.0 | 894.0 ms ± 11.6 | +7.5 % |
| `--userout`, 36 fields | 1.248 s ± 0.023 | 1.220 s ± 0.010 | −2.2 % |

Two separate effects, and they are worth keeping apart:

1. **A baseline shift of +5.7 % on a path the migration did not change.**
   `--userfields id` writes one `"%.1f"` per record through the same
   `std::fprintf` in both binaries. `.text` grew 3.7 % (764 258 → 792 186 B)
   and `results_show_userout_one` grew 9 %, so this is code layout, not
   writer cost. Bisected to phase 3, which touches nothing on this path.
2. **`fprint_integer` is not inlined.** In
   `results_show_blast6out_one` GCC emits `call <void
   fprint_integer<int>(FILE*, int)>` nine times rather than inlining the
   digit loop nine times: the function went from 305 to 529 bytes and from
   5 stdio calls to 24. Nine out-of-line calls at ~14 ns is the ~125 ns per
   record measured. It is *only* the call: with the same code inlined, a
   microbenchmark of the same ten-field tail gives

   | form | ns/record |
   |---|---|
   | one `std::fprintf`, 10 fields (`dev`) | 173.5 |
   | split into one call per field, inlined | 137.1 |
   | the record batched into one `fwrite` | 95.7 |

   so the split is 1.27× *faster* when inlined, and batching is 1.8× faster.
   `--userout` shows the inlined-shaped win (−2.2 %) because it was already
   one `fprintf` per field, so the split removed 36 format parses per record
   and added no calls.

**Resolved, and both numbers above are wrong.** The plan says a regression
here "is a signal to batch that line differently"; `utils/print_record.hpp`
is that batching, and it lands on this branch. Two measurement errors had to
be found first, and they were pulling in opposite directions — see
"Corrected measurements" at the end of this document. In short: the
benchmark used above was 75 % DUST masking rather than writing, which
diluted every writer effect, and a code-alignment artefact was adding ~34 ms
of pure noise on top. With both removed, every converted writer is faster
than `dev`.

### Final numbers (phase 9, all nine phases landed)

Release builds of `dev` and of the branch tip, same host, `hyperfine`.

| what | `dev` | branch | |
|---|---|---|---|
| `--otutabout --biomout --mothur_shared_out`, 24 MB table | 5.414 s ± 0.047 | 5.312 s ± 0.049 | |
| the same run without the tables (control) | 4.667 s ± 0.067 | 4.662 s ± 0.074 | identical |
| → the three table writers, differenced | 0.747 s | 0.650 s | **1.15×** |
| `--usearch_global`, real data, `--alnout --blast6out --uc` | 15.141 s ± 0.250 | 14.803 s ± 0.063 | neutral |
| `--cluster_size`, real data, `--centroids --uc` | 7.636 s ± 0.093 | 7.621 s ± 0.105 | neutral |
| `--derep_fulllength --output --uc`, 287 k seqs | 579.4 ms ± 5.5 | 584.6 ms ± 6.2 | neutral |
| `--cut`, 58 k seqs | 78.6 ms ± 1.4 | 79.9 ms ± 1.6 | neutral |
| `--fastq_chars` | 870 µs ± 130 | 843 µs ± 80 | neutral |
| `--blast6out`, 300 k records, one-sequence db | 842.0 ms ± 1.8 | 921.4 ms ± 3.1 | **+9.4 %** |

Every real command is neutral or better; the table writers are measurably
faster; the one regression is the `--blast6out`/`--uc` record on a synthetic
that does nothing but write, and its cause is known and written up under
phase 4 (`fprint_integer` is emitted out of line when one function calls it
nine times). That is the single open decision on this branch.

### Verification actually performed

- **byte-identical output**: 441 comparisons over a 144-run option matrix
  (`byte_identity.sh`), at every commit, at `--threads 1` and `--threads 8`,
  against both a DEBUG and a RELEASE `dev` baseline. One deliberate
  difference, in its own commit: the `cli.cc` out-of-bounds crash.
- **outside the matrix**, because it cannot reach them: every reachable
  diagnostic (11 error paths, diffed against `dev`), `--sff_convert` on a
  3 MB SFF slice with and without `--sff_clip`, `--clusters` over a 5-cluster
  run, a 1500 nt alignment at `--rowlen` 0/3/17/60/200, and the library-only
  `centroid_label` truncation at 8/1022/1023/1024/1025/4000 characters
  against a `dev`-built `libvsearch.a`.
- **test suite**: 9757 PASS, 0 FAIL. (`--derep_id --relabel --sizeout` also
  flakes under concurrent build load, like the two cases `CLAUDE.md` names;
  241 PASS in isolation.)
- **api_examples**: 40 PASS, 0 FAIL against the branch's release library.
- **primitives**: `to_decimal` against `std::snprintf` over 400 k+ values
  including `UINT64_MAX` and `INT64_MIN`, and `fprint_integer(value, width)`
  against `"%*d"` for every width 0–23, under
  `-fsanitize=undefined,address -D_GLIBCXX_DEBUG`.
- **`-Wformat=2`**: on, and clean.
- **cross-compiles**: MinGW x86_64, POWER (ppc64le), MIPS64el and aarch64 all
  build, with only the two pre-existing warnings (libstdc++'s `-Warray-bounds`
  false positive in `attributes.cpp`, and AltiVec compound literals in
  `align_simd.cpp`). MinGW was the interesting one: it is the reason the `PRI`
  macros existed, and it no longer needs them.
- **cppcheck** clean on every touched file. **clang-tidy** reports one
  `cert-err33-c` on `print_view.hpp`'s literal overload; it is a check
  limitation, not a defect — the same `static_cast<void>(std::fwrite(...))`
  is *not* reported outside a template, verified with a two-function probe,
  and the existing `fprint(View)` overload has carried it since before this
  branch.

## Corrected measurements (after landing `print_record.hpp`)

The phase-4 and phase-9 numbers above are superseded. Two errors in how they
were taken, both found while investigating the 34 ms:

**1. The "writer-dominated" benchmark was 75 % DUST masking.**
`--search_exact` masks every query with DUST by default, and `callgrind` puts
74.65 % of that run's instructions in `core/mask.cpp`'s `wo()`. So every
writer effect quoted above is diluted by a factor of about four, in both
directions. Adding `--qmask none --dbmask none` drops the run from 1.058 G to
284 M instructions and makes the printf machinery visible in the profile
(`__printf_buffer` + `__printf_buffer_write`, 12 % combined).

**2. A ~34 ms code-alignment artefact, not a defect.** Bisected to `0628b38c`
and, within it, to `core/fastx.cpp` alone — whose only changes there are in
`warn()` and `report_stripped_warning()`, both of which run once per command
and neither of which this benchmark calls. `callgrind` settles it: the two
builds execute **1,058,298,400 and 1,058,258,396 instructions**, and the
*slower* one executes 40 k fewer. There is no extra work.

The mechanism: `warn()` grew 54 bytes, which shifted everything after it in
link order by 0x150 bytes, and `wo()` — byte-identical, size `0x19e` in both
— went from `0x697c0` (64-byte aligned) to `0x69910` (`mod 64 == 16`). A
74-%-of-runtime loop lost its cache-line alignment. Confirmed causally:
rebuilding both with `-falign-functions=64` puts `wo()` at `mod 64 == 0` in
both and collapses the gap from 33 ms to 3.7 ms, i.e. to noise.

This is worth remembering when benchmarking vsearch at all: any change that
alters the size of any function can move a hot loop's alignment and produce a
few percent in either direction with no change in work done. `-falign-functions=64`
costs a little `.text` and makes such measurements repeatable; whether to
adopt it is a separate decision, not taken here.

### The writers, measured properly

300 k records, one-sequence database, `--qmask none --dbmask none`,
`--threads 1`, release. The control's writer path is byte-for-byte identical
in all three binaries.

| output | `dev` | branch (unbatched) | branch + `Record` | vs `dev` |
|---|---|---|---|---|
| control (`--userfields id`) | 150.3 ms ± 2.7 | 153.2 ms ± 2.4 | 150.1 ms ± 1.2 | — |
| `--blast6out` | 181.2 ms ± 3.2 | 223.3 ms ± 2.3 | **161.6 ms ± 3.0** | **1.12× faster** |
| `--uc` | 176.1 ms ± 3.4 | 194.8 ms ± 1.2 | **167.3 ms ± 3.9** | **1.05× faster** |
| `--samout` | 234.7 ms ± 4.8 | 292.3 ms ± 2.1 | **190.8 ms ± 4.0** | **1.23× faster** |
| `--alnout` | 439.1 ms ± 5.8 | 571.3 ms ± 5.3 | **335.5 ms ± 7.2** | **1.31× faster** |
| `--userout`, 36 fields | 586.3 ms ± 4.2 | **506.3 ms ± 4.2** | 507.3 ms ± 4.3 | **1.16× faster** |

Subtracting the control to leave the writer alone: `--blast6out` 30.9 → 11.5 ms
(**2.7×**), `--samout` 84.4 → 40.7 ms (**2.1×**), `--alnout` 288.8 → 185.4 ms
(1.56×), `--uc` 25.8 → 17.2 ms (1.5×), `--userout` 436.0 → 357.2 ms (1.22×,
from the split alone — it is not batched).

So the unbatched split was worse than first reported (+23 % on `--blast6out`,
+30 % on `--alnout`, not +9 % and +15 %), and the batched form is better than
first reported. `--userout` needed no buffer: it was already one `fprintf` per
field, so the split removed 36 format parses per record and added no calls.
