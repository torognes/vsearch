# DONE 2026-08-08 — eliminate `&parameter` at call sites

**Status: implemented on `tmp_20260806093624`, 10 commits (87435bd3..2b6701ed),
awaiting human review.** What actually happened is recorded in "Outcome" at
the bottom; the body below is the plan as approved, left intact so the two can
be compared.

Issue (closed 2026-08-07 15:06, Torbjørn-era code):

> fix all function calls using `&parameter`
>
> This was introduced by Torbjørn and it overcomplicates reference usage.
> A simpler solution is to use `T& parameter` in the function declaration.

Same refactoring was done in swarm over the past two days (see
`~/Documents/src/swarm`, e.g. `0e220ec line_buffer: pass read_one_line's
in/out parameters by reference` and `83582ed search: take the always-live
kernel vectors by reference`).


## What this is really about

`&x` at a call site is the *symptom*. The cause is a function that declares
`T * param` where it means "a single object that is always there". The fix is
always on the **declaration** side — change `T *` to `T &` — and the `&` at
every call site then disappears on its own. Callers that already held a
pointer gain a `*` instead, which is the correct trade: one explicit
dereference at a site that provably has an object, in exchange for deleting
the whole "can this be null?" question from the callee.

So this plan is organised by *callee*, not by call site.


## Census

The mechanical count is misleading, because `&` in vsearch is overwhelmingly
the bitwise operator (masks in `align_simd.cpp`, `bitmap.cpp`,
`cpu_features.cpp`, `unique.cpp`, `maps.cpp`), a lambda capture, a
reference-to-array declaration (`char const (&literal)[Size]`), or
address-of-function (`decltype(&CityHash64)`, `&std::fclose`). After
filtering those out:

| group | what | functions | call sites | `&` args removed |
|---|---|---|---|---|
| **A** | scalar out-parameters | 7 | 16 | **37** |
| **B** | whole-object pointer parameters | 10 | 30 | **19** |
| **C** | must stay pointers | — | ~60 | 0 |
| **D** | optional / measure first | 3 | 7 | 6 |

Groups A and B together remove **56** `&` arguments and, more importantly,
delete 17 nullable-looking parameters that were never nullable.

(Group B excludes `xstat`/`xfstat` — see B8, dropped by decision.)


## Group A — scalar out-parameters (do these first)

These are unambiguous: a `T *` that exists only so the callee can write one
value back. Every call site passes `&local`. No call site passes `nullptr`.

### A1. `LinearMemoryAligner::alignstats` — the big one

`core/linmemalign.hpp:172`

```cpp
auto alignstats(char const * cigar, View<char> a_sequence, View<char> b_sequence,
                int64_t * nwscore, int64_t * nwalignmentlength,
                int64_t * nwmatches, int64_t * nwmismatches,
                int64_t * nwgaps) -> void;
```

Five out-params × four call sites = **20** of the 37 `&` args in group A:
`core/chimera.cpp:2198`, `commands/allpairs_global.cpp:444`,
`core/cluster.cpp:782`, `core/searchcore.cpp:834`.

**DECIDED — return a struct, not references.**

```cpp
struct AlignStats {
  int64_t score = 0;
  int64_t alignmentlength = 0;
  int64_t matches = 0;
  int64_t mismatches = 0;
  int64_t gaps = 0;
};

auto alignstats(char const * cigar, View<char> a_sequence,
                View<char> b_sequence) const -> AlignStats;
```

Rationale: five adjacent parameters of the *same* type is exactly the
"easily swappable function parameters" hazard CLAUDE.md calls out, and the
plain-reference form does not fix it — five `int64_t &` in a row are just as
swappable as five `int64_t *`, minus the `&` that currently makes a
misordered argument at least visually odd. The struct makes them named
fields, lets the four call sites drop five local declarations each, and
probably lets the method become `const`. It is a bigger diff (~4 call sites
× ~12 lines) but a strictly better endpoint.

Caveat: at `core/cluster.cpp:782` and
`commands/allpairs_global.cpp:444` the five locals are *shared* with the
`search16` path (the SIMD aligner fills them in the fast path, `alignstats`
only in the fallback). So the call sites become
`auto const stats = lma.alignstats(...); nwscore = stats.score; ...` rather
than a clean replacement, unless D1 lands too. Still a net win; just do not
expect the locals to vanish in this phase.

### A2. cigar run-length readers — `char const **` → `char const * &`

`utils/cigar.hpp:74` and `:81`

```cpp
auto read_runlength(char const * first_character, char const ** first_non_digit) -> long long;
auto find_runlength_of_leftmost_operation(char const * first_character,
                                          char const ** first_non_digit) -> long long;
```

Call sites: `utils/cigar.cpp:167`, `core/searchcore.cpp:377`, `:408`,
`:645`, `core/linmemalign.cpp:730`, `core/results.cpp:876`, and
`core/msa.cpp:376`.

`core/msa.cpp:376` is the one that pays for this phase on its own. Today it
manufactures a pointer-to-pointer just to satisfy the signature, then reads
back through two levels:

```cpp
auto ** next_operation = &position_in_cigar;
auto const runlength = find_runlength_of_leftmost_operation(position_in_cigar, next_operation);
auto const operation = **next_operation;
```

With `char const * & first_non_digit` that collapses to a plain local and a
single `*`.

Inside `read_runlength` the body changes from `*first_non_digit =
end_of_digits;` to `first_non_digit = end_of_digits;`. The `std::strtoll`
call one line above keeps its `&end_of_digits` — that is a C library API
(group C).

**Comments — DECIDED.** Reword `utils/cigar.cpp:132-134` only; leave
`:146` (`// duplicate: msa.cc`) exactly as it is.

The `:132-134` comment explains why the local `end_of_digits` exists. The
reasoning survives the change, but its last clause names the caller-facing
type, which is what this phase alters:

```cpp
// strtoll's endptr is a char** even though it only ever points back into the
// read-only input; take it in a local and assign through the caller's
// reference-to-pointer-to-const, so callers holding a const buffer do not
// have to launder it themselves
```

Only the "hand the caller a pointer to const" clause changes — the
`strtoll`-is-a-`char**` observation stays, because the `std::strtoll` call
below still takes `&end_of_digits` and that is still why the local is
needed at all.

### A3. `buffer_filter_extend` — `bool * ok, char * illegal_char`

`core/fastq.cpp:239`. Two call sites, `:414` and `:525`, both
`&ok, &illegal_char`. → `bool & ok, char & illegal_char`.

### A4. `tax_parse` — `int * tax_start, int * tax_end`

`core/tax.cpp:76`, anonymous namespace, single call site `:154`.
Smallest possible change; good first commit to validate the workflow.

### A5. `merge_sym` — `char * sym, char * qual`

`core/mergepairs.cpp:219`, anonymous namespace, single call site `:316`.
Body uses `* sym = …` throughout; becomes `sym = …`.

### A6. `scan_matches` — `int * best_start, int * best_len`

`core/chimera.cpp:494`, single call site `:610`. Note the same function's
`int const * matches` parameter is a genuine **array** and must stay a
pointer — convert only the last two.


## Group B — whole-object pointer parameters

Here the parameter is a real object, not an out-slot. Each entry below has
been checked for `nullptr` call sites; all are non-null at every site.

### B1/B2. `align_trim` and `search_acceptable_aligned` — `struct hit *`

`core/searchcore.hpp:259` and `:262`. Do these two together: they are
called back-to-back on the same `hit` at three of four sites.

- `core/searchcore.cpp:867`, `:870` — pass `&hit` (a local `struct hit`)
- `core/cluster.cpp:823`, `:836` — pass `hit` (pointer into an array,
  dereferenced unconditionally on the lines around it)
- `commands/allpairs_global.cpp:487`, `:490` — same shape
- `commands/search_exact.cpp:202` — passes `hp`

Net: two `&` disappear, five call sites gain a `*`. That looks like a wash
on character count and is not — it removes "may be null" from two
signatures in the hottest path in the program.

`search_acceptable_aligned` already takes `struct searchinfo_s const &` as
its first parameter, so this also makes the signature internally
consistent.

### B3/B4. the hit comparators

`core/searchcore.cpp:141` `hit_compare_byid_typed`, `:190`
`hit_compare_bysize_typed`, `commands/allpairs_global.cpp:123`
`allpairs_hit_compare_typed` — all `struct hit const *` pairs.

Call sites: `core/searchcore.cpp:974`, `:984`, `:1008`, `:1018`, `:1055`;
`commands/allpairs_global.cpp:501`.

**The one place to be careful.** Four of those sites read:

```cpp
if ((best == nullptr) or (hit_compare_byid_typed(&hit, best) < 0))
```

`best` is genuinely nullable. This is safe because `or` short-circuits, so
after conversion the call becomes `hit_compare_byid_typed(hit, *best)` and
`*best` is only evaluated when `best != nullptr`. Verify the short-circuit
survives at each of the four sites before committing — a later refactor
that turns the guard into an `if/else if` would make this a null deref.
Consider an `assert(best != nullptr)` inside the `else` arm to make the
contract explicit, per CLAUDE.md.

The two comparators are `inline` in the `.cpp`; the sort predicates at
`:1055` and `:501` become the natural `return hit_compare_byid_typed(lhs,
rhs) < 0;`.

### B5. `cluster_query_init` / `cluster_query_exit` — `struct searchinfo_s *`

`core/cluster.cpp:230` and `:278`. Twelve call sites, all in `cluster.cpp`
(no header declaration exists — these are file-local in practice):

- `:329`, `:334`, `:344`, `:345` — pass `&si` (the four that motivated the issue)
- `:1053`, `:1056`, `:1133`, `:1136` — pass `si_p.data()` / `si_m.data()`
- `:1705`, `:1711`, `:1948`, `:1953` — pass `cs->si.get()` / `cs->si_minus.get()`

The `.data()` sites want `si_p.front()`, not `*si_p.data()` — and that
change is worth making for its own sake, because `.data()` on a
`std::vector` passed to a single-object parameter is precisely the
ambiguity this refactoring exists to remove. **Confirm the vectors are
non-empty at those four sites** (they are sized from the thread count
upstream, but check) and add an `assert` if it is not locally obvious. The
`.get()` sites become `*cs->si`.

**Bonus — NOT NEEDED; the plan was wrong here.** Phase 7b was scheduled on
the premise that these two functions had external linkage. They do not:
both sit inside the anonymous namespace that spans `cluster.cpp:175-290`,
which the census missed because it looked for a `static` keyword and a
header declaration and found neither. `nm` settles it:

```
$ nm -C src/core/cluster.o | grep cluster_query
0000000000000642 t (anonymous namespace)::cluster_query_exit(searchinfo_s&)
00000000000002f4 t (anonymous namespace)::cluster_query_init(searchinfo_s&, ...)
```

Lowercase `t` is a local text symbol. Phase 7b dropped, no commit.

### B6. `query_exit` — `struct searchinfo_s *`

`core/chimera.cpp:2014`, single call site `:2089` passing `&a_search_info`.

### B7. `Dbhash::search_first` / `search_next` — `struct dbhash_search_info_s *`

`core/dbhash.hpp:112` and `:115`. Four call sites: `core/dbhash.cpp:172`,
`:175`, `commands/search_exact.cpp:217`, `:221`. All pass `&info` on a
local. Both methods are already `const`, so this is contained.

Check at implementation time whether `core/dbhash.hpp` is in the include
closure of `vsearch_api.h`; if it is, bump `VSEARCH_API_VERSION_MINOR`
(currently 0.18.0). Per the standing ABI posture this is not a reason to
delay or reorder anything.

### B8. `xstat` / `xfstat` — **DROPPED by decision**

`os/system.hpp:85-86`, five call sites (`commands/udbinfo.cpp:85`,
`core/fastx.cpp:356`, `core/udb.cpp:156`, `:212`, `core/getseq.cpp:170`).

Keeping the C shape: the pointer form deliberately mirrors POSIX
`stat(2)`/`fstat(2)`, which these wrap, and matching the underlying API is
a real readability argument. These five `& fs` / `& file_status` call sites
stay as they are, and phase 10 disappears from the schedule. Nothing else
in the plan depended on B8.

Recorded here so a future sweep does not re-open it — this is a deliberate
keep, not an oversight. Effectively a group C entry.


## Group C — leave alone, and why

Record these in the review so the next sweep does not re-litigate them.

1. **C library / OS APIs.** `std::strtoll`, `std::strtod`, `std::fread`,
   `std::frexp`, `std::localtime`, `std::memcpy`, `getopt_long_only`,
   `regcomp`/`regfree`/`regexec` (`core/otutable.cpp`), `SHA1_*`/`MD5_*`
   (`utils/sequence_digest.cpp`), `sysinfo`, `sysctl`, `sysctlbyname`,
   `getrusage`, `posix_memalign`, `GlobalMemoryStatusEx`, `GetSystemInfo`,
   `GetProcessMemoryInfo`, and the bzip2 function pointers in
   `os/dynlibs.cpp`. Not ours to change.

2. **`chimera_process_query`'s `struct chimera_cli_state_s * cli`**
   (`core/chimera.cpp:2100`). Looks like a target — `:2321` passes
   `&state` — but it is **nullable by design**: the library path at `:2998`
   passes `nullptr` so the detection core writes to `ci->result_out`
   instead of to files, and the contract is documented at `:234`. Keep the
   pointer. (Same trap as the `*_one` writers in the container-propagation
   campaign.)

3. **`search16`'s output parameters** (`core/align_simd.hpp:101`). `CELL *
   pscores`, `unsigned short * paligned`, … are **arrays** sized by the
   `sequences` argument. Two of the four callers happen to pass `N == 1`
   and therefore `&scalar` (`core/cluster.cpp:760-766`,
   `commands/allpairs_global.cpp`), but the other two pass `.data()` on
   real vectors (`core/searchcore.cpp:787`, `core/chimera.cpp:2167`).
   Converting to references would be wrong. See D1 for the correct fix.

4. **`v_store(VECTOR_SHORT * ptr, …)`** (`core/align_simd.cpp:177`, `:238`,
   and the x86/generic backends). Wraps `vec_st` / `vst1q_s16` /
   `_mm_store_si128`, all of which take pointers; the pointer form also
   carries the alignment expectation. Keep.

5. **Address-of-function and reference-to-array.**
   `decltype(&hash_cityhash64)`, `decltype(&CityHash64)`, `&std::fclose`,
   `&std::vsnprintf`, `char const (&literal)[Size]`. Not call arguments at
   all; they matched the grep, nothing more.

6. **Bitwise `&`.** The bulk of raw grep hits. Ignore.


## Group D — optional, measure first, separate branch

Do **not** bundle these with A/B. They are listed so the census is complete.

- **D1. `search16` array parameters → `Span`.** The right fix for C3: ten
  parameters, four call sites, and it would let the two `N == 1` callers
  stop taking addresses of scalars. Also unlocks the clean version of A1b.
  Sizeable API change to the SIMD entry point; wants its own plan.
- **D2. `aligncolumns_first` / `aligncolumns_rest`'s `& h_min, & h_max`**
  (`core/align_simd.cpp:1731`, `:1980`). Genuine scalar out-params of
  `VECTOR_SHORT` type, so technically group A — but they sit in the
  innermost alignment loop and the type is a SIMD register across four
  backends. Only with a callgrind instruction-count check on all
  architectures.
- **D3. `commands/fastx_syncpairs.cpp:169-170`.** `for (auto * pair : {&
  outfiles.synced_fwd, & outfiles.synced_rev, …})` — an
  array-of-pointers-to-members idiom, not a call argument.
  `std::reference_wrapper` would express it better. Cosmetic, low value.


## Phasing

Branch from `dev` per CLAUDE.md (note: the tree is currently on
`tmp_20260806093624`, so `git checkout dev` first):

```sh
git checkout dev
git checkout -b "tmp_$(date +%Y%m%d%H%M%S)"
```

One commit per numbered item, in this order. Each is independently
revertible and independently testable. All four open questions are
resolved, so nothing blocks the start.

| # | phase | items | risk |
|---|---|---|---|
| 1 | single-call-site out-params | A4, A5, A6 | trivial |
| 2 | `buffer_filter_extend` | A3 | trivial |
| 3 | cigar readers (+ the `msa.cpp` double-indirection, + the `:132-134` reword) | A2 | low |
| 4 | `alignstats` → `AlignStats` struct return | A1 | medium (4 sites, hot path) |
| 5 | hit comparators | B3, B4 | **low but attention** — short-circuit |
| 6 | `align_trim` + `search_acceptable_aligned` | B1, B2 | low |
| 7 | `cluster_query_init` / `_exit` signatures | B5 | low (12 sites, mechanical) |
| 7b | those two into an anonymous namespace | B5 bonus | low; **own callgrind run** |
| 8 | chimera `query_exit` | B6 | trivial |
| 9 | `Dbhash::search_*` | B7 | low; check API closure |

Phases 1-3 and 8 are safe to do back-to-back without rebuilding between
each. Phases 4-7 touch `searchcore.cpp` / `cluster.cpp` / `chimera.cpp` —
build and test after each. Phase 7b must land as its own commit after 7,
never squashed into it.

Nine commits (ten with 7b), removing 56 `&` arguments and 17 nullable-looking
parameters.


## Verification

Per phase:

```sh
cppcheck --force --enable=warning,style --language=c++ --std=c++11 <file>
bash run_clang_tidy.sh <file>          # before and after; it is not a diff
```

After phases 4 and 7, and once at the end:

```sh
make clean && ./configure --enable-debug CFLAGS="-O0 -ggdb3" CXXFLAGS="-O0 -ggdb3"
make ARFLAGS="cr" -j
(cd ~/Documents/src/vsearch-tests/ ; bash run_all_tests.sh ../vsearch/bin/vsearch | grep "FAIL")
```

Two known-flaky scripts under concurrent build load (`getseqs --labels`
log, `mergepairs --threads 1024`) — re-run in isolation against a baseline
binary before believing a FAIL.

At the end, all four cross-compiles (mingw, POWER, RISC-V target as
configured in CLAUDE.md) — phase 10 in particular changes a file that only
one of them compiles.

### Performance

Expectation: **exactly zero** instruction-count change. A `T &` parameter
and a `T *` parameter are the same thing at the ABI level; the only real
codegen difference would come from the aliasing/null-check assumptions the
optimizer can now make, which should be neutral-to-positive.

So the acceptance criterion is the instruction count, **not** wall-clock:

```sh
valgrind --tool=callgrind --callgrind-out-file=cg.out ./bin/vsearch <args>
callgrind_annotate cg.out | grep -E '^ *[0-9,]+ ' | head -8
```

A few-percent hyperfine delta on `usearch_global` / `cluster_fast` with an
unchanged instruction count is code *placement*, not work — every function
in this plan changes size slightly, which shifts link order. Do not chase
it. If a real delta appears (instruction count actually moved), A1b in
phase 4 is the only item in this plan that changes the amount of work done,
so start there.

Benchmarks worth running at the end, since the touched code is the
alignment and search core: `usearch_global.sh`, `cluster_size.sh`. Add
`--qmask none --dbmask none` if measuring anything output-shaped, or DUST
masking will dominate the profile.


## Outcome

Ten commits on `tmp_20260806093624`, `87435bd3..2b6701ed`, 17 files,
+269/−277 lines. Every planned phase landed except 7b (see B5, which the
plan got wrong). One commit was added that the plan did not foresee.

| phase | commit | note |
|---|---|---|
| 1 | `87435bd3` | tax_parse, merge_sym, scan_matches |
| 2 | `491ad7c3` | buffer_filter_extend |
| 3 | `cfe59b0d` | cigar readers + msa.cpp + the `:132-134` reword |
| 4 | `66dd6479` | alignstats → AlignStats |
| 4b | `2b6701ed` | **not in the plan** — see below |
| 5 | `93636cb9` | hit comparators |
| 6 | `49663393` | align_trim + search_acceptable_aligned |
| 7 | `eaacac74` | cluster_query_init/_exit |
| 7b | — | **dropped: the plan's premise was false** |
| 8 | `9841e06c` | chimera query_exit |
| 9 | `9bdfe5d7` | Dbhash::search_first/_next |

### Two places the plan was wrong

1. **Phase 7b was based on a false premise.** `cluster_query_init` and
   `cluster_query_exit` already have internal linkage — they sit inside the
   anonymous namespace spanning `cluster.cpp:175-290`. The census concluded
   otherwise because it searched for a `static` keyword and a header
   declaration and found neither. `nm -C src/core/cluster.o` shows both as
   lowercase `t`. No commit; nothing to do.

2. **The census undercounted, twice.** The first table said group A was 32
   `&` args and group B 17; recounting during the write-up gave 37 and 24
   (19 after dropping B8). The per-item breakdowns were right throughout —
   only the summary totals were wrong. Actual removal, measured from the
   diff: 50 lines carrying a `&` call argument.

### The commit the plan did not foresee (4b)

cppcheck flagged four new `unreadVariable`/`variableScope` findings in
`chimera.cpp` immediately after phase 4. They were real: unlike the three
other `alignstats` call sites, chimera's `else` branch writes `ci->nw*`
straight from `ci->snw*` and never reads the locals, so once `alignstats`
returned a struct the four locals became pure pass-throughs. Deleting them
also cleared a *pre-existing* `nwcigar` variableScope finding.

This is the case the plan's A1 caveat said would not be available ("do not
expect the locals to vanish in this phase") — it was wrong about chimera
specifically, right about the other three.

### Verification

- **Build:** clean, zero warnings, at every phase.
- **Test suite:** `run_all_tests.sh` → **0 FAIL**.
- **Library:** api_examples `make test` → **40 PASS / 0 FAIL**, an
  identical PASS/FAIL set to the pre-campaign build. Matters here because
  the library session path calls `cluster_query_init` and the chimera
  detection core.
- **Differential:** a purpose-built 1000-query / 40-reference corpus through
  `usearch_global` (alnout, userout with cigar+caln, blast6out, samout, uc,
  fastapairs, matched/notmatched), `cluster_fast` (msaout, consout, profile,
  centroids, uc), `allpairs_global`, `uchime_denovo`, `uchime_ref`
  (uchimealns) and `search_exact` — **byte-identical** to the reference
  binary after every phase, once the command-line banner line is stripped.
- **cppcheck:** no new findings; one pre-existing finding removed. Four new
  `useStlAlgorithm` *suggestions* appeared in `searchcore.cpp` — these are
  false positives that the reference form newly exposed to pattern
  matching. The flagged loops carry `best = &hit`, a side effect
  `std::any_of` cannot express; the honest algorithm would be
  `std::min_element`, which is a loop-modernizing question, not this one.
  Left alone deliberately.
- **clang-tidy:** no new findings in any category, across all 13 touched
  files. `bugprone-easily-swappable-parameters` went **11 → 10**, which is
  phase 4 doing exactly what it was chosen to do.
- **Cross-compiles:** mingw, POWER and mips64el all build, each with a
  warning set identical to the pre-campaign baseline.
- **Performance (the acceptance criterion):** callgrind instruction counts,
  release `-O3`, pre-campaign vs HEAD:

  | command | before | after | delta |
  |---|---|---|---|
  | `usearch_global` | 2,327,568,160 | 2,327,563,222 | −0.0002% |
  | `cluster_fast` | 2,341,866,439 | 2,341,846,088 | −0.0009% |
  | `uchime_denovo` | 3,752,742,965 | 3,752,730,248 | −0.0003% |

  Work-neutral, as predicted — `T &` and `T *` are the same thing at the
  ABI level. Note the `usearch_global` figure is 95% `aligncolumns_*`, which
  dilutes the changed code roughly fourfold; the other two rows are the
  informative ones. No wall-clock benchmark was run, deliberately: at this
  instruction delta any hyperfine reading would be measuring code placement.

### The A2 second comment, resolved afterwards (`a6e1c01f`)

The A2 decision had left `// duplicate: msa.cc` (`utils/cigar.cpp:146`)
untouched. On a follow-up instruction it was deleted rather than reworded,
because there is nothing left to reword *to*: the marker came from
`3172ba45` ("factorize cigar string functions"), which created
`utils/cigar.cpp` and noted that `msa.cc` still held its own copy of the
run-length reader. That copy is gone. `utils/cigar.cpp` is now the only
implementation in the tree — the sole other `std::strtoll` is
`attributes.cpp`'s abundance parser, which is unrelated — and `msa.cpp` is
one of five ordinary callers. It was a task marker whose task is finished.

### Left for the maintainer

- Group D is untouched and still worth doing, D1 (`search16` → `Span`)
  most of all — it is what would let the three remaining `alignstats`
  call sites drop their locals the way chimera's did in 4b.
- Unrelated to this campaign, noticed while checking the above:
  `msa.cpp` walks cigars two different ways — `parse_cigar_string()` at
  `:192`, but a hand-rolled `find_runlength_of_leftmost_operation()` loop
  at `:370-378`. Not a duplicate implementation, just an inconsistent use
  of the two APIs `utils/cigar.hpp` offers. Left alone.


## Decisions (settled 2026-08-08, before any code was written)

1. **A1 — `alignstats` returns an `AlignStats` struct**, not five
   references. References would not fix the swappable-parameter hazard that
   five same-typed out-params create. Phase 4.
2. **A2 — reword `utils/cigar.cpp:132-134` only.** Its closing clause names
   the caller-facing type, which this phase changes. `:146`
   (`// duplicate: msa.cc`) stays untouched. Phase 3.
3. **B8 — dropped.** `xstat`/`xfstat` keep their pointer parameter; the C
   shape deliberately mirrors `stat(2)`. Five `&` call sites stay. Phase 10
   removed from the schedule.
4. **B5 bonus — yes, as a separate commit (7b).** `cluster_query_init` /
   `cluster_query_exit` move into an anonymous namespace *after* their
   signature change, with an independent callgrind run so any delta stays
   attributable.

No open questions remain. Implementation can start at phase 1.
