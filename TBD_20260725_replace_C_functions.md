# Remaining C-style functions in vsearch (2026-07-25)

Inventory of the C standard-library functions still called from the
vsearch sources, taken after the `qsort` → `std::sort` migration was
verified complete. The **printf family is excluded** throughout
(`fprintf`/`snprintf`/… are the project's chosen output mechanism and
are not refactoring targets).

Counts are call sites with comments and string literals stripped, over
`src/**/*.{cc,cpp,hpp,h}` excluding `src/vendored/`. `assert()` (85
sites) is excluded as a target: `CLAUDE.md` actively asks for it
("make implicit contracts explicit with assert()").

Total: **144 call sites** (73 + 36 + 17 + 15 + 3), plus ~20 POSIX/OS
calls. **91** after the 2026-07-27 `strlen`, `strcmp`, `sscanf`, `ldiv` and
`strstr` passes (see below).


## Status of `qsort`

Complete. Zero call sites remain anywhere in the repository (sources,
`src/vendored/`, tests, man pages, docs). The seven surviving *mentions*
were stale comments and include justifications, removed or updated in
commit `ed226b7c`; two more of the same kind (`memchr`, an unused
`<cstring>`) in `f1ff2dea`.


## By family

### `<cstring>` — 73 sites

| Function | Sites | Notes |
|---|---|---|
| `std::strlen` | 48 → 22 | see "The `strlen` root cause" below; 26 removed 2026-07-27 |
| `std::strcmp` | 13 → 3 | see "`strcmp` and the signedness trap" below; 10 removed 2026-07-27 |
| `std::memcpy` | 6 | **keep** — SIMD type-punning |
| `std::memcmp` | 2 | **keep** — `fastx.cpp` magic-byte sniffing; see below |
| `std::memset` | 1 | `align_simd.cpp` raw SIMD buffer |
| `std::memmove` | 1 | `align_simd.cpp` cigar shift |
| `std::strstr` | 1 → 0 | **done 2026-07-27** (`f706b6ae`) — it was also a bug |
| `std::strcspn` | 1 | **keep** — `cli.cc`, truncating optarg at a separator; see below |

**Do not touch the SIMD `memcpy`.** The five in `align_simd.cpp` and the
one in `arch/ppc64le/increment_counters.cpp` implement aliasing-safe
unaligned loads (`std::memcpy(&value, ptr, n)`), which *is* the correct
and portable C++11 spelling for that operation — see the issue #589
comment at `align_simd.cpp:460`. Replacing them with casts would
reintroduce the bug that comment records. `memset`/`memmove` there
operate on raw SIMD scratch buffers for the same reason.

### `strcmp` and the signedness trap — DONE (13 → 3)

An earlier revision of this document said the sort comparators "become
`View<char>` comparisons once the headers are carried as views (`View`
already has `operator<`, doing a `std::lexicographical_compare`)".
**That was wrong as written**, and following it would have silently changed
output. It is now true, because `View` was fixed — see "The fix went into
`View`, not into a helper" below.

`std::strcmp` compares bytes as **`unsigned char`**. `View::operator<` is a
`std::lexicographical_compare` over **`char`**, which is signed on x86-64 and
on the Windows target but unsigned on ARM and PowerPC Linux. Any header
carrying a high-bit byte — UTF-8 or Latin-1, which sequence headers routinely
do — would therefore sort differently from the historical output *and*
differently from one architecture to the next. Measured over 2 692 881 pairs
from a byte alphabet straddling `0x80`: `View::operator<` disagrees with
`strcmp` on 1 378 400 of them with signed `char`, and on none with
`-funsigned-char`. Concretely, with tied length and abundance so the header is
the only tie-break:

```
current (strcmp)     View::operator< would give
  >Frad                >Fréd
  >Frzd                >Frad
  >Fréd                >Frzd
```

This is the same class of bug as the `<cctype>` UB recorded below, and it
would also have bitten any future `std::set<View<char>>` (see the TODO at the
top of `span.hpp`).

### The fix went into `View`, not into a helper

The first attempt (`bf5154d8`) added a `utils/header_order.hpp` with
`header_less` and `header_compare`, leaving `View::operator<` alone. That
worked, but it left the trap armed for the next caller, so it was replaced.

The deciding observation: **`std::string` already orders as `unsigned char`**,
through `std::char_traits<char>::lt`. So of the three, `View<char>` was the
only one out of step — with `std::strcmp` *and* with the type it is modelled
on. That is a defect in `View`, not something callers should have had to route
around.

`utils/element_order.hpp` (`d1169e95`) holds the element comparison as a
trait: the primary template defers to the element type's own `operator<`, and
the explicit specialization for `char` reads both bytes as `unsigned char`.
`View::operator<` and `Span::operator<` keep their
`std::lexicographical_compare` and pass the trait as the comparison object, so
only the element step changed. Both classes also gained `compare()` — the
three-way form with `std::strcmp`'s sign convention, for callers that must
tell "equal" from "greater" in one pass — and `operator!=`, which neither had
(C++11 does not synthesise it from `operator==`).

Only `char` is specialized. `signed char` and `unsigned char` say what they
mean and are left to the primary template, the same line `std::char_traits`
draws.

`header_order.hpp` was then deleted (`91cb6032`): its two functions had become
the view's own operators spelled longhand. Perf-neutral — `--sortbysize` over
389 MB measures 1.00 ± 0.01 against `dev`, so passing the trait as a
comparison object costs nothing next to a bare `operator<`.

### Where the 13 sites went

They split three ways, not two:

- **Ordering** (7) — `db.cpp` ×3, `sortbylength:127`, `sortbysize:122`,
  `derep_prefix:327`, `derep.cpp:245`. Now `a < b` and `a.compare(b)` over
  `Database::header_view()`. The three `db.cpp` comparators also lost their
  raw `buffer + header_p` pointer arithmetic to a `header_of` view helper.
- **Header equality** (3) — `derep.cpp:652`, `:671`, `searchcore.cpp:612`.
  `View::operator==` is signedness-independent, and it compares sizes before
  bytes, so it short-circuits where `strcmp` walked to the terminating null.
  Both `derep` sites sit inside hash-collision probe loops and `searchcore`'s
  is in `search_acceptable_aligned`, so all three are on hot paths;
  `--derep_id` measured 4 % faster.
- **`"-"` stdin detection** (3) — `cli.cc:4448`, `fastx.cpp:342`,
  `open_file.cpp:89`. **Keep.** These compare an `argv` pointer against a
  one-character literal, which is `strcmp` used for exactly what it is for;
  a `std::string` temporary would allocate and a hand-rolled
  `s[0] == '-' and s[1] == '\0'` reads worse.

Verified byte-identical against `dev` on a dataset built around the
signedness boundary (`0x80`, `0xA9`, `0xC3`, `0xE9`, `0xFF` beside their
ASCII neighbours, with tied lengths and abundances so the header decides the
order) across `--sortbysize`, `--sortbylength`, `--derep_fulllength`,
`--derep_id`, `--derep_smallmem`, `--derep_prefix`, `--cluster_fast`,
`--cluster_size`, `--uchime_denovo` and `--usearch_global --self/--selfid`:
31/31 outputs identical, both for the helper version and for the
specialization that replaced it.

### `<cstdio>` non-printf — 36 sites

| Function | Sites |
|---|---|
| `std::sscanf` | 11 → 0 |
| `std::fread` | 9 |
| `std::fclose` | 5 |
| `std::fputc` | 3 |
| `std::fflush` | 2 |
| `std::ferror` | 2 |
| `std::fopen`, `std::fgets`, `std::fputs`, `std::fwrite` | 1 each |

The `FILE *` layer is **deliberate and already RAII-wrapped**
(`utils/open_file.hpp` holds `std::unique_ptr<std::FILE, deleter>`;
`open_file.cpp:158` checks `fflush`/`ferror` before `fclose` so a write
error is never silently dropped). This is not a refactoring target: the
file-handling refactor that produced it was completed deliberately, and
`std::fstream` was already evaluated and chosen *against* for the
record-by-record readers.

`std::sscanf` (11) **is** a genuine target, in two clusters:

- ~~**CIGAR run-length parsing** (5)~~ — **done 2026-07-27** (`529efe67`).
  `searchcore.cpp` ×3, `linmemalign.cpp`, `results.cpp`, all routed through
  `utils/cigar.cpp` as planned. One thing the plan did not anticipate:
  `find_runlength_of_leftmost_operation` ends with `std::max(runlength, 1LL)`,
  and that clamp would have made `results.cpp`'s
  `if (run < 0) fatal(...)` unreachable, silently rewriting a negative or zero
  run as 1. The primitive was therefore split — `read_runlength()` returns what
  `strtoll` read, `find_runlength_of_leftmost_operation()` is that plus the
  clamp — and `build_sam_strings` uses the raw form so its guard keeps working.
  Checked against the old idiom on `12M`, `M`, `1M`, `0M`, `-5M`, `7D3I`,
  `999999M`, `""`, `12` and `007M`: identical run and consumed width in every
  case. Neutral on performance end-to-end (the primitive is 7.5× faster, but
  the alignment dominates). Newly covered by 14 checks added to
  `scripts/usearch_global.sh` in vsearch-tests (branch
  `tmp_usearch_global_20260727162751`), which pin the implicit run length of 1
  at both ends of a CIGAR across `caln`, `samout`, `qilo`, `qihi`, `tilo` and
  `alnout` — untested before, and the only reachable observable of the clamp.
  A negative or zero run length is *not* reachable from the CLI, since every
  CIGAR on these paths is produced by the in-tree aligner; the `results.cpp`
  guard stays defensive, and the equivalence table above is its only test.
- ~~**CLI option parsing** (6 in `cli.cc`)~~ — **done 2026-07-27**
  (`299111b8`, `c01d30b0`). Not `std::stoi` as the old note suggested: it
  throws, and exceptions are library-only here, so the replacement is
  `strtoll`/`strtod` with an endptr.

  **A live bug surfaced first** (`299111b8`). `args_getlong` and
  `args_getdouble` tested `sscanf`'s return against `0`, but `sscanf` returns
  `EOF` when the input ends before any conversion. An empty argument was
  therefore read as zero: `--id ""` exited 0 and behaved as `--id 0.0`,
  accepting every hit, and `--topn ""` reported the range message for
  `--topn 0`. Fixed with `ret != 1` as a separate commit, still on `sscanf`,
  so the migration on top stayed behaviour-neutral.

  **The migration also removes undefined behaviour.** An integer that does not
  fit is UB for the scanf family (C11 §7.21.6.2p10), and the `errno == ERANGE`
  the old checks relied on to catch it is a glibc extension — the mingw build
  links a different runtime. `strtoll` clamps and sets `ERANGE` by contract.

  **One deliberate narrowing, agreed with the maintainer:**
  `--length_cutoffs` now takes exactly `shortest,longest,increment` with no
  whitespace. `sscanf` tolerated whitespace wherever a `%d` conversion began,
  so `1,*, 10` parsed while `1, *, 10` did not — a format-string artefact no
  manual page describes. `1,*, 10` is now rejected. Six checks in
  `scripts/fastq_eestats2.sh` pin the grammar.

### `<cctype>` — 17 sites

`std::tolower` 6 · `std::toupper` 5 · `std::isalnum` 5 · `std::isupper` 1.
All 17 now go through `utils/ascii_case.hpp` (`fefea7b4`, `e44277a1`);
`<cctype>` is included in that one header and nowhere else.

Twelve of the seventeen passed an unguarded `char` — fixed in `fefea7b4`;
see the UB section below. (17 calls over 16 lines:
`compare_strings_nocase.cpp:77` holds two.)

### `<cstdlib>` — 15 sites

| Function | Sites | Notes |
|---|---|---|
| `std::exit` | 9 | **keep** — `fatal()` paths |
| `std::ldiv` | 3 → 0 | **done 2026-07-27** (`989a7734`) |
| `std::strtoll` | 2 | already flagged for C++17 |
| `std::free` | 1 | `xfree`, pairs with `posix_memalign` |

- `std::exit` (9): `utils/fatal_exit.cpp` ×3, `utils/fatal_throw.cpp` ×3,
  `cli.cc:4125`/`:4231`, `fastq_mergepairs.cpp:163` (an unreachable guard
  after a `noreturn` `fatal()`). These are the intended process-exit
  points; in a library session `fatal()` throws `VsearchError` instead.
  Not a target.
- ~~`std::ldiv` (3)~~ — **done 2026-07-27** (`989a7734`). All three used only
  `.quot`, the remainder being recomputed with `%` on the next line. The
  `static_cast<long>` was worse than redundant: `long` is 32-bit on the
  `x86_64-w64-mingw32` target (confirmed with a `static_assert` through the
  cross-compiler; 64-bit on Linux), so it narrowed a `std::size_t` / `uint64_t`
  there. Reaching it needs ~2.1 billion records, so it was latent, but plain
  division on the original unsigned type removes the exposure. It also removed
  the second round of casting that `ldiv`'s *signed* quotient forced at every
  subscript — `-Wuseless-cast` now rejects putting those back — and
  `<cstdlib>` from `derep.cpp`. Medians verified identical to `dev` at 0, 1, 2,
  3, 4, 5, 8, 9, 100 and 101 records.
- `std::strtoll` (2): `utils/cigar.cpp:133`, `core/attributes.cpp:213`.
  `attributes.cpp:212` already carries the note that
  `std::from_chars` is the C++17 replacement. **Blocked** — the project
  targets C++11 and must keep building on GCC 4.9.
- `std::free` (1): `os/posix/system.cc:106`, inside `xfree`. It pairs
  with `posix_memalign` in `xmalloc` (`:85`), which POSIX guarantees is
  `free`-able. C++11 has no aligned allocation, so this pair cannot be
  expressed in the standard library before C++17
  (`std::aligned_alloc` / over-aligned `new`). **Keep.**

### `<ctime>` — 3 sites

`std::time`, `std::localtime`, `std::strftime`, all in
`utils/timestamp.cpp:78-80`. **Keep** — as the comment there already
states, `localtime` + `strftime` is the only portable C++11 way to turn
an instant into a calendar string; `<chrono>` gained calendar support in
C++20. The elapsed-time measurement was already moved to
`std::chrono::steady_clock`; only the wall-clock formatting remains.

### POSIX / OS-specific — ~20 sites

Not C-standard, but C-style, and listed for completeness. Well confined:
almost all live in the `src/os/` backends behind the `x*` wrappers
(`xmalloc`, `xfree`, `xfstat`, `xlseek`) and the `dynlib::` loader
family, which is exactly where the OS split put them.

- `src/os/posix/`: `posix_memalign`, `std::free`, `sysconf`, `fstat`,
  `stat`, `lseek`, `ftello`, `dlopen`/`dlsym`/`dlclose`
- `src/os/windows/`: `GlobalMemoryStatusEx`, `GetProcAddress`,
  `FreeLibrary`, `_fseeki64`
- `src/os/{linux,macos,freebsd}/`: `sysconf`, `sysctl`, `sysctlbyname`

Two leaks outside `src/os/`, both minor:

- `fileno()` at four `core/` sites (`udb.cpp:154`, `getseq.cpp:120`,
  `fastx.cpp:319`, `:420`, `:654`) plus `dup`/`close` in `fastx.cpp` and
  `dup`/`fdopen` in `open_file.cpp:117-119`. These are the deliberate
  `FILE *` → fd boundary: `udb` needs `fstat` on the *open descriptor*
  (never on the path — `stat`-on-path is what broke FreeBSD `--db <()`,
  fixed in `fa2e8aac`), and `fastx` needs the fd to hand to zlib.
- `isatty(fileno(stderr))` at `cli.cc:4921`, for progress-bar
  suppression. No standard equivalent exists.
- `write(STDERR_FILENO, ...)` at `vsearch.cc:245`. **Must stay.** It is
  the `new_handler`, which runs with memory exhausted and therefore must
  not allocate; `std::fprintf` may.


## The `strlen` root cause — PARTLY DONE (48 → 22)

48 sites, the single largest group, and worth stating precisely because
the obvious diagnosis is wrong.

**Status (2026-07-27).** Twenty-six of the 48 are gone, in four
commits; see "What was done" at the end of this section for what
remains and why. The diagnosis below held up exactly: not one of the
removed sites needed a new `Database` accessor.

**The lengths are already stored.** `Database` records `headerlen` and
`seqlen` per record at parse time and exposes them as
`getheaderlen()` / `getsequencelen()` (`db.hpp:196-204`), and it already
offers `header_view()`, `sequence_view()`, `quality_view()` and
`mutable_sequence()` returning `View<char>` / `Span<char>`
(`db.hpp:213-243`). 25 call sites already use them.

**Not one of the 48 `strlen` calls is applied to a `Database`
accessor.** Every one is applied to a `char const *` *function
parameter*. The pattern, e.g. `core/results.cpp`:

```cpp
                                 char const * query_head,   // :105
                      query_head,                           // :125
                      static_cast<int>(std::strlen(query_head)),  // :126
```

The length is known at the top of the chain, discarded by passing a bare
pointer, then recomputed at the leaf — an O(n) scan per output line to
recover a number the database already has in a struct field.

So this is **not** a `Database` problem, and adding more accessors will
not fix it. It is a call-chain signature problem: the work is to change
the *parameters* of the output and reporting helpers from
`char const * head` to `View<char> head`, following the same thread that
already converted `searchinfo_s::query_head` to `View<char>`. The
distribution by file shows where to aim:

| File | Sites |
|---|---|
| `core/results.cpp` | 8 |
| `core/getseq.cpp` | 6 |
| `commands/fastq_mergepairs.cpp` | 6 |
| `cli.cc` | 4 |
| `core/fastq.cpp`, `core/cluster.cpp` | 3 each |
| `userfields`, `compare_strings_nocase`, `search`, `msa`, `chimera`, `usearch_global`, `search_exact`, `allpairs_global` | 2 each |
| `vsearch.cc`, `showalign`, `fasta` | 1 each |

`core/results.cpp` is the natural first target: 8 sites, and it is
almost entirely output formatting, so the blast radius is small.

Two of the 48 are not part of this pattern and should be left alone:
`vsearch.cc:207` (summing `argv` lengths — genuinely a C string from the
OS) and `utils/compare_strings_nocase.cpp:114-115` (constructing views
*from* C strings at an API boundary, which is the correct direction).

### What was done (2026-07-27)

Four commits, each verified byte-identical against `dev` on purpose-built
datasets (non-ASCII headers, `--notrunclabels`, `--strand both`,
`--output_no_hits`, header-stripping and relabel options) as well as by the
test suite (0 FAIL), all three cross-compiles, the release build and
`api_examples` (39 PASS).

1. **`results.cpp` and its four callers — 16 sites.** The `results_show_*`
   helpers now take `View<char>` for `query_head`, `qsequence` and
   `qsequence_rc`. The separate `qseqlen` argument disappears from every
   helper that also receives the sequence (`alnout`, `userout`, `qsegout`);
   `blast6out` and `uc` keep theirs, having no sequence to derive it from.
   The header is printed with `"%.*s"` instead of `"%s"` — equivalent,
   because every caller's view comes from `Database::header_view()` or from
   `searchinfo_s::query_head_v` (filled with `head_len + 1` bytes and viewed
   at `head_len`). `build_sam_strings` took the same treatment, which also
   turned its `qpos`/`tpos` into `std::size_t`.

   The eight `fasta_print_general` header lengths in `usearch_global`,
   `allpairs_global`, `search_exact` and `cluster` fell out at the same
   time: those four callers already held a `View<char> query_head_view` and
   immediately called `.data()` on it.

2. **`fastq_mergepairs` — 6 sites.** Not a signature problem: the two header
   buffers in `merge_data_s` are `std::vector<char>` sized to the longest
   header seen so far, so `size()` says nothing about the current pair.
   `read_pair()` already reads both lengths to size those buffers; they are
   now stored as `fwd_header_length`/`rev_header_length`, matching the
   existing `fwd_length`/`rev_length` sequence fields.

3. **`fastq_print` deleted — 2 sites.** It had no caller anywhere in the
   tree and `core/fastq.hpp` is not among the module headers `vsearch_api.h`
   pulls in, so it was not reachable by a library user either. Its FASTA
   counterpart `fasta_print` does have a caller (`msa.cpp:503`) and stays.

4. **`msa` — 2 sites.** `msa_target_s::cigar` became a `View<char>` borrowed
   from the `clusterinfo` `std::string` whose `size()` the caller has. This
   needed `find_runlength_of_leftmost_operation`'s endptr out-parameter to
   become `char const **`; `std::strtoll` still wants a `char **`, so the
   wrapper keeps that in a local and the laundering stays in the one place
   that talks to the C API.

Performance (release, hyperfine vs a `dev` release build): neutral to
slightly better — `usearch_global` with `--alnout/--blast6out/--uc/--userout`
and with `--samout` both within noise, `cluster_size --msaout` ~3 % faster,
`fastq_mergepairs` ~3 % faster.

### The `strstr` in `getseq.cpp` was a bug — DONE

`std::strstr()` never returns `nullptr` for an empty needle, so with
`--label_word ""` the loop in `test_label_match()` had no exit other than a
successful boundary check: it walked past the end of the header until it found
two adjacent non-alphanumeric bytes, which the zero-filled tail of the 8 KB
`FastxBuffer` always supplies. An empty `--label_word` therefore matched every
record. The read stays inside the buffer's capacity, so ASan does not flag it —
a logic bug, not a memory-safety one.

The `--label_words` branch was already immune, using a `std::search` bounded by
`header_end`; a human comment there records why. `f706b6ae` gives
`--label_word` the same bound (keeping the end-of-header position, which an
empty needle may legitimately match), and `8a0aee28` then collapses the two
now-identical loops into one `matches_delimited()` helper — which is where the
four `*(hit - 1)` / `*(hit + wlen)` sites the `<cctype>` pass deferred finally
go.

Cross-checked against usearch 9.2.64, 10.0.240 and 11.0.667: they also accept
an empty `-label_word` and apply the ordinary rule to it.

**The delimiter class differs from usearch, deliberately — settled, do not
"fix" it.** usearch delimits words with *anything that is not a letter*, so
digits delimit and `-label_word sample` matches `sample1`. vsearch delimits
with *anything that is not alphanumeric*, so `sample` is not a word of
`sample1`. Confirmed with the maintainer on 2026-07-27: vsearch's definition is
intended, and it is the one `vsearch.1:1561` documents. The difference is
visible on abundance-annotated headers, which is the common case, so it is the
sort of thing that looks like a bug to anyone diffing the two tools — hence
this note.

Nothing pinned that choice, so three checks were added to
`scripts/fastx_getseqs.sh`: a digit does not delimit (`>sample1` +
`--label_word sample` → no match), neither does a letter (`+ "1"` → no match),
and the whole alphanumeric run does match (`+ "sample1"`). They hold against
v2.31.0 too. Changing the delimiter class now means changing the manual page
and those tests first, which is the intent.

## Two option-validation bugs found while probing the CLI parsers

Neither is a C-function question, but both surfaced from the `sscanf` work and
are recorded here so the reasoning is not lost.

### `--id ""` was read as `--id 0.0` — FIXED (`299111b8`)

`args_getlong`/`args_getdouble` compared `std::sscanf`'s return against `0`,
but `sscanf` returns `EOF` when the input ends before any conversion. An empty
argument therefore left the output variable at its initialiser and the `"%n"`
width at `0`, so the trailing-garbage check (`0 < strlen("")`) did not fire
either. `--id ""` exited 0 accepting every hit; `--topn ""` reported the range
message for `--topn 0`. `ret != 1` is the correct test.

### `--topn -1` silently ignored the option — FIXED (`8b14e45c`)

The manual page specifies a *positive* integer and `0` was rejected, but the
guard tested `== 0`, so negatives passed and then read as **"no limit"** three
separate accidental ways: `deck.size() > static_cast<unsigned long>(-1)` is
false (`sortbysize.cpp:192`), `std::min(selected, static_cast<uint64_t>(-1))`
is `selected` (`derep.cpp:155`), and `selected == opt_topn` never holds for a
rising count (`derep_prefix.cpp:402`). `< 1` is the correct test, and it was
the only integer option in `cli.cc` guarded with `== 0` — the others use
`<= 0` or `< 0`.

usearch accepts `0` and negatives alike as "no limit" and never errors, so this
widens a divergence the manual page already justifies rather than creating one.

**Four tests in vsearch-tests asserted the opposite of their own names**
(`success`/`failure` the wrong way round after `--topn "-1"`, in
`derep_fulllength.sh`, `derep_id.sh`, `derep_prefix.sh`, `fastx_uniques.sh`)
and were green only because of the bug. Corrected in the same round.

### `--gapopen ""` is *not* a bug

Worth recording because it looks like one. The manual page sets out an
initialize-then-override model — "vsearch always initializes the six gap
opening penalties using the default parameters (20I/2E). The user is then free
to declare only the values he/she wants to modify" — so an empty declaration
modifies nothing. Measured: identical raw score to no option at all and to an
explicit `20I/2E` (50, against 66 for `4I/2E`). usearch accepts an empty
`-gapopen` the same way while rejecting a malformed one. Pinned by five checks
in `scripts/usearch_global.sh` rather than changed.


### The 22 that remain, and why

| Sites | Where | Why it stays |
|---|---|---|
| 4 | `cli.cc:155`, `:162`, `:415`, `:430` | `argv` strings; trailing-garbage checks. Belongs with the `sscanf` item below, not here |
| 5 | `getseq.cpp` ×5 | **done 2026-07-27** (`f706b6ae`, `8a0aee28`): `strstr` gone, the pointer-arithmetic pass done with it. The five left all build a length from an `argv` string or an `fgets` buffer at the boundary — the right direction |
| 2 | `search.cpp:258`, `:413` | library API boundary — building a length *from* a user's C string, the correct direction |
| 2 | `chimera.cpp:2855`, `:2863` | same, plus `:2855` is a deliberate length-vs-`strlen` consistency check on caller-supplied data (S18) |
| 2 | `compare_strings_nocase.cpp:112-113` | already flagged above as the correct direction |
| 2 | `fasta.cpp:490`, `fastq.cpp:653` | `opt_label_suffix`, an `argv` string |
| 1 | `cluster.cpp:1307` | `opt_clusters`, an `argv` string |
| 2 | `userfields.cpp:124`, `:149` | a table of string literals and a parse cursor over one |
| 1 | `vsearch.cc:207` | `argv`, already flagged above |

So of the 48, the "length known at the top of the chain, discarded, then
recomputed at the leaf" class is now empty except for `getseq.cpp`, which
is deliberately deferred to the pointer-arithmetic pass over that same
function. The rest are C strings arriving from `argv` or from a library
caller, where computing the length once at the boundary is the right thing.

Not attempted here: the 10 comparator `strcmp`, which the ordering below
expects to fall out of this change. They did not — those comparators sort
`Database` headers inside `std::sort` lambdas and are a separate thread
from the output call chains.


## Does exposing `View`/`Span` leak internal types onto library users?

Short answer: **yes, and it already does today** — but the leak is
bounded, and users are not forced to use them.

### It is already in the public header

`vsearch_api.h:165` includes `core/db.hpp`, and `db.hpp:65-66` includes
`utils/span.hpp` and `utils/view.hpp`. Any consumer of the library API
already gets `View` and `Span` declared. Converting more call chains to
views does not open a new hole; it widens use of one that is already
open.

### Users are not forced

The raw accessors were kept *alongside* the view ones:
`getheader()`/`getheaderlen()` still exist next to `header_view()`. The
view accessors are purely additive and opt-in. As long as the pointer
accessors are not deleted, a library user can ignore `View` entirely.
**This is the property to preserve** — it should be an explicit rule for
the migration, not an accident.

Interop is also one line in either direction, since `View` is a
pointer+length pair with `data()`/`size()`:

```cpp
auto const header = db.header_view(seqno);
std::string owned(header.begin(), header.end());   // to a std::string
some_c_api(header.data(), header.size());          // to a C API
```

### The real hazard is the global namespace, not the type

`View` and `Span` are declared with **no namespace at all**
(`view.hpp:79`, `span.hpp:85`). A library user with their own `View` or
`Span` class — hardly exotic names — gets an ambiguity or a hard
collision. That is the actual problem, and it is orthogonal to how many
of our own call sites use them.

This is what the already-agreed `vsearch::` namespace does for. It
should land **before** the view migration widens, not after: moving
`View`/`Span` into `vsearch::` later is a source-breaking change for
every consumer that has started naming the type.

### Why the usual objection does not apply here

Exposing a header-only class template in a public header is normally an
ABI hazard. Here it is not, for two reasons already settled:

- Distribution is **source + static-only, C++-only**, with the C façade
  deferred by design. There is no shared-library ABI to keep stable, so
  the template's layout is not a compatibility surface.
- `View`/`Span` are dependency-free and self-contained, including only
  `<algorithm>`, `<cassert>`, `<cstddef>`, `<iterator>` and `<limits>`,
  and they mirror a documented subset of C++20 `std::span` — so a user
  reading the header learns nothing surprising, and the eventual C++20
  migration is a typedef.

### Recommendation

1. Put `View`/`Span` in `namespace vsearch` **first**. Highest value,
   and it gets more expensive the longer it waits.
2. Then widen the view migration through the call chains, starting with
   `core/results.cpp`.
3. Keep the raw pointer + `*len()` accessors on `Database`
   indefinitely, and say so in `LIBRARY_API.md`, so the view types stay
   opt-in.
4. Follow the existing precedent for genuinely internal types: the
   comment at `vsearch_api.h:158-163` documents why `searchcore.hpp` is
   deliberately *not* in the public header (`searchinfo_s`, `struct hit`).
   `View`/`Span` are the opposite case — vocabulary types, intended to
   cross the boundary — and `LIBRARY_API.md` should state that
   distinction rather than leave it implied.


## Latent UB: unguarded `<cctype>` calls — FIXED (`fefea7b4`)

All 12 sites below were routed through a new `utils/ascii_case.hpp`
(`to_upper` / `to_lower` / `is_alnum` / `is_upper`) in commit
`fefea7b4`. The six files each dropped their now-unused `<cctype>`.
Tests 0 FAIL; Windows, POWER and mips64el cross-compiles clean. The
section is kept for the rationale and for the five sites deliberately
left alone.

**12 of the 17 `<cctype>` call sites passed a plain `char`.**

`std::toupper`, `std::tolower`, `std::isalnum` and friends take an `int`
whose value must be representable as `unsigned char`, or be `EOF`.
Where `char` is signed (the default on x86-64 and on the Windows and
POWER targets; note ARM and PowerPC Linux default `char` to *unsigned*,
so this is also an architecture-dependent divergence), a byte with the
high bit set converts to a **negative** `int`, and the call is
**undefined behaviour**.

This matters for vsearch because sequence *headers* are arbitrary bytes:
a FASTA description carrying UTF-8 or Latin-1 — an accented author name,
a non-ASCII locality — is enough to reach it. Sequence data is normally
ASCII, so the header paths are the exposed ones.

In practice glibc's `tolower` indexes a table offset by 128 and tolerates
`-128..255`, which is why this has never been observed to misbehave. It
is still formally UB, it is not portable to other libcs, and it is
exactly the class of thing a hardened or sanitizing build flags.

### Unguarded sites

| Location | Call |
|---|---|
| `core/mask.cpp:160` | `std::toupper(seq[i])` |
| `core/msa.cpp:107` | `std::toupper(nucleotide)` — param is `char const` |
| `core/tax.cpp:157` | `std::tolower(header.data()[offset])` — header path |
| `core/chimera.cpp:1011` | `std::tolower(ci->paln[f][i])` |
| `core/chimera.cpp:1298` | `std::tolower(ci->paln[0][i])` |
| `core/chimera.cpp:1303` | `std::tolower(ci->paln[1][i])` |
| `core/chimera.cpp:1539` | `std::tolower(ci->diffs[i])` |
| `core/getseq.cpp:270` | `std::isalnum(*(hit - 1))` |
| `core/getseq.cpp:272` | `std::isalnum(*(hit + wlen))` |
| `core/getseq.cpp:319` | `std::isalnum(*(hit - 1))` |
| `core/getseq.cpp:321` | `std::isalnum(*(hit + wlen))` |
| `commands/fastx_mask.cpp:132` | `std::isupper(seq[j])` |

`core/getseq.cpp` is the most exposed of these: it classifies bytes
*adjacent to* a match inside a header (`--label_word` boundary
checking), so it is applied to arbitrary header content by construction.
It also does pointer arithmetic (`*(hit - 1)`, `*(hit + wlen)`) that
`CLAUDE.md` asks to replace with `std::next`, so the two cleanups belong
in one pass.

### The five that were already correct — since converted too (`e44277a1`)

These were never UB; they already cast through `unsigned char` inline.
They were converted anyway so that one idiom remains in the tree rather
than two, and `<cctype>` is now included in exactly one place
(`utils/ascii_case.hpp`).

| Location | Note |
|---|---|
| `commands/sff_convert.cpp:584`, `:588` | inline `static_cast<unsigned char>` |
| `core/otutable.cpp:196` | inline cast inside the lambda |
| `utils/compare_strings_nocase.cpp:73-77` | casts into named `lhs_unsigned`/`rhs_unsigned` (`:75-76`), plus an `assert((lhs >= 0) or (lhs == EOF))` at `:74` |

`compare_strings_nocase.cpp` was the best model: it makes the contract
explicit with an `assert` as well as casting. That `assert` is kept in
the converted form.

### The fix that was applied

A shared helper, in preference to twelve inline casts — one place to
state the contract, and it reads better at each call site. Landed as
`src/utils/ascii_case.hpp`, registered in `src/Makefile.am`:

```cpp
inline auto to_upper(char const character) -> char {
  return static_cast<char>(std::toupper(static_cast<unsigned char>(character)));
}
inline auto to_lower(char const character) -> char {
  return static_cast<char>(std::tolower(static_cast<unsigned char>(character)));
}
inline auto is_alnum(char const character) -> bool {
  return std::isalnum(static_cast<unsigned char>(character)) != 0;
}
inline auto is_upper(char const character) -> bool {
  return std::isupper(static_cast<unsigned char>(character)) != 0;
}
```

Returning `bool` from the `is_*` wrappers also removed the `!= 0` /
`== 0` comparisons that the `isalnum` and `isupper` call sites needed in
order to satisfy the "avoid implicit bool conversion" rule.

Not marked `noexcept`: the `<cctype>` functions are not declared
`noexcept` by the standard, so the wrappers promise no more than what
they call, and the comparable small helpers in `utils/` (`seqcmp`,
`round_up`, `string_normalize`) carry no `noexcept` either.

Note on scope: these functions are **locale-dependent**. vsearch never
calls `setlocale`, so it runs in the `"C"` locale and the behaviour is
plain ASCII case folding. If pure-ASCII behaviour is what is wanted —
and for nucleotide data it is — a locale-free implementation would be
more honest still, but that is a separate decision and was deliberately
not bundled with the UB fix.

### Still open at these sites

- The four `core/getseq.cpp` sites keep their pointer arithmetic
  (`*(hit - 1)`, `*(hit + wlen)`), which `CLAUDE.md` asks to replace with
  `std::next` / `std::prev`. Left untouched to keep the UB fix reviewable
  on its own; worth doing in a follow-up pass over that function.
- ~~The five already-guarded sites still spell the cast inline.~~
  **Done** in `e44277a1`: all 17 sites now go through the helper, and
  `<cctype>` appears in exactly one file.


## Suggested ordering

1. ~~**`<cctype>` UB fix**~~ — **done**, `fefea7b4`.
2. **`vsearch::` namespace for `View`/`Span`.** Cheapest now, and it
   gates item 3.
3. ~~**`View` through the output call chains**, `core/results.cpp`
   first.~~ — **done for the output chains** (48 → 22 `strlen`, after
   `fastq_print` was deleted as dead code); see
   "What was done" above. Note this was taken **before** item 2, at the
   maintainer's request. The gating argument in item 2 is about
   *external* consumers naming the type, and every site converted here
   is internal, so the later move into `vsearch::` is unaffected — but
   `View`/`Span` are now named in `results.hpp`, `msa.hpp` and the new
   `utils/element_order.hpp` as well, so item 2 has grown by three
   headers. The 10 comparator `strcmp` did **not** fall out of this
   change; they were done separately (13 → 3), and the way the plan
   assumed would have changed output — see the signedness section.
4. ~~**`sscanf` → the existing cigar helper** (5 CIGAR sites)~~ — **done**,
   `529efe67`; see the `<cstdio>` section for the clamp that had to be split
   out first. The `cli.cc` option parsing (6 sites) is still open and still
   wants its own pass: several rely on `%n` to reject trailing garbage, and
   `std::stoll`/`std::stod` throw, so the replacement is `strtoll` with an
   endptr rather than the `sto*` family.
5. ~~**`ldiv` → `/` and `%`** (3 sites)~~ — **done 2026-07-27** (`989a7734`).
   Not purely cosmetic after all: `static_cast<long>` narrowed a `std::size_t`
   on the Windows target, where `long` is 32-bit (confirmed with a
   `static_assert` through the cross-compiler). Removing `ldiv` removed the
   cast, the second round of casting its signed quotient forced at each
   subscript, and `<cstdlib>` from `derep.cpp`.

Explicitly **not** targets: the SIMD `memcpy`/`memset`/`memmove`, the
`FILE *`/RAII layer, `std::exit` on the `fatal()` paths, `<ctime>` in
`timestamp.cpp`, `posix_memalign`/`std::free` in `xmalloc`/`xfree`, the
`write(2)` OOM handler, and `assert()`.

Two more joined that list on inspection (2026-07-27), both of which had looked
like one-line cleanups:

- **the two `memcmp` in `fastx.cpp`.** Unifying them on the `std::equal` used
  higher up in the same file *breaks gzip detection*: there the buffer is a
  `char const *` while `magic_gzip` is `std::array<unsigned char, 2>`, so
  `std::equal` compares `0x8b` as 139 against the same byte as −117 wherever
  `char` is signed. Measured: `memcmp` matches, `std::equal` does not,
  `std::equal` through a `reinterpret_cast<unsigned char const *>` does — and
  that cast reads worse than the `memcmp`. Same trap as the `strcmp`
  comparators. Only the repeated literal `2` was replaced, by
  `magic_gzip.size()` (`db44dc30`).
- **the `strcspn` in `cli.cc`.** `optarg[std::strcspn(optarg, "; \t\r\n\v\f")] = '\0'`
  is one call. `std::find_first_of` needs an added `std::strlen` for the end
  iterator and an off-by-one to keep the separator array's own NUL out of the
  separator set. Strictly worse; kept. Blocked on C++11:
`std::from_chars` for the two `strtoll` sites, `std::aligned_alloc`
for `xmalloc`.
