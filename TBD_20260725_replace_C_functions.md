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
calls.


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
| `std::strlen` | 48 | see "The `strlen` root cause" below |
| `std::strcmp` | 13 | sort comparators + `"-"` stdin detection |
| `std::memcpy` | 6 | **keep** — SIMD type-punning |
| `std::memcmp` | 2 | `fastx.cpp` magic-byte sniffing |
| `std::memset` | 1 | `align_simd.cpp` raw SIMD buffer |
| `std::memmove` | 1 | `align_simd.cpp` cigar shift |
| `std::strstr` | 1 | `getseq.cpp` label search |
| `std::strcspn` | 1 | `cli.cc:4053`, truncating optarg at a separator |

**Do not touch the SIMD `memcpy`.** The five in `align_simd.cpp` and the
one in `arch/ppc64le/increment_counters.cpp` implement aliasing-safe
unaligned loads (`std::memcpy(&value, ptr, n)`), which *is* the correct
and portable C++11 spelling for that operation — see the issue #589
comment at `align_simd.cpp:460`. Replacing them with casts would
reintroduce the bug that comment records. `memset`/`memmove` there
operate on raw SIMD scratch buffers for the same reason.

`std::strcmp` splits into two groups:

- **Sort comparators** (10): `db.cpp` ×3, `sortbysize`, `sortbylength`,
  `derep_prefix`, `derep.cpp` ×3. These compare header C-strings inside
  `std::sort` lambdas. They become `View<char>` comparisons once the
  headers are carried as views (`View` already has `operator<`, doing a
  `std::lexicographical_compare`), but not before.
- **`"-"` stdin detection** (2): `fastx.cpp:342`, `open_file.cpp:89`.
  Trivially replaceable, low value.
- **Self-hit label compare** (1): `searchcore.cpp:613`.

### `<cstdio>` non-printf — 36 sites

| Function | Sites |
|---|---|
| `std::sscanf` | 11 |
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

- **CIGAR run-length parsing** (5): `searchcore.cpp` ×3,
  `linmemalign.cpp:733`, `results.cpp:823`. All are
  `sscanf(p, "%" PRId64 "%n", &run, &len)` — read an integer, advance a
  pointer. `utils/cigar.cpp` already does the same job with
  `std::strtoll` and documents the endptr contract, so the natural fix
  is to route these five through the existing cigar helper rather than
  to reimplement the parse.
- **CLI option parsing** (5 in `cli.cc` + 1 helper): `cli.cc:118`, `:153`,
  `:160`, `:230`, `:414`, `:428`. `cli.cc:150` already carries a
  `// refactoring: std::stoi(), faster than sscanf()` note. Careful:
  `:153`/`:160` parse a 3-field and a 2-field comma/`*` grammar
  (`--length_cutoffs`), and several sites rely on `%n` to detect trailing
  garbage — the replacement has to keep rejecting `"12abc"`. `cli.cc:415`
  and `:430` show the existing pattern for that check.

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
| `std::ldiv` | 3 | median computation |
| `std::strtoll` | 2 | already flagged for C++17 |
| `std::free` | 1 | `xfree`, pairs with `posix_memalign` |

- `std::exit` (9): `utils/fatal_exit.cpp` ×3, `utils/fatal_throw.cpp` ×3,
  `cli.cc:4125`/`:4231`, `fastq_mergepairs.cpp:163` (an unreachable guard
  after a `noreturn` `fatal()`). These are the intended process-exit
  points; in a library session `fatal()` throws `VsearchError` instead.
  Not a target.
- `std::ldiv` (3): `sortbylength.cpp:147`, `sortbysize.cpp:143`,
  `derep.cpp:167`, all computing a median index as
  `std::ldiv(size, 2)`. Plain `/` and `%` would do, and would drop the
  `long` round-trip and its narrowing cast. Low risk, low value.
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


## The `strlen` root cause

48 sites, the single largest group, and worth stating precisely because
the obvious diagnosis is wrong.

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
3. **`View` through the output call chains**, `core/results.cpp` first.
   Removes most of the 48 `strlen`, and the 10 comparator `strcmp` fall
   out of the same change.
4. **`sscanf` → the existing cigar helper** (5 CIGAR sites), then the
   `cli.cc` option parsing (6 sites) separately, since that one needs
   care over trailing-garbage rejection.
5. **`ldiv` → `/` and `%`** (3 sites). Cosmetic; bundle with something.

Explicitly **not** targets: the SIMD `memcpy`/`memset`/`memmove`, the
`FILE *`/RAII layer, `std::exit` on the `fatal()` paths, `<ctime>` in
`timestamp.cpp`, `posix_memalign`/`std::free` in `xmalloc`/`xfree`, the
`write(2)` OOM handler, and `assert()`. Blocked on C++11:
`std::from_chars` for the two `strtoll` sites, `std::aligned_alloc`
for `xmalloc`.
