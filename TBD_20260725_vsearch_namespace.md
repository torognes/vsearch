# Putting the library in `namespace vsearch` (2026-07-25)

Decision record and plan. Follow-on to the library-distribution decision
(source + static-only, C++-only, C façade deferred by design) and to
`TBD_20260725_replace_C_functions.md`, whose `View`/`Span` discussion
raised the question.

**Question:** should every function and type exposed through the library
be declared in `vsearch::`, or is it enough to namespace `Span`/`View`,
which are helpers with an obvious collision risk?

**Answer: all library code goes in `namespace vsearch`.** Namespacing
only the helpers would fix the least dangerous part of the problem.


## Evidence: what the public surface exposes today

Everything reachable from `vsearch_api.h` (which includes `vsearch.hpp`,
`utils/fatal.hpp`, and the ten `core/*.hpp` module headers). Nothing is
in a namespace today except four small internal ones (`fatal_detail`,
`dynlib`, `log_file`, `compression`).

| Name | Kind | Header | Collision risk |
|---|---|---|---|
| `fatal` (3 overloads) | function | `utils/fatal.hpp:105-107` | **severe** — costed out below |
| `dust`, `hardmask`, `dust_all`, `hardmask_all`, `dust_single` | functions | `core/mask.hpp` | **high** |
| `Parameters` | type | `vsearch.hpp:85` | **high** |
| `Database` | type | `core/db.hpp:104` | **high** |
| `View`, `Span` | types | `utils/view.hpp:79`, `utils/span.hpp:85` | high |
| `Dbindex` | type | `core/dbindex.hpp:81` | moderate |
| `Masking` | enum struct | `core/mask.hpp:68` | moderate |
| `DynamicLibraries` | type | `vsearch.hpp:83` | moderate |
| `MergePairs`, `MergeResult`, `MergeInput`, `MergeError`, `QualityTables` | types | `core/mergepairs.hpp` | moderate |
| `seqinfo_s`/`seqinfo_t`, `chimera_result_s`, `cluster_result_s`, `derep_result_s`, `search_result_s`, and the `*_session_s` forward declarations | types | `core/*.hpp` | low (suffixed) |
| `chimera_*`, `cluster_*`, `derep_*`, `search_*` (28 functions) | functions | `core/*.hpp` | low (subsystem-prefixed) |
| `fprint_kmer` | function | `core/dbindex.hpp` | moderate |
| `fprint` (3 overloads), `fprint_integer` (2 overloads), `fprint_spaces` | functions | `utils/print_view.hpp` | **high** — see below |
| `fprint` (3 overloads), `fprint_integer` (2 overloads), `fprint_spaces`, all taking a `Record &` | functions | `utils/print_record.hpp` | **high** — same group |
| `Record`, `OutputRecord` | class template, alias | `utils/print_record.hpp` | moderate |
| `decimal::to_decimal`, `decimal::Buffer`, `decimal::max_width` | function, type, constant | `utils/decimal_digits.hpp` | low (already in `decimal`) |
| `vsearch_api_version`, `vsearch_api_version_string` | functions | `vsearch_api.h` | low (prefixed) |

### Why `fatal()` settles it

`fatal` is a bare, three-overload function named after something almost
any C-ish project has. It is worse than any of the type names, because
**functions collide more dangerously than types**:

- A type clash is a loud, immediate compile error at a line that names
  the type.
- An unqualified `fatal("...")` in consumer code can be silently captured
  by vsearch's overload set through ordinary overload resolution and ADL,
  or become ambiguous in a way that reports against the wrong line.

`dust()`, `hardmask()` and `dust_all()` have exactly the same shape.
Against those, `View` is comparatively benign: it breaks loudly.

### The `fprint*` group

`DONE_20260804_c_style_elimination.md`'s Decision 4 deliberately put five
new global-namespace overloads next to the existing `fprint(View)` rather
than qualifying ~300 call sites mid-migration. `utils/print_record.hpp`,
added late in that migration to batch the hot record writers, then
overloaded the same six names again on a `Record &` sink — so this pass has
**twelve `fprint*` overloads across two headers** to move, plus the `Record`
class template and its `OutputRecord` alias. The two headers must move
together: the overload sets participate in the same name lookup, and
splitting them would make a call resolve differently depending on which
header a translation unit happens to include. That is a bigger but purely
mechanical move, and it **must not be dropped**: these are bare functions
whose first parameter is a `std::FILE *`, which is exactly the shape a
consumer's own output helper has, so they collide the quiet way `fatal()`
does, not the loud way `View` does. The `decimal` namespace is already
there and only needs nesting inside `vsearch`.

So the helper-only option addresses the loud failures and leaves the
quiet ones in place, which is the wrong way round.


## Worked example: `fatal()`, measured

`fatal()` is the largest single name in the surface, so it is the useful
one to cost out. **306 call sites** (comments and string literals
stripped), across:

| Location | Sites |
|---|---|
| `core/` | 104 |
| `cli.cc` | 93 |
| `commands/` | 79 |
| `utils/` | 13 (three of which are the definitions themselves) |
| `os/` | 10 |
| `parameters.cpp` | 6 |
| `arch/` | 1 |

### Almost none of those 306 sites change

This is the central practical point, and it generalises to every name in
the table above. A call such as

```cpp
      fatal("Invalid FASTA - header must start with > character");  // core/fasta.cpp:295
```

is **byte-identical before and after**, because `core/fasta.cpp` is itself
wrapped in `namespace vsearch`: unqualified lookup finds `vsearch::fatal`
in the enclosing namespace. No qualification, no `using`, no edit — at any
of the 104 sites in `core/`, the 79 in `commands/`, the 13 in `utils/` or
the 6 in `parameters.cpp`.

The migration is therefore **one wrap per file, not 306 edits**. That is
the main reason the whole-library scope is affordable, and it is worth
stating explicitly in the review request so the diff's size is not
mistaken for its risk.

### What does change: three boundary cases

1. **`os/` (10 sites)** — only if `os/` stays outside the namespace (open
   question 2 below). Each call then needs `vsearch::fatal(...)`, or one
   `using vsearch::fatal;` per file to leave the call sites untouched.
   This is the cheapest concrete argument for putting `os/` inside.
2. **`main()`** cannot live in a namespace, so `vsearch.cc`'s wrap must
   stop short of it. Its single `fatal()` call is in `flush_stdout()`
   (`:323`), not in `main()`, so the split is clean: wrap the helpers,
   leave `main()` outside, and have it call `vsearch::dispatch_command()`.
3. **Consumer code** — the visible break, e.g.
   `api_examples/example_recover.cc:48`:
   `catch (VsearchError const & error)` becomes
   `catch (vsearch::VsearchError const & error)`, and every example's
   `struct Parameters parameters;` / `Database db;` gains `vsearch::`.
   The ten example files are a good canary for the `LIBRARY_API.md`
   migration note: they are the only in-tree consumer code.

### A trap specific to `fatal()`: three TUs must move in lockstep

`fatal_detail::exit_or_throw` is declared in `utils/fatal.hpp` but
**defined in one of two alternative translation units** —
`fatal_exit.cpp` (CLI, `-fno-exceptions`) or `fatal_throw.cpp` (library,
`-fexceptions`) — selected by the source lists in `Makefile.am`, not by
the preprocessor.

If the header and `fatal.cpp` are wrapped but one of those two TUs is
missed, **there is no compile error**: `fatal.cpp` resolves to
`vsearch::fatal_detail::exit_or_throw` (declared, never defined) while the
missed TU still defines `::fatal_detail::exit_or_throw`. The result is an
undefined reference at *link* time, in whichever configuration links the
missed file — so the CLI can link clean while `libvsearch.a` consumers
fail, or the reverse. Move `fatal.hpp`, `fatal.cpp`, `fatal_exit.cpp` and
`fatal_throw.cpp` in one commit, and build **both** configurations before
pushing.

Any other name whose definition is split across build-selected TUs
deserves the same care.

### Unrelated fix that came out of reading this

Auditing `fatal()` for the above turned up that its three *declarations*
carried no `__attribute__((noreturn))` while all nine definitions
(`fatal.cpp`, `fatal_exit.cpp`, `fatal_throw.cpp`) did, and while
`exit_or_throw`'s declarations in the same header did. Callers therefore
could not know `fatal()` never returns. Fixed separately — it is nothing
to do with namespacing, but it explains why the header now carries three
more attribute lines than the snippets in this document show.

### C++11 constraint on the mechanical wrap

`namespace vsearch::fatal_detail { ... }` is a **C++17** nested namespace
definition and will not compile here. Every nested case needs the
two-level form:

```cpp
namespace vsearch {
namespace fatal_detail {
  // ...
}
}
```

This applies to all four existing sub-namespaces (`fatal_detail`,
`dynlib`, `log_file`, `compression`) and matters because the project must
keep building on GCC 4.9.

### Where the wrap goes in a header

After `#pragma once` **and after the `#include`s** — includes must stay
outside, or the wrap would declare standard-library names inside
`vsearch::`. Existing doc comments stay exactly where they are.


## Why not piecemeal

In order of weight:

1. **A half-namespaced API is worse than either extreme.** Consumers
   cannot predict whether a given name needs `vsearch::`; the
   documentation has to carry an explicit list; and every future public
   name needs a per-name judgment call that will, over time, be made
   inconsistently.
2. **Each partial move is its own source-breaking wave.** Namespacing
   `View` now and `Parameters` later breaks consumers twice. One wave is
   strictly cheaper — and now is the cheapest it will ever be: the
   library is pre-1.0, distribution is source + static (no shared-library
   ABI), and the C façade is deferred, so there is no second boundary to
   keep in step.
3. **"Public only" is not yet an enforceable boundary.** The curated
   `include/` directory does not exist, and `vsearch_api.h` transitively
   pulls in ten module headers. Until that curation lands, "which names
   are public" cannot be answered mechanically, so a public/internal
   namespace split cannot be checked — only remembered.


## Scope

**All library code in `namespace vsearch`**, rather than "everything
public". It is the only rule that is checkable today without first
solving the curation problem, and it costs nothing extra.

- Wrap declarations in the headers and definitions in the `.cpp` files.
  Wide, mechanical, low-risk.
- **Nest the four existing namespaces**: `vsearch::dynlib`,
  `vsearch::log_file`, `vsearch::compression`, `vsearch::fatal_detail`.
  They fit naturally and need no renaming.
- Anonymous namespaces inside `.cpp` files compose fine with an enclosing
  named namespace and stay exactly as they are — that is already the
  established pattern in this codebase.
- **No `using namespace vsearch;` in any header.** There are none in the
  tree today (verified), which is the correct state; adding one to a
  public header would defeat the entire exercise. Inside internal `.cpp`
  files it is a legitimate way to stage a large diff, but wrapping the
  definitions properly is barely more work and leaves no cleanup debt.


## What a namespace cannot do

- **The seven `VSEARCH_API_*` macros are unaffected.** Macros ignore
  namespaces: `VSEARCH_API_VERSION`, `VSEARCH_API_VERSION_MAJOR`,
  `_MINOR`, `_PATCH`, `_STRING`, `VSEARCH_API_STRINGIFY` and
  `VSEARCH_API_STRINGIFY_`. They are already correctly prefixed — leave
  them exactly as they are. This is the one place where the prefix
  convention remains the right tool, and it should be documented as
  deliberate rather than as an oversight.
- **A future C façade must live outside the namespace**: `extern "C"`
  functions with `vsearch_` prefixes, by necessity. Recorded here so the
  decision is not relitigated when C support returns to the table — the
  two do not conflict, they occupy different layers.


## Deliberately kept separate

Once the namespace lands, two names carry a redundant prefix:

- `VsearchError` → `vsearch::Error`
- `vsearch_api_version()` → `vsearch::api_version()`

**Do not bundle these with the namespace move.** A rename is a different
kind of break from a namespace move: the namespace move is mechanically
fixable by consumers (add a qualification or one `using`), whereas a
rename requires them to find and edit each use. Mixing the two makes the
migration note harder to write and consumers' fixes harder to automate.
Do the namespace as a pure, mechanical move first; revisit the redundant
prefixes afterwards, as an optional follow-up that can be declined.


## Suggested implementation order

1. **`utils/` leaves first** (`view.hpp`, `span.hpp`, `fatal.hpp`,
   `seqcmp`, `maps`, `ascii_case`, …). Few dependencies, so the blast
   radius is small and the pattern gets settled before the wide files.
2. **`core/`**, then **`commands/`**, then `cli.cc` / `vsearch.cc`.
3. **`os/`** last: it holds the `x*` wrappers and the `dynlib` backends,
   and it is the layer where `extern "C"` system declarations live, so it
   needs the most care about what must *not* move inside the namespace.
4. `vsearch_api.h` and `LIBRARY_API.md` updated together, including a
   short migration note for consumers and an explicit statement that the
   `VSEARCH_API_*` macros are intentionally un-namespaced.
5. Bump the API version. The consumer-visible break is source-level, so
   this is a major bump by the project's own versioning rules.

Per-step verification, matching the usual checklist: build with
`--enable-debug`, `cppcheck` the touched files, the full test suite
(expect 0 FAIL), the `api_examples` net against a **release** build (a
debug `libvsearch.a` fails there on the `_GLIBCXX_DEBUG` ABI mismatch),
and the Windows / POWER / mips64el cross-compiles.


## Open questions for the maintainer

1. **Namespace name.** `vsearch` is the obvious choice. Worth confirming
   there is no intent to use something shorter (`vs`) — which would
   conflict with the "avoid identifiers shorter than three characters"
   guidance — or a versioned inline namespace, which would be premature
   given there is no ABI to protect.
2. **Should `os/` really be inside?** The argument for is one uniform
   rule. The argument against is that `os/` is a portability shim whose
   job is to wrap platform C APIs, and those declarations must stay
   outside regardless, so the boundary inside those files is already
   mixed. Either answer is defensible; a decision up front avoids churn.
3. **Timing against the `View` migration.** The namespace move should
   land *before* the call-chain `View` migration described in
   `TBD_20260725_replace_C_functions.md` widens, so that the wider use of
   `View` is written qualified from the start rather than being requalified
   later.
