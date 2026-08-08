# TBD 2026-08-08 — msa.cpp uses both cigar-walking APIs, inconsistently

Noticed while removing the stale `// duplicate: msa.cc` marker from
`utils/cigar.cpp` (commit `a6e1c01f`). Not part of that work, and not a
live bug — see "Severity" below before scheduling this.


## The finding

`utils/cigar.hpp` offers two ways to walk a cigar string:

```cpp
auto find_runlength_of_leftmost_operation(char const * first_character,
                                          char const * & first_non_digit) -> long long;

auto parse_cigar_string(View<char> cigar_string) -> std::vector<std::pair<Operation, long long>>;
```

`core/msa.cpp` uses **both**, on **the same cigars**, in two functions that
`msa()` calls one after the other:

| | walk A — `find_max_insertions_per_position` (`msa.cpp:181`) | walk B — `compute_and_print_msa` (`msa.cpp:337`) |
|---|---|---|
| API | `parse_cigar_string()` (`:188`) | hand-rolled `find_runlength_of_leftmost_operation()` loop (`:370-378`) |
| switches on | `Operation` enum, three cases, no `default` | raw `char` (`'D'`/`'M'`/`'I'`), plus `default: break;` |
| allocates | one `std::vector<std::pair<...>>` per target | nothing |
| `"12M1"` (trailing digits, no operation letter) | `fatal("ill-formed CIGAR string")` — `cigar.cpp:169-172` | **no guard**; reads `*next_operation` with `next_operation == cigar_end` |
| unknown letter (e.g. `'X'`) | `assert` fires in debug; in release `convert_to_operation` falls through to `return Operation::match` (`cigar.cpp:93-102`) | silently skipped by `default: break;`, so `qpos`/`tpos` never advance |

Both functions iterate `targets.drop(1)` over the same `msa_target_s::cigar`
views. `msa()` calls A at `:595` and B at `:609`, unconditionally, with no
early return between them. So **every cigar in every cluster is parsed
twice, two different ways, and the two ways disagree about malformed
input.**

The disagreement is currently unreachable: A runs first and has already
fataled (trailing digits) or normalised (unknown letter) the same bytes
before B sees them. B's handling of both cases is therefore dead code whose
deadness depends on an unstated precondition — *"A ran first, over this
same data"* — that nothing in B's signature, comments or call site records.

Both helpers are file-local (`nm -C src/core/msa.o` shows lowercase `t`;
anonymous namespace spans `msa.cpp:180-578`), so Option 2 is contained
entirely to `msa.cpp`. Option 3 also touches `utils/cigar.{hpp,cpp}`.


## Severity — read this before scheduling

**This is a maintainability fix, not a bug fix.** The cigars reaching
`msa()` are produced by vsearch's own aligners (`hit.nwalignment`, from the
SIMD aligner or `LinearMemoryAligner`), and both emit only `M`/`I`/`D` with
well-formed run lengths. Neither malformed case is reachable from user
input today. Nothing is miscomputed, and no output is wrong.

What is worth fixing is the coupling: B is correct only because A ran
first, and that is invisible at B's call site. A future reordering, an
early return added to `msa()`, or a reuse of `compute_and_print_msa()`
elsewhere would silently move B's divergent paths from unreachable to
reachable.

Option 1 from the original discussion (point walk B at
`parse_cigar_string`, change nothing else) is deliberately **not** written
up here: it fixes the inconsistency but takes allocation from 1/target to
2/target, making the code strictly more expensive for no other gain. Both
options below are cheaper than today.


## Option 2 — parse each cigar once, hand the result to both walks

*Recommended if this gets scheduled at all.* Removes the inconsistency and
the redundant second parse together; the code ends up doing **less** work
than today (1 parse/target instead of 2).

### Shape

`msa()` parses once, up front, and both helpers become consumers:

```cpp
// core/msa.cpp, in msa()

using ParsedCigar = std::vector<std::pair<Operation, long long>>;

/* one entry per target, indexed by target position so the indices line up
   with `targets` exactly; index 0 is the centroid, which carries no cigar
   and stays empty */
std::vector<ParsedCigar> parsed_cigars(targets.size());
for (std::size_t i = 1; i < targets.size(); ++i) {
  parsed_cigars[i] = parse_cigar_string(targets[i].cigar);
}

auto const max_insertions = find_max_insertions_per_position(parsed_cigars, centroid_length);
...
compute_and_print_msa(targets, parsed_cigars, max_insertions,
                      profile, aln_v, fp_msaout, db, parameters);
```

**Size the vector to `targets.size()`, not `targets.size() - 1`.** Indexing
a `drop(1)`-shaped array alongside a full-length one is a parallel-array
off-by-one waiting to happen; leaving index 0 empty costs one unused
`std::vector` (24 bytes) per cluster and removes the skew entirely.

Walk A then loses its own parsing and its outer loop shrinks to the
consumption it always was:

```cpp
auto find_max_insertions_per_position(std::vector<ParsedCigar> const & parsed_cigars,
                                      int const centroid_len) -> std::vector<int> {
  std::vector<int> max_insertions(static_cast<std::vector<int>::size_type>(centroid_len + 1));

  // index 0 is the centroid, which carries no cigar string
  for (std::size_t i = 1; i < parsed_cigars.size(); ++i) {
    auto position = 0LL;
    for (auto const & a_pair : parsed_cigars[i]) {
      // ... unchanged body ...
    }
  }
  return max_insertions;
}
```

Walk B drops its hand-rolled loop and switches on `Operation`:

```cpp
// replaces msa.cpp:368-378
for (auto const & a_pair : parsed_cigars[target_index]) {
  auto const operation = a_pair.first;
  auto const runlength = a_pair.second;

  switch (operation) {
  case Operation::deletion:   // was case 'D'
    ...
    break;
  case Operation::match:      // was case 'M'
    ...
    break;
  case Operation::insertion:  // was case 'I'
    ...
    break;
  }
}
```

The `default: break;` goes away: the switch over a three-value scoped enum
is exhaustive, which is also why walk A never needed one. (CLAUDE.md's
"avoid switch missing default case" is about switches on open types; GCC
does not warn on an exhaustive enum switch, and walk A is the precedent
already in this file.)

Note walk B's loop currently iterates `targets.drop(1)` by reference and
has no index, so it needs one to reach `parsed_cigars`. Either enumerate
with an index, or — cleaner — bundle the parsed cigar into the iteration.
Worth deciding at implementation time; the index is the boring choice.

### Cost

Peak memory rises from one parsed cigar at a time to all of them at once:
roughly `total_operations × 16` bytes per cluster (each
`pair<Operation, long long>` is 16 bytes after padding — `char` + 7 padding
+ `int64_t`). A 10 000-member cluster averaging 100 operations per cigar is
about 16 MB. Bounded, predictable, and freed when `msa()` returns, but it
is a genuine trade against today's transient single-cigar footprint.

If that ceiling is unwelcome, Option 3 removes it rather than trading it.

### Verification

`--msaout`, `--consout` and `--profile` output must be **byte-identical**.
The differential harness from the `&parameter` campaign already covers
exactly this (`cluster_fast` with all three outputs, `--threads 1`); reuse
it. Add a callgrind instruction count on a `--msaout` run to confirm the
removed second parse actually shows up — it should be a small but real
*reduction*, and if it is not, that is worth understanding before merging.


## Option 3 — a non-allocating visitor over cigar operations

The endpoint. Removes the inconsistency **and** all allocation in both
walks, and leaves `utils/cigar.hpp` with one walking implementation instead
of two. Larger change; worth it only if the per-target allocation ever
shows up in a profile.

### What "visitor" means here

The two APIs differ in *who drives the loop*.

`parse_cigar_string` is **pull** iteration: it materialises the entire
sequence of `(operation, runlength)` pairs into a `std::vector`, hands the
container back, and the caller pulls elements out with its own loop. The
vector exists only so the caller has something to iterate — it is a
temporary buffer between two loops that could have been one.

A visitor (also called *internal iteration*, or a *push* walker) inverts
that: the walker keeps the loop, and the caller supplies a callable that
the walker invokes once per operation. Nothing is stored. Each
`(operation, runlength)` pair lives in registers or on the stack for the
duration of one call and is gone. The caller still writes what looks like a
loop body; it just writes it as a lambda instead.

The trade is control flow. Inside a normal `for` loop a caller can `break`,
`continue`, or `return` from the enclosing function. Inside a callback it
can do none of those — `return` exits only the callback. Walkers that need
early exit usually have the callback return `bool` ("keep going?"). **Both
msa.cpp callers consume every operation to the end**, so the simple `void`
form is sufficient here and should be preferred; do not add the `bool`
protocol speculatively.

### The declaration

```cpp
// utils/cigar.hpp

/* Walks a cigar string and invokes 'visit(operation, runlength)' once per
   operation, in order, without allocating -- the push counterpart to
   parse_cigar_string(), which materialises the same sequence into a vector
   for callers that need to hold it. Prefer this one for a single ordered
   pass; the callback cannot break out of the walk. Rejects an ill-formed
   cigar (trailing digits with no operation letter) exactly as
   parse_cigar_string() does, because that one is now written in terms of
   this. */
template <typename Visitor>
auto for_each_cigar_operation(View<char> const cigar_string, Visitor && visit) -> void {
  auto const * position = cigar_string.begin();
  auto const * const cigar_end = cigar_string.end();

  while (position < cigar_end)
    {
      // Consume digits (if any), return the position of the first char
      // (M, D, or I), store it, move cursor to the next byte.
      // next_operation aliases the read-only cigar buffer and is only read.
      char const * next_operation = nullptr;
      auto const run = find_runlength_of_leftmost_operation(position, next_operation);
      // do not dereference if outside of cigar_string! (= missing operation!)
      if (next_operation >= cigar_end) {
        // fail if ill-formed (ex: '12M1'), we could also silently skip over
        fatal("ill-formed CIGAR string");
      }
      // operations: match (M), insertion (I), or deletion (D)
      auto const operation = *next_operation;
      position = std::next(next_operation);
      visit(convert_to_operation(operation), run);
    }
}
```

It has to be a template, so it has to live in the header, and that has
consequences worth stating plainly:

- **`convert_to_operation` must move out of `cigar.cpp`'s anonymous
  namespace** (`cigar.cpp:91-102`) and into `cigar.hpp` as an `inline`
  function. It is currently deliberately TU-local — that comment
  (`// anonymous namespace: limit visibility and usage to this translation
  unit`) is a decision being reversed, so it needs sign-off, and the
  comment needs updating rather than silently deleting.
- **`cigar.hpp` gains includes**: `"fatal.hpp"`, `<cassert>` (for
  `convert_to_operation`'s assert) and `<iterator>` (`std::next`). `cigar.cpp`
  already includes all three, so nothing new enters the build — but a
  header pulling in `fatal.hpp` widens what every cigar consumer sees.
- **C++11 has no generic lambdas.** Every call site must spell the parameter
  types: `[&](Operation const operation, long long const runlength) { ... }`,
  not `[&](auto op, auto run)`. This is the single most likely thing to trip
  up an implementer used to C++14.

### `parse_cigar_string` becomes a caller

This is what makes the option worth doing: one walking implementation, and
the pull API keeps working for anyone who genuinely needs to hold the
result.

```cpp
// utils/cigar.cpp

auto parse_cigar_string(View<char> const cigar_string) -> std::vector<std::pair<Operation, long long>> {
  std::vector<std::pair<Operation, long long>> parsed_cigar;
  for_each_cigar_operation(cigar_string,
                           [&parsed_cigar](Operation const operation, long long const runlength) -> void {
                             parsed_cigar.emplace_back(operation, runlength);
                           });
  return parsed_cigar;
}
```

### Every `parse_cigar_string` caller in the tree is single-pass

Checked while writing this up, and it is the strongest argument for Option
3. There are seven call sites, and **not one of them holds the returned
vector for anything**; every single one is `parse → range-for once →
discard`:

| site | shape |
|---|---|
| `utils/cigar.cpp:186` (`print_uncompressed_cigar`) | parse → loop → print each letter |
| `core/msa.cpp:188` (walk A) | parse → loop → accumulate `max_insertions` |
| `core/showalign.cpp:167` | parse → loop inline |
| `core/showalign.cpp:386` | parse → loop → `putop()` per operation |
| `core/chimera.cpp:451` | parse → loop → switch on `Operation` |
| `core/chimera.cpp:846` | parse → loop → switch on `Operation` |
| `core/chimera.cpp:884` | parse → loop → switch on `Operation` |

So the vector `parse_cigar_string` returns is, at every existing call site,
a temporary buffer between two loops that could have been one. That is not
an argument for deleting the pull API — a caller that needs random access
or a second pass may appear — but it does mean Option 3's visitor is the
shape six of the seven sites actually want, and that the allocation it
removes is not just msa.cpp's.

It also widens Option 3's blast radius: reimplementing
`parse_cigar_string` on top of the visitor changes the code path under all
seven, so all seven need verifying even though only two are being rewritten.

`print_uncompressed_cigar` (`cigar.cpp:185`) is the natural first follower —
it parses into a vector purely to loop over it once and print, which is the
exact pattern the visitor exists to remove:

```cpp
auto print_uncompressed_cigar(std::FILE * output_handle, View<char> const cigar_string) -> void {
  for_each_cigar_operation(cigar_string,
                           [output_handle](Operation const operation, long long const runlength) -> void {
                             auto const letter = convert_from_operation(operation);
                             for (auto i = 0LL; i < runlength; ++i) {
                               fprint(output_handle, static_cast<char>(letter));
                             }
                           });
}
```

That would need `convert_from_operation` exposed too, or the lambda can
keep calling it while it stays in `cigar.cpp` — `print_uncompressed_cigar`
is defined there, so no move is required for this one.

### Both msa.cpp walks become visitors

Walk A:

```cpp
auto find_max_insertions_per_position(View<struct msa_target_s> const targets,
                                      int const centroid_len) -> std::vector<int> {
  std::vector<int> max_insertions(static_cast<std::vector<int>::size_type>(centroid_len + 1));

  // the centroid (targets.front()) carries no cigar string, hence drop(1)
  for (auto const & a_msa_target : targets.drop(1)) {
    auto position = 0LL;
    for_each_cigar_operation(a_msa_target.cigar,
                             [&](Operation const operation, long long const runlength) -> void {
                               switch (operation) {
                               case Operation::match:
                               case Operation::insertion:
                                 position += runlength;
                                 break;

                               case Operation::deletion:
                                 assert(runlength <= std::numeric_limits<int>::max());
                                 max_insertions[static_cast<std::vector<int>::size_type>(position)] =
                                   std::max(static_cast<int>(runlength),
                                            max_insertions[static_cast<std::vector<int>::size_type>(position)]);
                                 break;
                               }
                             });
  }
  return max_insertions;
}
```

Note the signature does **not** change — walk A keeps taking `targets`.
That is the practical advantage of Option 3 over Option 2: no plumbing, no
parallel array, no index threading. The body changes; the interfaces do not.

Walk B, replacing `msa.cpp:368-423`:

```cpp
for_each_cigar_operation(target.cigar,
                         [&](Operation const operation, long long const runlength) -> void {
                           switch (operation) {
                           case Operation::deletion:
                             for (auto j = 0; j < runlength; ++j) { ... }
                             for (auto j = runlength; j < max_insertions[...]; ++j) { ... }
                             is_inserted = true;
                             break;
                           case Operation::match:
                             ...
                             break;
                           case Operation::insertion:
                             ...
                             break;
                           }
                         });
```

The `[&]` capture picks up `is_inserted`, `qpos`, `tpos`,
`position_in_alignment`, `target_seq`, `target_abundance`, `profile`,
`aln_v` and `max_insertions` by reference, which is what the loop body
already mutates in place. Neither walk uses `continue` or breaks out of the
outer `while`, so nothing is lost in the conversion — verified by reading
`msa.cpp:190-204` and `:380-422`.

### Cost

- Allocation in both walks: gone. Peak memory strictly lower than today,
  and much lower than Option 2.
- `cigar.hpp` gains a template and three includes, and `convert_to_operation`
  loses its deliberate TU-local scoping.
- Deeper indentation in both walks. Walk B's body is already long; wrapping
  it in a lambda pushes it further right and may argue for extracting the
  three cases into named helpers first.
- Inlining is a question, not a given. The visitor is a header template so
  the compiler *can* inline through it, but walk B's callback is large. If a
  callgrind run shows the walk got slower, that is the reason — check before
  concluding the option failed.

### Verification

Same as Option 2 (byte-identical `--msaout`/`--consout`/`--profile`), plus
coverage for all seven `parse_cigar_string` sites, because reimplementing it
on the visitor puts every one of them on new code:

- `--userfields aln` — `print_uncompressed_cigar` (`results.cpp:457`, case
  22) writes the *uncompressed* alignment. Note this is `aln`, **not**
  `caln`: `caln` (case 23) is the compressed cigar and does not go through
  this path.
- `--alnout` — `showalign.cpp` drives both of its sites.
- The chimera commands — `uchime_denovo`, `uchime2_denovo`,
  `uchime3_denovo`, `uchime_ref`, `chimeras_denovo` — cover
  `chimera.cpp:451/846/884`.
- `--msaout`/`--consout`/`--profile` — walk A and walk B.

Callgrind on a `--msaout` run and on a `--userout ... aln` run. The
allocation removal should show as a reduction at every one of the seven.


## Recommendation

If this is scheduled: **Option 3**, and the seven-call-site census above is
what decides it. Option 2 removes one redundant parse in one file, at the
cost of a parsed-cigar array threaded through two functions and indexed in
parallel with `targets` — the kind of coupling this codebase has been
removing, not adding. Option 3 removes the allocation at *six other* call
sites too, leaves `utils/cigar.hpp` with a single walking implementation,
and needs no signature change in `msa.cpp` at all: the bodies change, the
interfaces do not.

Option 2 is the right pick only if putting a template in `cigar.hpp` and
un-scoping `convert_to_operation` are judged too invasive for the payoff.

**But the honest recommendation is to schedule neither with any urgency.**
Nothing is wrong today. This is worth doing when someone is next in
`msa.cpp` for another reason, and worth doing *then* rather than as its own
campaign.


## Open questions

1. **Is `convert_to_operation` leaving its anonymous namespace acceptable?**
   Option 3 requires it. The existing comment states the TU-local scoping as
   an intent, so this reverses a documented decision.
2. **Should `cigar.hpp` depend on `fatal.hpp`?** Option 3 requires that too.
   It makes every cigar consumer transitively see the fatal-error API.
3. **Is the ill-formed-cigar `fatal()` the behaviour we want in walk B?**
   Both options give walk B the guard it currently lacks. Since these cigars
   are vsearch-generated, the guard should never fire — but if it ever did,
   both options turn a silently-tolerated malformed cigar into a hard exit.
   That is the correct direction, and it is a behaviour change on a path
   that is currently unreachable.
