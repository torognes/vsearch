# Plan: raw loops → range-for and standard algorithms (2026-08-07)

> **DONE (2026-08-08), branch `tmp_20260806093624`, 13 commits from
> `511ed257`.** All ten phases landed. The plan text below is unchanged
> from what was approved; what actually happened, and the four places
> where the plan was wrong, are recorded in *Outcome* at the end.

## Goal

The last three campaigns moved vsearch's *storage* to the standard
library (`manual_memory`, `replace_C_functions`) and its *signatures* to
the vocabulary types (`view_span_propagation`, `propagate_containers`,
landed earlier today). The **bodies** did not always follow.

A function that receives a `View<struct hit>` and then walks it with

```c++
for (int64_t t = 0; t < toreport; t++)
  {
    auto const & hit = hits[static_cast<std::size_t>(t)];
```

has paid the whole cost of the new type and collected none of the
benefit: the bound is still a second variable that can disagree with the
view, the cast is back, and the reader still has to check that `t` is not
used for anything else before believing the loop is a simple traversal.

This campaign closes that gap. It changes no signature that the container
campaign already settled; it rewrites what happens *inside*.


## The rule this campaign applies

> **A loop should name what it does.** Prefer, in order:
>
> 1. a standard algorithm, when one says it in a word
>    (`std::count_if`, `std::max_element`, `std::fill_n`, `std::equal`,
>    `std::find_if`, `std::transform`);
> 2. a range-for over the *exact* range, sliced with `first()` / `drop()`
>    / `subspan()` rather than guarded by an `if` inside the body;
> 3. an index loop — but only when **the index is itself data**.

The third clause is what keeps this from becoming a mechanical sweep.
There are three legitimate reasons to keep an index, and all three are
already documented at their sites:

- the index is **printed** (`--uc` cluster numbers, `derep.cpp:362`);
- the index **selects behaviour** (`results.cpp:1081` — the SAM flag
  marks every hit after the first as secondary, and says so in a comment
  added when that loop deliberately kept its index);
- the index is a **coordinate in a flat 2-D array** (`chimera.cpp`'s
  `cand * query_len + qpos`, `msa.cpp`'s `profsize * i + nucleotide`) —
  though see Phase 7, where the right fix is to stop having
  a flat array with hand-computed coordinates.

Everything else in the census below is a traversal wearing an index.


## Census (2026-08-07, `src/`, excluding `vendored/`)

| | count |
|---|---:|
| index-style `for (a; b; c)` loops | **309** |
| range-for loops | 91 |
| subscripts of the form `[static_cast<...>(index)]` | **379** |

Where the index loops are:

| file | index loops | note |
|---|---:|---|
| `core/chimera.cpp` | 53 | flat 2-D arrays + `count` members — **own campaign** (Decision 3) |
| `core/align_simd.cpp` | 39 | **out of scope**, as in the C-style plan |
| `commands/sintax.cpp` | 17 | Phase 4 |
| `core/cluster.cpp` | 15 | Phase 8 |
| `core/msa.cpp` | 12 | Phase 7 |
| `core/results.cpp` | 11 | Phases 1, 3 |
| `commands/fastq_eestats2.cpp` | 11 | not surveyed; own plan |
| `core/udb.cpp`, `core/linmemalign.cpp` | 10 each | index-is-data / own plan |
| `core/mask.cpp` | 8 | Phase 10 (hot) |
| `core/mergepairs.cpp` | 8 | not surveyed; own plan |
| `core/derep.cpp` | 7 | index-is-data (see Leave alone) |
| `commands/search_exact.cpp` | 7 | Phases 1, 2 |
| `core/search.cpp`, `core/searchcore.cpp` | 6 each | Phases 2, 5, 6 |
| `commands/usearch_global.cpp` | 6 | Phases 1, 2 |
| `core/dbindex.cpp` | 5 | index-is-data (see Leave alone) |

`core/chimera.cpp` alone holds 158 of the 379 cast-laden subscripts, and
`core/mergepairs.cpp` another 39. Those two are the tail of this
distribution and are treated accordingly: chimera is surveyed here and
run as its own campaign, mergepairs is out of scope entirely.


## What "already done" looks like

Five bodies landed in the correct shape, in the same commits that changed
their signatures. They are the templates the phases below copy, and each
should be read before the phase that cites it:

| site | shape |
|---|---|
| `utils/string_normalize.cpp:69`, `utils/reverse_complement.cpp:69` | `std::transform` over `View` → `Span`, with a `.rbegin()/.rend()` for the reverse case. No loop at all. |
| `core/mask.cpp:251` `hardmask()` | `std::transform` in place over a `Span<char>` |
| `core/msa.cpp:158` `find_max_insertions_per_position` (`ec5cc3cf`) | `for (auto const & target : targets.drop(1))` — the "first element is special" case handled by slicing, not by an `if` in the body |
| `core/unique.cpp:274` `Uniquer::count_shared` (`93e28644`) | range-for over the returned `View<unsigned int>` |
| `core/results.cpp:709` `results_show_alnout` (`32e786ec`) | range-for over `View<struct hit>` |

Note what `find_max_insertions_per_position` does *not* do: it does not
range-for over the whole view and `continue` on the centroid. It slices.
That distinction is the substance of Phase 1.


---

# Phases

Ordered cheapest-and-safest first, one family per commit, each verified
byte-identical before the next — the ordering the C-style migration used
and that worked.

## Phase 1 — the top-hits prefix (7 sites)

The single most repeated loop in the tree. Every writer that reports a
query's hits carries the same three lines:

```c++
auto const top_hit_id = hits.front().id;
for (auto const & hit : hits)
  {
    if ((parameters.opt_top_hits_only != 0) and (hit.id < top_hit_id)) { break; }
```

at `core/results.cpp:592` (lcaout), `:741` and `:758` (alnout, twice),
`:1085` (samout), `commands/usearch_global.cpp:217`,
`commands/search_exact.cpp:274`, `commands/allpairs_global.cpp:189`.

`--top_hits_only` does not filter the hits; it **truncates** them at the
first one that is not tied with the best. That is a prefix, so it is a
`View`, and `std::find_if` computes it:

```c++
/* core/searchcore.hpp, next to struct hit */

/* The hits to report for this query: all of them, or -- under
   --top_hits_only -- only the leading run tied with the best identity.
   Seven writers used to spell this as a break in the middle of their
   loop. find_if stops at the first hit that fails the test, exactly as
   the break did, so the result is unchanged even for an unsorted list. */
inline auto top_hits(View<struct hit> const hits, bool const top_hits_only)
  -> View<struct hit> {
  if ((not top_hits_only) or hits.empty()) { return hits; }
  auto const best_id = hits.front().id;
  auto const * const first_worse =
    std::find_if(hits.begin(), hits.end(),
                 [best_id](struct hit const & candidate) -> bool {
                   return candidate.id < best_id;
                 });
  return hits.first(static_cast<std::size_t>(std::distance(hits.begin(), first_worse)));
}
```

Each site then loses its `top_hit_id` local, its `if`, and its `break`.

The three commands gain a second thing. Each of them already spells the
`--maxhits` clamp three or four times:

```c++
/* commands/usearch_global.cpp:170,175,185,194 */
auto const toreport = std::min(state.parameters.opt_maxhits, static_cast<int64_t>(hits.size()));
...  make_view(hits).first(static_cast<std::size_t>(toreport)),   /* x3 */
...  for (int64_t t = 0; t < toreport; t++)                        /* and then again by index */
```

One named local — `auto const reported = top_hits(make_view(hits).first(...), ...)`
— replaces the three copies and the index loop at once. The same shape is
at `commands/search_exact.cpp:236-270` and
`commands/allpairs_global.cpp:158-185`.

Note `results_show_lcaout` also carries a `tophitcount` counter that
exists only to re-derive this same prefix for its second pass
(`core/results.cpp:647`, `hits.first(tophitcount)`); with `top_hits()` at
the top of the function both passes iterate the same named view and the
counter becomes `top.size()`.

**Risk: low.** Pure restructuring on the output path, no arithmetic
changes. Fully covered by the differential verification.

**Deliverable:** 7 breaks removed, 3 duplicated clamp expressions
removed, 1 new inline helper.

---

## Phase 2 — the two-strand loop (9 sites)

```c++
for (int s = 0; s < number_of_strands(parameters.opt_strand); s++)
  {
    struct searchinfo_s * strand_si = (s != 0) ? si_minus : si_plus;
```

at `core/search.cpp:278,325,423,474`, `core/cluster.cpp:869,958,1867`,
`commands/usearch_global.cpp:376`, `commands/search_exact.cpp:430,489`,
`commands/sintax.cpp:440,560`. Twelve occurrences of a loop over a
two-element sequence that is not stored as a sequence, so the body has to
reconstruct the element from the counter.

The vocabulary types already handle this; the array is a local, so
nothing dangles:

```c++
std::array<struct searchinfo_s *, 2> const strands {{si_plus, si_minus}};
for (auto * const strand_si :
     make_view(strands).first(static_cast<std::size_t>(number_of_strands(parameters.opt_strand))))
  {
```

**Two sites must keep the counter and are excluded**:
`commands/search_exact.cpp:489` and `commands/sintax.cpp:440` assign
`si->strand = s` / use `s` as a strand index into `all_seqno` and
`boot_count`, so there the index *is* data. They keep the index loop and
get a one-line comment saying why, in the manner of `results.cpp:1081`.

A `strands_of(si_plus, si_minus, opt_strand)` helper returning the view
was considered and is **not** proposed: the array has to be a named local
at the call site anyway for its lifetime, so the helper would save
nothing and hide the lifetime.

**Risk: low.** No arithmetic; `number_of_strands()` is called once per
loop instead of once per iteration, which is if anything cheaper.

---

## Phase 3 — `core/results.cpp`: the taxonomy levels

Three separate defects behind one function, `results_show_lcaout`
(`:559`), plus the `tax_split` signature it shares with sintax.

**3a — `tax_split`'s two out-arrays.** `core/tax.hpp:64`:

```c++
auto tax_split(int seqno, int * level_start, int * level_len, struct Database const & db) -> void;
```

Three call sites (`core/results.cpp:602,652`, `commands/sintax.cpp:180`),
all of which pass `.data()` of two `std::array<int, tax_levels>` locals
that were declared next to each other for the purpose. This is exactly
the shape Phase 4 of the container campaign fixed for
`get_hex_seq_digest_*` — a `T *` out-parameter that is always the whole
of one fixed-size array — and it was missed because the grep for it
looked for `Span`, not for `int *`.

The maintainer's own brief is at `core/results.cpp:600`:

```c++
std::array<int, tax_levels> new_level_start {{}};  // refactoring: std::array<struct a_level{.start, .length}, tax_levels>
```

Carry it out. Two parallel arrays with the same length and the same
index become one array of a two-field struct:

```c++
struct TaxLevel { int start = 0; int length = 0; };  /* a window into the header */
auto tax_split(int seqno, std::array<TaxLevel, tax_levels> & levels,
               struct Database const & db) -> void;
```

The wrong-size buffer stops compiling instead of running off the end, and
the two arrays can no longer be filled to different depths.

> **Considered and rejected: `std::array<View<char>, tax_levels>`.**
> Every consumer immediately converts the (start, length) pair into
> `db.header_view(seqno).subspan(start, length)`, so the views look like
> the honest type — and `commands/sintax.cpp:162` has *already* made that
> move locally, storing `std::array<std::array<View<char>, tax_levels>, bootstrap_count>`.
> But `results_show_lcaout` keeps its candidate levels across iterations
> and re-derives them from a *different* `seqno` each time
> (`:614-620`), so a view stored in round 1 would be compared against a
> header fetched in round 3. The offsets are the safe currency inside
> `tax_split`; the views are the right thing for the caller to build
> immediately afterwards, which is what sintax does. Keep the split.

**3b — the level-match loops → `std::equal`.** `results.cpp:619` and
`:657` are the same eight lines twice:

```c++
auto match = true;
for (std::size_t j = 0; j <= k; ++j)
  {
    auto const cand_level = db.header_view(...).subspan(cand_level_start[k][j], cand_level_len[k][j]);
    auto const query_level = db.header_view(...).subspan(new_level_start[j], new_level_len[j]);
    if (cand_level != query_level) { match = false; break; }
  }
```

"Do the first k+1 levels agree" is `std::equal` with a predicate, over
the `TaxLevel` arrays 3a introduces:

```c++
auto const match = std::equal(cand_levels[k].begin(), std::next(cand_levels[k].begin(), k + 1),
                              new_levels.begin(), same_level_name);
```

where `same_level_name` is a lambda capturing the two headers. The flag,
the `break` and eight casts go with it, twice.

**3c — the `userout` field separator.** `results.cpp:362`:

```c++
for (std::size_t c = 0; c < userfields_requested.size(); ++c)
  {
    if (c != 0) { fprint(output_handle, '\t'); }
    auto const field = userfields_requested[c];
```

The index exists only to suppress the leading tab. The body is a
90-case switch, so the clean fix is to lift it into
`print_userfield(output_handle, field, ...)` and write

```c++
if (not userfields_requested.empty()) { print_userfield(output_handle, userfields_requested.front(), ...); }
for (auto const field : make_view(userfields_requested).drop(1)) {
  fprint(output_handle, '\t');
  print_userfield(output_handle, field, ...);
}
```

**Confirmed in scope (Decision 2), against the recommendation this plan
first made.** It is a large mechanical move of a switch into a new
function, so it gets its own commit, moved verbatim — no field's
formatting is touched in the same commit that moves it, or the
differential run cannot tell a transcription slip from an intended
change. The `bool first` alternative was rejected: it removes the index
without making the loop say anything more than it already does.

Two things to watch while moving it. The switch reads `hit` (possibly
null), both query views, the reverse-complement view, the database and
the parameters, so `print_userfield` takes the same argument list the
enclosing function has — resist trimming it "while we are here". And
`userfields_requested` is `parameters.opt_userfields`, whose element type
must be checked before `make_view()`: `drop(1)` on an empty container is
already safe (it clamps), but the `front()` call above it is not, hence
the `empty()` guard.

**Risk: low-medium.** 3a touches a header and three call sites in two
commands; 3b changes the shape of a comparison that decides taxonomy
output; 3c is the widest-blast-radius commit in Phases 1-9 by line count,
though not by subtlety. `--lcaout`, `--userout` (all 19 fields) and
`--sintax` all need the differential run.

---

## Phase 4 — `sintax_analyse`: the last raw pair, and two algorithms

`commands/sintax.cpp:149`:

```c++
static auto sintax_analyse(struct sintax_state_s & state,
                    char const * query_head,
                    int const strand,
                    int const * all_seqno,
                    int const count) -> void
```

`(int const *, int)` is row 1 of the container campaign's table, still
here because the census counted `struct hit *` and `char *` parameters
and this one is `int *`. The call site (`:524-528`) already has the
container:

```c++
sintax_analyse(state, query_head, best_strand,
               all_seqno[static_cast<std::size_t>(best_strand)].data(),
               boot_count[static_cast<std::size_t>(best_strand)]);
```

→ `make_view(all_seqno[best_strand]).first(boot_count[best_strand])`, and
the body's `for (auto i = 0; i < count; i++)` (`:173`) becomes a
range-for over the view.

Then two algorithms inside the same function:

**`:225` is `std::max_element`.**

```c++
for (auto i = 0; i < count ; i++) {
  if (cand_matchcount[cand_i] > level_matchcount[level]) {
    level_best[level] = i;
    level_matchcount[level] = cand_matchcount[cand_i];
  }
}
```

Strictly-greater keeps the first of a tie, which is what `std::max_element`
returns, so the two agree element for element. Careful: `level_best`
starts at `-1` and must stay `-1` when every count is zero — so the
rewrite is `max_element` plus the explicit `> 0` test, not `max_element`
alone. That subtlety is why this one deserves its own commit and its own
line in the message.

**`:206-219` is `std::find_if`.** The inner `for (j = 0; j <= i; j++)`
with a `break` at the first match is "find the first still-included
earlier candidate with the same name at this level"; `find_if` over
`make_view(cand_level_name).first(i + 1)` says it, and the `break`
disappears.

**`:182` and `:255`/`:280`** iterate `for (k = 0; k < tax_levels; k++)`
and immediately compute `static_cast<std::size_t>(k)` — three loops whose
only use of the counter is to index one `std::array`. Plain range-fors.

**Risk: medium.** `--sintax` output is a consensus vote; a
first-of-a-tie change is invisible on most inputs and glaring on the
right one. Verify against a corpus with deliberate ties, not only the
test suite.

---

## Phase 5 — `core/search.cpp`: filling the result span

The library batch path, re-typed this morning in `8e6fbabd`, still fills
its `Span<struct search_result_s>` with a hand-rolled counter — twice,
in near-identical copies at `:302` (`search_session_single`) and `:452`
(the batch worker):

```c++
int count = 0;
for (auto const & h : hits)
  {
    if (static_cast<std::size_t>(count) >= results.size()) { break; }
    auto & r = results[static_cast<std::size_t>(count)];
    r.target = h.target;
    ...  /* ten fields */
    ++count;
  }
```

The bound is known before the loop, so it is a `std::min` and a
`std::transform`:

```c++
auto const reported = std::min(hits.size(), results.size());
std::transform(hits.begin(), std::next(hits.begin(), static_cast<std::ptrdiff_t>(reported)),
               results.begin(), to_search_result);
return static_cast<int>(reported);
```

with `to_search_result` a file-local lambda (it needs `si->db` and the
query length, both fixed for the query) — which also **deletes the
second copy**, since the two loops differ only in where the span comes
from.

**Risk: low.** The library batch APIs are the least-exercised code in
the tree, so this needs the release-build `api_examples` run (a debug
`libvsearch.a` segfaults the suite on the `_GLIBCXX_DEBUG` ABI
mismatch), not just the CLI suite.

---

## Phase 6 — `core/searchcore.cpp`: `align_delayed`

`:773` and `:798`, both:

```c++
for (auto x = static_cast<std::size_t>(searchinfo->finalized); x < hits.size(); ++x)
  {
    auto const & hit = hits[x];
```

The starting offset is exactly what `Span::drop()` is for:

```c++
for (auto const & hit : hits.drop(static_cast<std::size_t>(searchinfo->finalized)))
```

The second loop needs a mutable `struct hit *`; `make_hits_span()`
returns `Span<struct hit>`, whose `begin()` is `Type *`, so `auto & hit`
gives it directly and the `&hits[x]` goes away.

**Risk: low, but this is inside the search hot path.** `drop()` clamps
with `std::min`, so it adds one comparison per *loop*, not per iteration.
Callgrind the instruction count anyway — see Verification; this is
precisely the size of change that the 2026-08-04 alignment trap made
look like a 4 % regression.

---

## Phase 7 — `core/msa.cpp`: the profile columns

**7a — two fills.** `compute_and_print_consensus` (`:415`) opens by
censoring both ends of the alignment:

```c++
for (auto i = 0; i < left_censored; ++i)
  { aln_v[static_cast<std::vector<char>::size_type>(i)] = '+'; }
for (auto i = alignment_length - right_censored; i < alignment_length; ++i)
  { aln_v[static_cast<std::vector<char>::size_type>(i)] = '+'; }
```

`std::fill_n(aln_v.begin(), left_censored, '+')` and a `std::fill` over
the matching tail. Four casts gone, and the second loop's bound
arithmetic becomes visible as a `last()` slice.

**7b — the profile as columns.** `profile` is a flat
`std::vector<prof_type>` of `profsize`-wide columns, and every access
spells the stride by hand:

```c++
profile[(static_cast<std::vector<prof_type>::size_type>(profsize) *
         static_cast<std::vector<prof_type>::size_type>(i)) + nucleotide]
```

— three times in ten lines at `:444-464`, and again in `update_profile`.
The maintainer's brief at `core/msa.cpp:573` proposes
`std::vector<std::array<prof_type, profsize>>` and marks it C++20; it
does not need to wait, because a column accessor gives the same call-site
shape in C++11:

```c++
/* the profsize counters of alignment position `index` */
inline auto column(std::vector<prof_type> & profile, int const index) -> Span<prof_type> {
  return make_span(profile).subspan(static_cast<std::size_t>(index) * profsize, profsize);
}
```

The "find most common of A, C, G, T" loop at `:444` is then
`std::max_element(col.begin(), col.begin() + 4)` — with the same
first-of-a-tie caveat as Phase 4, and the same `best_count == 0` fallback
to N kept explicitly. The three monster subscripts become `col[4]` and
`col[5]`.

**This is the one brief under Decision 1 that is rewritten rather than
deleted.** The accessor reaches the C++20 form's call-site shape, but it
does not *become* it: `std::vector<std::array<prof_type, profsize>>` also
makes the stride a type-level fact and lets the compiler see the column
width, which a `subspan` cannot. So the brief stays, narrowed to what is
still outstanding — the storage type, not the call sites, which this
phase settles.

**Risk: medium.** `--msaout`/`--consout`/`--profile` are exercised by the
cluster commands and were verified byte-identical this morning, so a
baseline is at hand. The consensus tie-break is the thing to watch.

---

## Phase 8 — `core/cluster.cpp`

**8a — `clusterinfo_v`, three loops.** `:1288`, `:1357`, `:1526` all run
`for (int i = 0; i < seqcount; i++)` over a vector constructed as
`std::vector<clusterinfo_t> clusterinfo_v(static_cast<std::size_t>(seqcount))`
(`:1246`) — so `seqcount` *is* `clusterinfo_v.size()`, and the loops
between them perform ten `static_cast<std::size_t>(i)` subscripts. Plain
range-fors.

One check before each: `Progress::update()` is called with the ordinal in
some of these loops. Where it is, keep an explicit counter or hoist the
progress update — do not silently drop the progress bar's resolution.

**8b — `evaluate_extra_hits`'s raw pair.** `:613`:

```c++
static auto evaluate_extra_hits(struct searchinfo_s * si,
                                struct searchinfo_s * si_plus,
                                const int * extra_list,
                                int const extra_count, ...)
```

Two call sites (`:962`, `:1871`), both passing `extra_list.data()` of a
`std::vector<int>` and a separate fill level → `View<int>`, and the
`for (int j = 0; j < extra_count; j++)` at `:642` becomes a range-for.
Same shape as Phase 4, found by the same widened grep.

**Risk: low-medium.** `--cluster_fast`/`--cluster_size`/`--cluster_smallmem`,
all three, differentially.

---

## Phase 9 — the maintainer's own algorithm briefs

Four `// refactoring:` comments in the tree name a standard algorithm
outright. They are cheap, they are already sanctioned, and they are
scattered — so they make one small commit each, or one commit together if
the maintainer prefers.

| site | brief | what it becomes |
|---|---|---|
| `core/filter.cpp:217` | `/* filter by n's */  // refactoring: std::count_if();` | `std::count_if` over `input_handle->sequence_view().subspan(start, length)` — which also removes the raw `get_sequence() + start` pointer arithmetic, since the reader has had a `sequence_view()` since the fastx-reader migration |
| `commands/fastq_convert.cpp:109` | a five-line brief spelling out subspan → subtract → clamp → add → clamp | one `std::transform` with the clamp chain in the lambda |
| `commands/sff_convert.cpp:605` | `// refactoring: soft_mask_read(transform(begin(), left_mask_end); transform(right_mask_start, end()))` | two `std::transform`s over subspans |
| `core/chimera.cpp:731` | `// refactoring: reset or initialize?` | answer it in the comment; the `std::fill` beneath it is already correct |

`core/filter.cpp` is on the `--fastq_filter` / `--fastx_filter` hot path
and needs the same callgrind treatment as Phase 6.

---

## Not a phase — `core/chimera.cpp`: the flat 2-D arrays

**Decision 3: this is its own campaign, not part of this branch.** The
survey is kept here in full so the next plan starts from it rather than
re-deriving it, and so the census above is not read as a sweep that
missed the biggest file in it.

53 index loops, 158 cast-laden subscripts — half the tree's total in one
file. The cause is a single design decision repeated: `chimera_info_s`
(`:129`) holds `match`, `insert`, `smooth`, `maxsmooth`, `scan_p`,
`scan_q` as flat `std::vector<int>`/`<double>` of
`cand_count * query_len` elements, plus ten parallel
`std::array<..., maxcandidates>` alignment-result arrays behind one
`cand_count`, plus `best_parents`/`best_start`/`best_len` behind
`parents_found`.

Every access recomputes the coordinate:

```c++
size_t const z = (static_cast<size_t>(i) * static_cast<size_t>(ci->query_len)) + static_cast<size_t>(qpos);
```

The fix is the same accessor idea as Phase 7b, one level up — a
`match_row(cand)` / `smooth_row(cand)` returning `Span<int>` — after
which `:761-778` ("count the wins") is a range-for, `:782-789` ("select
best parent") is `std::max_element` (again with the "stay -1 if all zero"
guard), and `:711-724` ("wipe out matches covered by the previous
parent") is a `std::fill` over a subspan of each row.

**Why it is separate.** Half the tree's problem in its least-covered
engine is not something to append to a nine-phase branch: the arrays are
indexed from four different functions, and `--uchime_denovo` output is a
scored classification, so an off-by-one changes a verdict rather than
crashing. Running Phases 1-9 first also tests the row-accessor idiom on
`msa.cpp`'s profile (Phase 7b), which is the same shape at a tenth of the
scale — so the chimera campaign can start from a proven idiom instead of
proposing one.

**Two things that campaign should settle before it starts**, both
visible from here: whether the ten parallel
`std::array<..., maxcandidates>` alignment-result arrays want to become
one array of a record (as `query_record_s` did for the library batch API
this morning), and whether `search16()`'s eight parallel array
parameters — which those ten arrays exist to feed — can be re-typed
without opening `align_simd.cpp`, which is excluded by standing
decision.

---

## Phase 10 — the hot paths: measure, then decide

Two conversions are correct, obvious, and **must not be assumed
free**. Both are in the innermost loop of a command users run at scale.

**10a — `core/mask.cpp` `dust_core` (`:144`).** `dust()` takes a
`Span<char>` (`:196`) and immediately unwraps it into `(char *, int)` for
`dust_core`, after which the file is raw:

```c++
for (auto i = 0; i < len; i++) { seq[i] = to_upper(seq[i]); }        /* :158 → std::transform */
for (auto j = a + i; j <= b + i; j++) { seq[j] = 'N'; }              /* :174 → std::fill over a subspan */
for (auto j = a + i; j <= b + i; j++) { seq[j] = local_seq[j] | 32U; } /* :181 → std::transform */
```

The outer window loop at `:165` mutates `i` inside the body
(`i += half_dust_window - b`) and **must stay an index loop** — say so in
a comment when the others go.

`wo()` (`:83`) additionally takes `(int len, char const * s, int * beg, int * end)`
— a length/pointer pair and two out-parameters — where
`View<char>` and a returned `struct DustRegion { int score; int begin; int end; }`
would say it. That is a signature change on the hottest function in the
file and belongs in this phase only if the measurement in 10a comes back
clean.

**Why the caution:** `wo()` measured at **74.65 % of an entire benchmark
run** during the 2026-08-04 campaign, and in the same campaign a
byte-identical `wo()` moved 4 % purely because an unrelated function grew
54 bytes and shifted it off a 64-byte boundary.

**10b — `core/unique.cpp` `count_hash` (`:169`).** The k-mer primer and
main loops (`:216`, `:228`) walk `s` from `seq.data()` to two computed
sentinels; they are `seq.first(wordlength - 1)` and the matching
`drop()`, and would read far better as two range-fors. `Uniquer::count()`
is the innermost k-mer loop of *every* search command, and the file
already carries a comment recording that widening a single loop bound
from 32 to 64 bits cost **+8 %** on `search_topscores`
(`core/searchcore.cpp:324`).

**Both 10a and 10b run measure-first, revert-if-negative (Decision 4).**
Landing them is not the goal; knowing whether they cost anything is. If
either costs more than noise, the commit becomes a comment at the site
recording the measurement, which is worth as much as the change.

---

# Leave alone (so a later reader does not read these as a missed sweep)

- **`core/results.cpp:1081`** — `results_show_samout` keeps its index and
  already explains why: the SAM flag marks every hit after the first as a
  secondary alignment (`0x100`), so the position is output.
- **`core/derep.cpp:248,284,320,358,383`, `commands/derep_prefix.cpp:367,389,421,458`**
  — the `--uc` and `--tabbedout` writers print the cluster ordinal
  (`fprint_integer(fp_uc, i)`) and feed `progress.update(i)`. Index is
  data.
- **`core/dbindex.cpp:157,185,200,221`** — loops over sequence numbers and
  hash slots, where the counter *is* the key. The one loop in that file
  that was a traversal (`add_sequence`'s k-mer walk) already became a
  range-for in `93e28644`.
- **`core/align_simd.cpp`, 39 loops** — same standing exclusion as the
  C-style migration (`DONE_20260804:291`): the SIMD kernels and their
  `#if 0` debug dumps. Not a target.
- **`commands/search_exact.cpp:489`, `commands/sintax.cpp:440`** — strand
  loops whose counter is a strand index used as data (Phase 2).
- **`core/mergepairs.cpp` (8 loops, 39 casts) and
  `commands/fastq_eestats2.cpp` (11 loops)** — real candidates, but
  neither received a container parameter in the recent campaigns, so
  neither is in this campaign's remit. They deserve a plan of their own.

---

# Signature residue found while surveying

Not part of this campaign's remit — this campaign changes bodies — but
found while reading and worth recording before it is forgotten. Each is
row 1 or row 2 of the container campaign's table, missed because that
census grepped for `struct hit *` and `char *`:

| site | shape | note |
|---|---|---|
| `core/tax.hpp:64` `tax_split` | two `int *` out-arrays | Phase 3a fixes it, because Phase 3b cannot be written without it |
| `commands/sintax.cpp:152` `sintax_analyse` | `int const *, int count` | Phase 4 fixes it, same reason |
| `core/cluster.cpp:615` `evaluate_extra_hits` | `int const *, int count` | Phase 8b |
| `core/derep.cpp:240,276,310,374,...` the report writers | `std::vector<struct bucket> const &, uint64_t clusters` | **row 2, not fixed by any phase here.** The hashtable is sorted so the filled buckets lead, and `clusters` is the fill level — the same "reused high-water-mark buffer" the campaign described. `View<struct bucket>` is the type. Left out because the loops themselves must keep their index (see Leave alone), so there is no body work to pair it with; it is a one-commit follow-up to the *container* campaign, not to this one. |
| `core/mask.cpp:83` `wo()` | `(int, char const *, int *, int *)` | Phase 10a, gated on the measurement |

---

# Verification

Per phase, in this order. No phase lands until the phase before it is
verified.

1. **Byte-identical differential.** `git worktree add --detach` a
   baseline (vsearch-tests already holds `dev`, so a second checkout is
   mandatory), run both binaries at **`--threads 1`** — otherwise the
   diff is reordering noise — strip the command-line first line, and
   `diff -r`. For `--samheader`, the command line also appears on the
   `@PG` line, which has to be filtered too.
   Per-phase commands: Phase 1 → `--usearch_global` with
   `--alnout --blast6out --uc --userout --lcaout --samout --top_hits_only --maxhits 3`,
   `--search_exact`, `--allpairs_global`; Phase 3 → add `--lcaout` against
   a taxonomy-annotated database; Phase 4 → `--sintax`; Phase 7 →
   `--cluster_size --msaout --consout --profile`; Phase 8 → all three
   cluster commands.

2. **Instruction count before wall-clock**, for Phases 5, 6, 9
   (`filter.cpp`) and 10:

   ```sh
   valgrind --tool=callgrind --callgrind-out-file=cg.out ./bin/vsearch <args>
   callgrind_annotate cg.out | grep -E '^ *[0-9,]+ ' | head -8
   ```

   Equal instruction counts across the change mean any wall-clock delta is
   code *placement*, not work — the trap that produced two wrong
   conclusions on 2026-08-04. Read the function list every time: it also
   tells you whether the benchmark measures what you think.
   Use `--qmask none --dbmask none` on anything that is not
   deliberately measuring DUST, or the masking dilutes the measurement
   about fourfold. `perf` does not work here (`perf_event_paranoid` is 4).

3. **`hyperfine`** on release builds, only after (2) says the work
   changed. Optionally build both sides with `-falign-functions=64` to
   remove the placement noise (it is not free — +2.5 % `.text`, 0-1.3 %
   *slower* — so it is an A/B aid, never a release flag).

4. **`api_examples` on a release build** for Phase 5:
   `make -C src libvsearch.a` is a separate explicit target the default
   `make` does not build, and a debug library segfaults the suite on the
   `_GLIBCXX_DEBUG` ABI mismatch.

5. **clang-tidy before/after on each touched file** — it reports current
   state, not a diff, so run it on the baseline commit too. Expect
   `pointer-arithmetic`, `identifier-length` and the
   `static_cast`-in-subscript families to fall, and
   `unchecked-container-access` to *rise* wherever a raw pointer subscript
   becomes a `View`/`Span`/`std::array` one — that check flags the bounds
   assertion it does not know about, as it did in `32e786ec` and
   `8e6fbabd`. Say so in the commit rather than treating it as a
   regression.

6. **cppcheck** on each touched `.cpp`, the four cross-compiles
   (mingw / POWER / RISC-V) and the full
   `vsearch-tests/run_all_tests.sh`. Two cases in that suite
   (`getseqs --labels` log, `mergepairs --threads 1024`) fail spuriously
   under concurrent build load — re-run the single script in isolation
   before believing a FAIL.

---

# Decisions taken (maintainer, 2026-08-07)

All four settled before the branch opened. Two of them changed the shape
of this plan; both are recorded with what they overruled, so a later
reader does not re-open a closed question.

1. **The carried-out `// refactoring:` comments are deleted, and their
   rationale is kept.** Each of the six briefs is implemented, the brief
   itself is removed, and whatever it *explained* survives as a comment
   on the new code — the `ec5cc3cf` treatment, which carried out
   `find_max_insertions_per_position`'s brief, deleted it, and moved its
   residue into a doc comment on the new parameter. Affected:
   `core/results.cpp:600` (3a), `core/msa.cpp:573` (7b),
   `core/filter.cpp:217`, `commands/fastq_convert.cpp:109`,
   `commands/sff_convert.cpp:605`, `core/chimera.cpp:731` (9).

   Note that two of these are marked for a future standard —
   `msa.cpp:573` says C++20, `os_byteswap.hpp` and friends say C++23.
   Only `msa.cpp:573` is touched, and only because a `column()` accessor
   reaches the same call-site shape in C++11; the brief is **not**
   deleted on the grounds that the C++20 form is still the better one
   when it becomes available. It is rewritten to say so.

2. **Phase 3c runs: the `userout` switch is lifted into
   `print_userfield()`.** This overrules the plan's own recommendation,
   which was to leave the loop indexed on the grounds that the move is
   large and buys little. It gets its own commit, and the switch moves
   verbatim — no field's formatting changes in the commit that moves it,
   or `--userout`'s 19-field differential cannot distinguish a
   transcription slip from an intended change.

3. **`core/chimera.cpp` is its own campaign**, not Phase 10 of this one.
   The survey stays in this document (see *Not a phase*, above) so the
   next plan starts from it. This branch is therefore **Phases 1-9 plus
   the measure-first Phase 10**, and nothing in it touches
   `chimera.cpp` except the one-line comment answer at `:731` under
   Decision 1.

4. **Phase 10 is measure-first, revert-if-negative.** `dust_core` and
   `count_hash` are converted, measured with callgrind *before*
   `hyperfine`, and kept only if the instruction count says the work did
   not change. Where it did change for the worse, the commit becomes a
   comment at the site recording the measurement and the raw loop stays
   — which is the more useful artefact either way, because it stops the
   next sweep from re-litigating it. `wo()`'s signature change (10a's
   second half) is gated on the first half measuring clean.

---

# Out of scope (worth their own plans)

- **`core/chimera.cpp`** — Decision 3. Surveyed in full above (*Not a
  phase*), which is where the next plan should start; the biggest single
  item in the census, and the only one excluded by a decision rather than
  by never having been re-typed.
- **`core/mergepairs.cpp`** (39 cast-laden subscripts) and
  **`commands/fastq_eestats2.cpp`** (11 index loops) — neither was
  re-typed by the container campaign; both are real candidates.
- **`core/align_simd.cpp`** — standing exclusion, as above.
- **`core/linmemalign.cpp`, `core/udb.cpp`** (10 index loops each) — the
  first is a DP matrix (coordinates are data), the second is a binary
  format walk.
- **`search16()`'s eight parallel arrays** (`core/align_simd.hpp`,
  called from `core/searchcore.cpp:784` and `core/chimera.cpp`) — the
  parallel-array shape Phase 5a of the container campaign fixed for the
  library API, still present on the SIMD boundary. It is a signature
  change into `align_simd.cpp`, which is out of scope by standing
  decision; noted so the exception is on record rather than looking like
  an oversight.
- **`Database`'s missing record range.** `results.cpp:978`,
  `usearch_global.cpp:751,781`, `search_exact.cpp:815,845` and
  `mask.cpp:269` all run `for (i = 0; i < db.getsequencecount(); ++i)`.
  Some of those genuinely need the index (they print it, or index a
  parallel `dbmatched` vector), but not all. A `db.records()` range would
  be a bigger and more interesting change than anything in this plan, and
  should not be smuggled in as a loop cleanup.

---

# Outcome (2026-08-08)

Branch `tmp_20260806093624`, 13 commits from `511ed257`, 19 files,
+810/-623. One commit per phase, each verified before the next, plus a
final commit for the clang-tidy findings the campaign itself introduced.

| phase | commit | |
|---|---|---|
| 1 | `16bb2cf1` | `top_hits()` + one named `to_report` view; 7 breaks, 3 duplicated clamps gone |
| 2 | `1fcc808f` | 9 strand loops as `std::array` slices, 3 keep the counter |
| 3a/3b | `4c01081c` | `TaxLevel`, `std::equal`; `taxonomic_fields.h` gains `#pragma once` |
| 3c | `372f223a` | the `--userout` switch moved verbatim into `print_userfield()` |
| 4 | `9357568d` | `sintax_analyse` takes `View<int>`; `max_element` + `group_of()` |
| 5 | `835a49c9` | `fill_results()`; the second copy of the loop deleted |
| 6 | `c828770b` | `align_delayed` over `drop()`; the pointer alias dropped too |
| 7 | `36d2ded6` | `column()` accessor pair; `fill_n`/`fill`/`max_element` |
| 8 | `d85df118` | `extra_list` as a `View`; one of three loops a range-for |
| 9 | `d0bc5902` | the four `// refactoring:` briefs carried out and deleted |
| 10a | `8e73eaa9` | DUST's three loops; `wo()` returns a `DustRegion` |
| 10b | `07ebeec1` | the two k-mer counters over `first()`/`drop()` |
| — | `6789d053` | the six clang-tidy findings this campaign introduced |

Census, measured with one command on both trees (`src/`, excluding
`vendored/`):

| | 511ed257 | now |
|---|---:|---:|
| index-style `for (a; b; c)` | 283 | **249** |
| range-for | 89 | **103** |
| `[static_cast<...>(...)]` subscripts | 424 | **404** |

The index-loop count falls by 34 rather than by the ~40 the phase list
implies, because eight loops the plan expected to convert legitimately
keep their counter (below), and because `chimera.cpp` -- 53 of the 283 --
is out by Decision 3.


## Where the plan was wrong

Four places, all found while implementing, none changing the shape of the
campaign:

1. **Phase 1's three commands cannot lose their index.** The plan says the
   named `top_hits()` view "replaces the three copies and the index loop at
   once". It replaces the copies, but each loop also tests `t == 0` to
   decide whether to write a `--uc` record (`--uc_allhits` writes them all),
   which is rule 3's *index selects behaviour*, the same category as the SAM
   secondary flag the plan cites. The loops keep an index -- now a
   `std::size_t` over the named view rather than an `int64_t` paired with a
   separate bound -- and say why.

2. **Phase 2 excludes three sites, not two.** `commands/sintax.cpp:560`
   assigns `si->strand = s` exactly as `commands/search_exact.cpp:489`
   does; the plan named only the latter and `sintax.cpp:440`.

3. **Phase 4's inner `find_if` is not an improvement.** The index the
   inner search finds *is* the result -- stored in `cand_match`, used to
   subscript `cand_matchcount` -- so a `std::find_if` would have to
   recover it from the iterator, or reach for the element's address inside
   the predicate to test `cand_included`. A named lambda with an early
   return (`group_of()`) removes the break from the middle of the nest
   without that. Phase 4's `max_element` landed as planned, including the
   `> 0` guard the plan flagged.

4. **Phase 4's "three loops over `tax_levels`" index four arrays, not
   one.** `:255` and `:280` index `level_best`, `level_matchcount`,
   `taxonomic_fields` and `cand_level_name`'s second subscript in lockstep,
   so a range-for over any one of them cannot serve. They keep an index
   (now `std::size_t`, so the per-iteration cast goes) and gain the
   assertion that `level_best` is not `-1` there -- the invariant that
   keeps that subscript in range, previously unstated. The same is true of
   Phase 8: two of its three `clusterinfo_v` loops drive a `Progress` bar
   with the ordinal, which the plan anticipated.

Two further things the plan did not know:

- **`utils/taxonomic_fields.h` had no include guard.** It only ever worked
  because nothing included it twice; Phase 3a's `tax.hpp` needs
  `tax_levels`, which made every consumer include it twice. Fixed in
  `4c01081c`.

- **The measurements went the other way twice.** Phase 9's `filter.cpp`
  and Phase 7's `msa.cpp` are *faster*, not neutral -- see below.


## Measurements (Decision 4 and verification step 2)

Instruction counts under callgrind, always with `--qmask none --dbmask
none` where DUST was not the subject, and always read against the function
list rather than the total alone.

| phase | benchmark | before | after | |
|---|---|---:|---:|---|
| 5 | `example_search` (library) | 238,394,636 | 238,389,869 | -0.002% |
| 6 | 400 queries x 600 refs | 134,318,123,914 | 134,318,058,085 | -0.00005% |
| 7 | `--cluster_size --msaout` | 79,572,313,201 | 79,568,162,231 | -0.005% |
| 9 | `--fastx_filter`, 400k reads | 26,745,914,103 | 24,262,258,963 | **-9.3%** |
| 10a | `--fastx_mask --qmask dust` | 8,777,348,569 | 8,766,963,331 | -0.118% |
| 10a | ... `--hardmask` | 8,468,600,301 | 8,457,964,473 | -0.126% |
| 10b | `count_bitmap` alone | 191,143,852 | 191,128,732 | -0.008% |
| 10b | `count_hash` alone | 92,551,438 | 92,549,838 | -0.002% |

Decision 4's revert-if-negative was never triggered: every measurement came
back at or slightly below the baseline, so both halves of Phase 10 are kept,
and `wo()`'s signature change (gated on 10a's first half) went in.

No `hyperfine` run: verification step 3 asks for one only after step 2 says
the work changed. Where it did (Phase 9), the cause is visible in the
profile and is not a wall-clock question -- see below.

Two of these numbers are worth more than their sign:

- **Phase 9's -9.3% is an inlining threshold, not an algorithm.**
  `fastq_get_qual` is 2,943,265,182 instructions as a separate function
  before the change and is inlined into `analyse()` afterwards; the
  `count_if` removing code from the caller is what flipped the decision.
  Real on this compiler, not to be banked on. That benchmark is also 42%
  libm `pow()`, from the expected-error computation, not filtering.

- **Phase 7's own share moved much more than its total.** The run is 98%
  SIMD `aligncolumns`; `msa.cpp`'s own instructions fall from 143.3M to
  127.1M, of which `update_profile` alone is 41,095,239 -> 29,972,519
  (-27%) -- the per-nucleotide stride multiply becoming one `subspan` per
  call.

A third is worth recording as a trap avoided: `wo()` is no longer a
separate function in the DUST profile at all -- it is inlined into
`dust_core` after Phase 10a, where it was 81.46% of the run before. Its
individual alignment sensitivity, the thing that produced a spurious 4%
on 2026-08-04, is therefore not comparable across this change. The totals
are, and they are equal.


## Verification, final state

- **Byte-identical output.** 441 comparisons of `byte_identity.sh` and 138
  runs of a campaign matrix written for this branch (`--lcaout`,
  `--top_hits_only`, `--maxhits`, `--uc_allhits`, `--sintax` at four
  cutoffs, the cluster MSA/consensus/profile writers, `--fastq_maxns`,
  `--fastq_convert`, the mask modes and three word lengths), against a
  release build of `511ed257`, both sides at `--threads 1`, after every
  commit.
- **Phase-specific corpora**, where the shared matrix was not sensitive
  enough: 400 real queries against the whole 5,000-sequence PR2 UTAX
  reference across three seeds and three cutoffs for the sintax tie-break;
  `--cluster_size` at three identities for the consensus tie-break; a
  `--clusters` run for the progress-bar ordinal; and a **synthetic
  15-read SFF** built for Phase 9, because the 1.1 GB SFF in
  `~/Documents/data` carries no clip information at all -- every base
  comes out uppercase, so it does not exercise the soft-masking that phase
  rewrites.
- **`vsearch-tests`:** 9,766 PASS, 0 FAIL.
- **`api_examples` on a release build:** 40 PASS, 0 FAIL, and
  `example_search` -- the only example driving both `search_session_single`
  and `search_batch` -- byte-identical.
- **Cross-compiles:** mingw, POWER, RISC-V toolchain, all clean; the only
  warnings are the pre-existing `attributes.cpp` GCC false positive and
  the PPC `-Wpedantic` compound-literal in `align_simd.cpp`, neither of
  which this branch touches.
- **cppcheck:** no new finding on any touched file.
- **clang-tidy**, run on the baseline commit as well and set-diffed:
  1,678 warnings -> 1,587. Six new findings were fixed (`6789d053`);
  eleven remain and are deliberate, the largest group being the category
  shift the plan predicted -- `core/tax.cpp`'s raw `int *` subscripts
  became `std::array` ones, so six reports moved from
  `pro-bounds-pointer-arithmetic` to `pro-bounds-avoid-unchecked-container-access`,
  a check that can see neither `std::array`'s fixed extent nor the debug
  build's assertions. That category still falls by 35 overall.


## Still open

Unchanged from *Out of scope* above, and now with a proven idiom behind
the first of them:

- `core/chimera.cpp` (Decision 3) -- the survey stands, and Phase 7b's
  `column()` accessor is the same shape at a tenth of the scale, so the
  row-accessor idiom that campaign needs is no longer a proposal.
- `core/mergepairs.cpp`, `commands/fastq_eestats2.cpp`.
- `core/derep.cpp`'s report writers, row 2 of the container campaign's
  table: still `(std::vector<struct bucket> const &, uint64_t clusters)`.
- `search16()`'s eight parallel arrays; `Database`'s missing record range.
