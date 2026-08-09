# TBD 2026-08-08 — stop passing std containers as raw pointers

## Status (2026-08-09)

**Phases 1-9 are done**, on `tmp_20260806093624` rather than the
`tmp_<ts>_references` branch the plan proposed (the maintainer asked for the
current branch). Eleven commits, `1f9b22aa`..`e3f17082`. Counted the same way
as the census above: `.data()` 236 -> 188, `.c_str()` 39 -> 14.

**C4 (`search16`) and C9 (`largeread`/`largewrite`) are not done** and still
want the separate branches the plan describes. So do groups D, E and C8, which
this plan does not schedule.

Four places where the plan was wrong or incomplete, corrected in the work:

1. **C6, `core/tax.cpp:167` is not a plain `.data()[i]` -> `[i]` conversion.**
   The level letter can be the header's last byte (`">s;tax=d"`), so the `':'`
   test reads one position *past* the view and worked only because Database
   NUL-terminates its headers. `View::operator[]`'s assert would have fired on
   that input. It now carries an explicit bounds test, which returns the same
   answer without the overread. This is trap 3, at a site the plan listed as
   safe.
2. **C6, the two digest writers must not use `make_view(hex_digest)`.** The
   array is `(2 * digest_length) + 1` bytes and its last byte is the `'\0'` the
   hex writer appends, so the whole-array view would have written a NUL into
   `--samout`'s `M5:` field and into `--relabel_md5`/`--relabel_sha1`. They use
   `.first(size() - 1)`.
3. **C7, `sff_convert`'s `write_report` keeps its `std::fputs`.** `index_kind`
   is a fixed 8-byte buffer whose *content* can be shorter: it is all-zero when
   the index block was never reached ("SFF index missing"), and the eight bytes
   are read raw from the file. A `View` over it would print eight NUL bytes in
   a reachable case. The parameter becomes the `std::array` by const reference,
   which removes the two `.data()` call sites without changing what is written.
4. **C2's `fastq_fatal` and `warn` are converted, not overloaded.** An overload
   pair would have had to duplicate a one-expression body or forward to it and
   allocate anyway; and `warn`, moving to internal linkage with a single caller,
   would have been left with a dead `char const *` form. `set_deferred_error`
   does get the overload the plan describes, with the shared body in a private
   `record_deferred_error(View<char>)`.

One thing not in the plan, added at the maintainer's request during the work:
**`derep()` no longer takes an input file name** beside the `Derep_mode` that
already determines it (`a84803a5`).

Verification actually run: debug build clean at every phase; `cppcheck` and
`clang-tidy` on the touched files with no new findings (and
`bugprone-easily-swappable-parameters` gone from `search_findbest2_*`);
`vsearch-tests` 9882 lines / 0 FAIL, test-name set identical to a baseline run
of the same suite against the pre-campaign binary; `api_examples` 40 PASS / 0
FAIL against a release `libvsearch.a`; the three cross-compiles with warning
sets identical to the pre-campaign baseline; every output of 29 command
invocations byte-identical to the pre-campaign binary (including the `.udb`
binary itself); and callgrind instruction counts within 0.005% on
`usearch_global`, `cluster_fast` and `uchime_denovo`, and -0.016% on
`dust_core`'s own benchmark.

One measurement note for whoever reads a `--log` diff: the "Max memory" line
varies by 0.1-0.2 MB between runs of the *same* binary, so it is not a
comparable field.

---


*"Avoid passing std data structures as raw pointers. Search and check if
calls to `.data()` and `.c_str()` could be replaced."*

Census date: 2026-08-08, branch `tmp_20260806093624` (HEAD `94d4d856`).

Companion to two things that finished the same day:

- **`DONE_20260808_ptr_parameter.md`** — the `&parameter` campaign. That one
  removed `T *` parameters that meant "one object, always there". This one is
  its container-side twin: `T *` (or `char const *` + a separate length) that
  means "the storage of a container the caller is holding".
- **swarm's T5** (`~/Documents/src/swarm/TBD_20260806_misc.md:593`, five
  commits `b9004de`..`bb51303`, **done**). Same issue statement, same week.
  This plan deliberately reuses swarm's A–E group taxonomy so the two trees
  read the same way, and carries over its lessons in "What swarm learned"
  below.

It also **owns group D1 of the `&parameter` plan** — `search16`'s array
parameters — which that campaign deferred with "still worth doing, D1 most of
all". See group C4.


## What this is really about

`.data()` and `.c_str()` at a call site are the *symptom*. The cause is a
callee that declares `char const *` (or `T *` plus an `int len`) where the
caller already holds a `std::vector`, `std::array`, `std::string`, `Span` or
`View` that knows its own extent. The fix is on the **declaration** side, and
the `.data()` then disappears on its own.

So, as with the `&parameter` plan, this document is organised by *callee*.

**The naive reading of the issue — sweep every `.data()` — would break code.**
230 `.data()` sites and 38 `.c_str()` sites do not all mean the same thing;
roughly a third of them are the abstraction working correctly, and a handful
are load-bearing in ways a sed would silently destroy. The split is therefore
the first job, and the "keep" groups are recorded here so a later sweep does
not re-litigate them.


## Census

`.data()`: **236** grep hits, **230** in code (6 in comments).
`.c_str()`: **39** grep hits, **38** in code (1 in a comment).

| group | what | sites | disposition |
|---|---|---|---|
| **A** | inside the View/Span/buffer/print abstraction itself | ~26 | keep — this *is* the abstraction |
| **B** | C / OS / third-party API boundaries | ~41 | keep — the pointer is the contract |
| **C** | vsearch's own functions taking raw pointers | **~91** | **the target** |
| **D** | pointer arithmetic over a container's storage | ~21 | 7 candidates; 14 are a measured design and stay |
| **E** | deliberate, and would be a regression to "fix" | ~14 | keep, and say why |

Group C is the work. It splits by callee into nine sub-groups, C1–C9, in
increasing order of risk. (C9, the UDB reader and writer, was group D at
census time and moved on the 2026-08-09 decision — see "Decisions".)


## Group A — keep: the abstraction's own plumbing

These implement `data()`, or consume it exactly one level below the
abstraction. Converting them is circular.

`utils/span.hpp:266` (`make_span`), `utils/view.hpp:251` (`make_view`),
`utils/print_view.hpp:104` (the `fwrite` inside `fprint(FILE *, View<char>)`),
`utils/print_record.hpp:136,178`, `utils/decimal_digits.hpp:174,195,213,233`,
`core/fastx.hpp:99,100,296,393,396`, `core/db.hpp:173,178,188`,
`core/db.cpp:143,172,483,507,531`, `core/bitmap.cpp:81`,
`core/dbindex.cpp:84,98`.


## Group B — keep: C, OS and third-party boundaries

Identical to the `&parameter` plan's group C1, and to swarm T5's group B. Not
ours to change:

`std::fread`, `std::fputs`, `std::fgets`, `std::snprintf`, `std::strftime`,
`std::strlen`, `std::memcpy` / `memmove` / `memset`, `std::char_traits::find`,
`getopt_long_only` (`cli.cc:3093`), `regexec` (`core/otutable.cpp:178,215,242`
— and the `pmatch_*.data()` arguments beside them), `SHA1` / `MD5`
(`utils/sequence_digest.cpp:128,130,158`), `CityHash64` / `CityHash128`
(`utils/cityhash.cpp:69,75`), and the eight `reinterpret_cast`s onto SIMD
types in `core/align_simd.cpp`.

One note for the reviewer: `std::string::assign(ptr, len)` and
`std::string(ptr, len)` are *also* pointer-taking C++ APIs, but unlike the
above they have an iterator-pair form, so they land in C6 rather than here.


## Group C — the target, by callee

### C1. `fatal()` — 13 sites that are a pure deletion

`fatal(std::string const &)` **already exists** (`utils/fatal.hpp:122`, added
by the C-style-output campaign) and its body is `fatal(message.c_str())`
(`utils/fatal.cpp:87`). Every caller that writes `fatal(message.c_str())` is
therefore taking the length off a string so that the overload it just selected
can measure it again with `std::strlen`.

Drop the `.c_str()`; nothing else changes.

| site | form |
|---|---|
| `parameters.cpp:92`, `:181` | named local `message` |
| `cli.cc:298`, `:4633` | named local `message` |
| `utils/open_file.cpp:132` | named local `message` |
| `commands/fastq_stats.cpp:137` | named local `message` |
| `core/fastq.cpp:235` | named local `message` |
| `core/eestats.cpp:78`, `:86` | temporary `("…" + std::to_string(x) + ")")` |
| `core/filter.cpp:93`, `:101` | temporary, same shape |
| `commands/fastq_mergepairs.cpp:147`, `:155` | temporary, same shape |

`utils/fatal.cpp:87` is the overload's own body and stays.

swarm hit the same thing at `db.cpp:759` and the disposition was the same — a
deletion, not a signature change (`0a362a9`).

### C2. the deferred-error / warning trio — 9 sites, one overload each

Same shape as C1, but the `std::string` overload does not exist yet. Each of
the three is a two-line addition:

| callee | declaration | `.c_str()` sites |
|---|---|---|
| `fastx_s::set_deferred_error` | `core/fastx.hpp:313` | `core/fasta.cpp:152,165`; `core/fastq.cpp:232`; `core/fastx.cpp:285,305,756` |
| `warn` | `core/fastx.cpp:248` | `core/fastx.cpp:328` |
| `fastq_fatal` | `core/fastq.cpp:219` | `core/fastq.cpp:429,546` |

`warn()` is defined at file scope in `core/fastx.cpp` with external linkage but
is declared in no header, and its only caller is 80 lines below it in the same
TU. **Decided 2026-08-09: the same commit moves it into that file's anonymous
namespace (`:104-116`)**, which is the linkage a single-caller, header-less
function should have.

Verify before moving, not after — the whole argument rests on there being no
other reference:

```sh
nm -C src/core/*.o src/commands/*.o src/utils/*.o | grep -w warn
```

Expect the definition as `T` in `fastx.o` and no `U` anywhere. If any TU
references it, stop and leave the linkage alone — the overload is still worth
adding on its own.

`set_deferred_error` is worth doing for more than the call sites.  Its body
(`core/fastx.cpp:182-195`) does a truncating copy into
`std::array<char, 512> errmsg` and computes the length with
`std::strlen(message)`. A `std::string const &` overload uses `message.size()`
instead — no `strlen`, and correct in the presence of an embedded `'\0'`. The
four literal callers (`core/fasta.cpp:287,305`, and two more) keep selecting
the `char const *` overload, exactly as the `fatal()` comment at
`utils/fatal.hpp:118-120` describes.

### C3. `Progress`'s prompt — 4 sites, and an undocumented lifetime requirement

`Progress(char const * prompt, …)` (`utils/progress.hpp:71`) **stores the
pointer** (`prompt_`) and re-prints it later (`:87`). Four callers hand it a
`std::string`:

`commands/derep_smallmem.cpp:258`, `core/udb.cpp:249`, `core/derep.cpp:569`,
`core/db.cpp:295`.

**Decided 2026-08-09: the member becomes `std::string prompt_` and the
constructor takes `std::string` by value.** This deletes the undocumented
lifetime requirement rather than documenting it. It costs one small allocation
per `Progress` — built once per *phase*, never per record, so unmeasurable —
and the 58 literal sites gain a temporary without changing in source.

Rejected: a `std::string const &` overload that still stores `.c_str()`. It
would remove the four call-site `.c_str()`s while making the
dangling-temporary case *easier* to reach.

### C4. `search16`'s array parameters — the `&parameter` plan's D1

`core/align_simd.hpp:101-110`:

```cpp
auto search16(s16info_s * searchinfo,
              unsigned int sequences,
              unsigned int const * seqnos,
              CELL * pscores,
              unsigned short * paligned,
              unsigned short * pmatches,
              unsigned short * pmismatches,
              unsigned short * pgaps,
              std::string * pcigar,
              struct Database const & db) -> void;
```

Eight parallel arrays whose extent is carried by `sequences`, five of them
same-typed and adjacent — the textbook swappable-parameter hazard. Four
callers:

| caller | form |
|---|---|
| `core/searchcore.cpp:787-796` | `.data()` × 7 on real vectors, count = `target_count` |
| `core/chimera.cpp:2167-2176` | `.data()` × 7 on `chimera_info` members |
| `commands/allpairs_global.cpp:409-418` | `.data()` × 7 on real vectors, count = `searchinfo.hit_count` |
| `core/cluster.cpp:758-767` | `1` and seven `& scalar` |

**The detail that matters, and that a naive conversion gets wrong:**
`sequences` is a *used* count, not the container size. `searchcore.cpp` passes
`target_count`, having just filled the first `target_count` entries of
tophits-sized vectors. So the conversion is **not** `make_span(v)` — it is
either `make_span(v).first(count)` at every call site, or keep `sequences` and
`assert(span.size() >= sequences)` in the callee. The former is better (one
source of truth for the extent) but it makes `cluster.cpp`'s N=1 caller build
seven one-element spans; that caller is inside the clustering hot loop.

`search16_qprep(s16info_s *, View<char> qseq)` (`:98`) is already in the new
style, so the file's own precedent points this way.

**This is the largest and riskiest item in the plan and should be its own
branch.** Sequencing note below.

### C5. per-thread `searchinfo_s` arrays — 21 sites, and a latent bug

The dominant pattern in the search commands:

```cpp
struct searchinfo_s * const si_plus  = state.si_plus.data();
struct searchinfo_s * const si_minus = state.si_minus.empty() ? nullptr : state.si_minus.data();
…
search_thread_init(si_plus + t, …);
if (si_minus != nullptr) { search_thread_init(si_minus + t, …); }
```

`state.si_plus` is a `std::vector<searchinfo_s>` with one element per thread.
The code flattens it to a pointer, re-derives element `t` by pointer
arithmetic, and encodes "no reverse strand" as a null pointer — three things
the vector already says better.

| file | sites |
|---|---|
| `core/cluster.cpp` | `:359,362` (`.data() + work.query_first + q`), `:896,897`, `:1829,1830` |
| `commands/usearch_global.cpp` | `:375,376`, `:424`, `:456`, `:489,490` |
| `commands/sintax.cpp` | `:459,460`, `:576,577`, `:692,693` |
| `commands/search_exact.cpp` | `:733,738` |
| `core/chimera.cpp` | `:2527` |

The callees mostly already take a *reference* to one `searchinfo_s` (the
`&parameter` campaign did that: `cluster_query_init`, `query_exit`), so most
sites become `state.si_plus[t]` with no signature change at all — the `.data()`
and the `+ t` both vanish.

Where a callee genuinely needs the whole array, `Span<searchinfo_s>` replaces
the pointer and `.empty()` replaces the `!= nullptr` test.

**Two sites in this group are not a style question.**

1. **`commands/search_exact.cpp:733,738` and `core/chimera.cpp:2527` store
   `.data()` of a *local* vector into a longer-lived state struct**
   (`state.si_plus = si_plus_v.data();`, `state.cia = cia_v.data();`). The
   pointer members outlive nothing today because the vectors happen to be
   declared in the enclosing scope, but the type says nothing about it —
   whereas `usearch_global.cpp` and `sintax.cpp` hold the *vector* in the
   state struct. Two commands, two shapes, one of them fragile. Make
   `search_exact` and `chimera` match the other two.

2. **`si_m`'s nullability is not checked by anything — it is implied by a
   second variable, reached through a third object.** This is the finding of
   the census, and it is not a style point.

   `search_findbest2_bysize` / `search_findbest2_byid`
   (`core/searchcore.cpp:964`) and `free_hit_alignments`
   (`core/cluster.cpp:864`) each take `struct searchinfo_s * si_m`. Their four
   call sites disagree about what that pointer is:

   | call site | what `si_m` is |
   |---|---|
   | `core/cluster.cpp:973,977,1032` | `opt_strand ? si_minus + i : nullptr` (`:958`) — **a real null** |
   | `core/cluster.cpp:1890,1894,1934` | same ternary (`:1874`) — **a real null** |
   | `core/cluster.cpp:1786` | `cs->si_minus.get()` — **null when the `unique_ptr` is empty** |
   | `core/cluster.cpp:1085,1089,1126` | `si_m.data()` on a `std::array<searchinfo_s, 1>` (`:1049`) — **never null, but the object is unusable**: `cluster_query_init` only runs on it when `opt_strand` (`:1052-1055`) |

   **Neither callee tests the pointer.** `free_hit_alignments` builds
   `std::array<searchinfo_s *, 2>{{si_p, si_m}}` and iterates
   `.first(number_of_strands(parameters.opt_strand))`;
   `search_findbest2_byid` guards its `si_m` walk with
   `if (parameters.opt_strand)` — and reads that `parameters` out of
   `*si_p->parameters` (`core/searchcore.cpp:967`). So the invariant is

   > `si_m` is dereferenceable **iff** `si_p->parameters->opt_strand`

   spread across two parameters and one object, asserted nowhere. It holds
   today. Nothing in the signature says it must.

   **This rules out the obvious fix.** A `searchinfo_s &` is not formable at
   three of the four call sites (there is no object), and at the fourth it
   would bind to an object that exists but was never initialised. The correct
   shape is **`Span<searchinfo_s>`** — empty means absent, the callee can
   `assert(si_m.empty() != parameters.opt_strand)`, and the `std::array<_, 1>`
   caller passes `make_span(si_m)` while the ternary callers pass
   `Span<searchinfo_s>{}`. That is a real change to a threaded hot path and
   should be its own commit within phase 6.

   **Record for whoever extends this:** `std::vector::data()` on an empty
   vector is *unspecified* — it may or may not be null. Any future site that
   uses `container.data() == nullptr` as an "absent" signal is relying on
   libstdc++ behaviour, not on the standard.

### C6. `View` / `Span` flattened back to `(pointer, length)` — ~15 sites

The callee already accepts, or trivially could accept, the sized type. These
undo work the View/Span propagation campaign (`DONE_20260728`) already did.

| site | current | replacement |
|---|---|---|
| `commands/fastx_syncpairs.cpp:184` | `std::string {header.data(), header.size()}` | `std::string {header.begin(), header.end()}` |
| `commands/fastx_syncpairs.cpp:207,208,210` | `record.x.assign(v.data(), v.size())` | `assign(v.begin(), v.end())` |
| `vsearch.cc:167` | `std::string(cores.data(), cores.size())` | `decimal::to_text(system_get_cores())` — the helper already exists (`utils/decimal_digits.hpp:210`) |
| `cli.cc:4629,4631` | `std::string(limit.data(), limit.size())` | `decimal::to_text(…)`, same |
| `core/tax.cpp:110,161,167` | `header.data()[offset]` | `header[offset]` — `View::operator[]` exists and asserts |
| `core/tax.cpp:84` | `header.data() == nullptr` | `header.empty()` — **but see trap 2** |
| `commands/sintax.cpp:466` | `si_plus[t].query_head.data()` | keep the `View<char>` |
| `core/chimera.cpp:641` | `std::fputs(ci->query_head.data(), stdout)` | `fprint(stdout, ci->query_head)` |
| `core/results.cpp:994` | `std::fputs(md5hex.data(), …)` | `fprint(…, make_view(md5hex))` |
| `utils/sequence_digest.cpp:175,183` | `std::fputs(hex_digest.data(), …)` | same |
| `core/otutable.cpp:341` | `char const * otu_name = it_otu.c_str();` then `std::fputs` | `fprint(output_handle, make_view(it_otu))`, local goes |
| `core/attributes.cpp:219` | `std::next(header.data(), …)` | `header.subspan(…)` |
| `core/searchcore.cpp:370` | `hit.nwalignment.c_str()` | see C8 (cigar) |

`fprint(std::FILE *, View<char>)` already exists at `utils/print_view.hpp:90`
and `fwrite`s the exact extent, so the `fputs` conversions also drop a
dependency on a `'\0'` being present.

### C7. one-off callees — 8 sites

| callee | site | note |
|---|---|---|
| `dust_core(char *, int, bool)` | `core/mask.cpp:223` | `dust()` already takes `Span<char>`; `dust_core` is `static` with two callers. **Hot path — see trap 5.** |
| `open_output_file(char const *)` | `core/cluster.cpp:1427` (+ `:1431` in the error message) | **Decided: build the name as a `std::string`** (`opt_clusters + decimal::to_text(clusterno)`) and keep `.c_str()` at the `std::fopen` boundary as an honest group-B site. This also retires the hand-rolled `std::copy` + manual `'\0'`, the `space_for_cluster_id` sizing dance, and the second `.data()` at `:1431`. Not a `View<char>` overload — `std::fopen` needs the terminator, so the overload would rebuild the string anyway. |
| `fasta_print(…, char const * header, char const * seq, uint64_t len, …)` | `core/msa.cpp:510` | **one caller in the whole tree.** Its body (`core/fasta.cpp:389-399`) already reassembles `View<char>{seq, len}` internally and `fputs`es the header, so taking two Views deletes code rather than adding an overload. The rest of the family (`fasta_print_general`, `core/fasta.hpp:90`) is already on Views. |
| `write_report(…, char const * index_kind)` | `commands/sff_convert.cpp:700,704` | `static` in the same file, 2 callers |
| `random_shuffle(T *, std::size_t, URBG &)` | `commands/shuffle.cpp:93` | `utils/random.hpp:151`; **1 caller in the tree**, passes `deck.data(), deck.size()`. `Span<T>` is a strict improvement and the cheapest commit in the plan. |
| `std::string = vector.data()` | `core/derep.cpp:1029` (`bp->seq = ds->seq_up.data()`) | an implicit `strlen` over the sequence, per unique record, in the dereplication insert path. The correctly-sized view is already in scope as `seq_up_v` (`:635`). Possibly a small *win*. |
| `.assign(ptr, len)` | `core/mergepairs.cpp:780,781` | iterator-pair form |

### C8. the cigar walkers — coordinate, do not duplicate

`LinearMemoryAligner::alignstats(char const * cigar, …)`
(`core/linmemalign.hpp:187`) is fed `.c_str()` at four sites
(`core/chimera.cpp:2193`, `core/cluster.cpp:782`,
`commands/allpairs_global.cpp:444`, `core/searchcore.cpp:834`), and
`LinearMemoryAligner::align()` returns `cigar_string.data()`
(`core/linmemalign.cpp:697`). `core/searchcore.cpp:370,694` take
`hit.nwalignment.c_str()` for the same reason.

`alignstats` was **just** reshaped by the `&parameter` campaign (phase 4,
`66dd6479`, five out-params → an `AlignStats` return), and
`TBD_20260808_msa_cigar.md` is an open plan about the two cigar-walking APIs.
Three plans touching one function is how comments go stale.

**Decided 2026-08-09: C8 is out of this campaign.** Fold the
`char const * cigar` → `View<char>` question into `TBD_20260808_msa_cigar.md`,
which is already the document about how cigars are walked. Note in passing
that `core/showalign.cpp:385` already carries the comment *"cigar string can be
trimmed (left and right): `cigar.size()` maybe != `std::strlen(cigar.data())`"*
— i.e. the codebase has already written down that the pointer form loses
information here. That is the argument for the View, and it belongs in that
plan.


### C9. `largeread` / `largewrite` — 16 sites, the UDB binary format

**Decided 2026-08-09: convert, with the extent spelled explicitly at every
call site, on its own branch.** The census had recommended leaving these
alone; the deciding argument the other way is that `Span`'s constructor
asserts, and a binary format parser is exactly where a bounds assert earns its
keep. It is the highest-risk item in the plan after C4 — the only one that can
produce a file that looks valid and is not — so it does not ride along with
the `.c_str()` deletions on the main branch.

```cpp
auto largeread (std::istream & input,  void       * buf, uint64_t nbyte, …);  // core/udb.cpp:100
auto largewrite(std::ostream & output, void const * buf, uint64_t nbyte, …);  // commands/makeudb_usearch.cpp:82
```

**The parameter cannot be `Span<char>`:** half the callers pass
`std::vector<unsigned int>` or `std::array<unsigned int, 50>`. It has to be a
template on the element type, taking the byte count from `size_bytes()`:

```cpp
template <typename Type>
auto largeread(std::istream & input, Span<Type> buf, uint64_t offset, Progress & bar) -> uint64_t;
```

Call sites then spell the extent in **elements**, not bytes, and the pervasive
literal `4` disappears.

**Classify every site before touching it.** The 16 split three ways, and the
middle group is where the risk lives:

| kind | sites | note |
|---|---|---|
| full range (`nbyte == size_bytes()`) | `udb.cpp:296,327,419`; `makeudb:185` | `kmercount`/`kmerindex` are `resize()`d on the line above; `sequence_lengths(seqcount)` reads `4 * seqcount` |
| **sub-range** | `udb.cpp:250,316,345,377`; `makeudb:182,189,205,234,243,257` | must be `.first(n)` |
| **no container at all** | `udb.cpp:408,455`; `makeudb:211,249,263` | three different cases, all of them expressible — see below |

**The five "no container" sites, resolved (decided 2026-08-09).** They looked
like one problem and are three, and none of them needs the second overload the
census assumed:

| site | becomes | what it gains |
|---|---|---|
| `makeudb:211` | `make_span(dbindex.kmerindex).subspan(kmerhash[i], kmercount[i])` | it *was* a container all along — a mid-vector offset. `subspan` asserts the offset; `+ kmerhash[i]` does not, and `kmerhash[i]` is an index value, not a constant. |
| `makeudb:263` | `db.sequence_view(i)` | already an exact match: `View<char>{getsequence(i), getsequencelen(i)}` |
| `udb.cpp:408,455` | one `Span<char>{db.data_.data(), datap_bytes}` formed after `udb_reserve`, then `.first(udb_headerchars)` and `.subspan(udb_headerchars, nucleotides)` | **the bounds assert lands where the length comes from the file.** `udb_read` is already a `friend` (`core/db.hpp:122`), so no accessor is needed. |
| `makeudb:249` | `db.header_view(i)`, then the terminator written separately | see below |

**`makeudb:249` — write the terminator, do not reach for it.** The site writes
`len + 1` bytes from `db.getheader(i)`: the header, plus a `'\0'` that it reads
*out of `Database`'s buffer*. That byte is real — `Database::add` appends it
(`core/db.cpp:222-226`, "the views need not be NUL-terminated, so add()
appends the terminator itself") — but reading it here couples the UDB writer to
`Database`'s internal layout.

The terminator is a property of the **UDB format**, not of `Database`
("headers (ascii, zero terminated, not padded)", `:245`). So it is written at
the write site, and `db.header_view(i)` is used unchanged. **No
`header_view_z()` accessor** — that would have made the coupling official
instead of deleting it.

Route the single byte through `largewrite` rather than calling
`out_stream.put()` directly, so it keeps both things the split would otherwise
drop: the `if (not output) { fatal(…); }` check
(`commands/makeudb_usearch.cpp:91-94`) and the `pos +=` accounting that every
later offset depends on.

```cpp
static constexpr char terminator = '\0';
…
pos += largewrite(out_stream, db.header_view(i), pos, progress_bar);
pos += largewrite(out_stream, View<char>{&terminator, 1}, pos, progress_bar);
```

Output is byte-identical: `out_stream` is a buffered `std::ofstream` (`:114`)
and this format section has no padding or alignment.

**Two sub-range sites are the argument for the whole phase:**

- **`core/udb.cpp:377`** reads `uint64_t{4} * seqcount` bytes into
  `std::vector<unsigned int> header_index(seqcount + 1)` declared two lines
  above (`:375`), because element `[seqcount]` is filled by hand afterwards
  (`:379`). `make_span(header_index)` would read one element too many from the
  UDB file and desynchronise every subsequent `pos`. Written as
  `.first(seqcount)` next to a `seqcount + 1` declaration, the asymmetry is
  visible instead of hiding inside a byte multiplication.
- **`commands/makeudb_usearch.cpp`** reuses **one** `std::vector<unsigned int>
  buffer(buffersize)` (`:168`) as scratch for six different fields
  (`50*4`, `1*4`, `4*elements`, `4*8`, `4*seqcount` twice). Six sites where
  `make_span(buffer)` would write the whole scratch vector to the output file.

**One thing the conversion must not hide.** The literal `4` is the UDB
format's field width, which is 4 bytes *by definition of the format*, not
because `sizeof(unsigned int)` happens to be 4. `size_bytes()` makes that
implicit, so the phase must add
`static_assert(sizeof(unsigned int) == 4, "UDB stores 32-bit fields");` next to
the reader and the writer. Without it, the format's contract survives only as
an accident of the platform — and the cross-compiles would not catch it,
since all four targets have 32-bit `int`.

**Verification for this phase specifically:** a UDB round trip
(`makeudb_usearch` → `udb_info` → a search against the `.udb`), byte-compared
against the reference binary, plus reading a UDB written by the *pre-change*
binary and vice versa. A byte count that shifts by one element produces a
valid-looking file that fails somewhere else entirely.

## Group D — pointer arithmetic over a container's storage

Some of these want a `Span`; several are a measured design and must not be
swept. Splitting them is a per-site judgement, so this group is listed for
completeness and scheduled last, or not at all.

**Candidates for a `Span`:**
`core/fastq.cpp:252`, `core/fasta.cpp:199,257`, `core/mask.cpp:204`,
`commands/makeudb_usearch.cpp:212`, `core/chimera.cpp:1333,2030`.

**`largeread` / `largewrite` — moved to group C on the 2026-08-09 decision.
See C9; it gets its own branch.**

**Keep — measured design (group E overlap):**
`core/unique.cpp:121,122,204,205,273,286` hoist `bitmap_.data()` /
`hash_.data()` / `list_.data()` into locals *outside* a hot loop, and
`core/align_simd.cpp:1072,1150,1248,1405,1438-1440,1535,1579` do the same for
the SIMD kernels. swarm measured this exact shape (`80df510`): inlining the
`.data()` at each use cost two prologue instructions and 16 bytes of
`search8`. These are the abstraction being paid for once, correctly.


## Group E — keep, and say why

1. **`core/derep.cpp:1073,1074`** — `results[count].header = b.header.c_str();`
   The `results` array is the **public library API** (`vsearch_api.h`); its
   fields are `char const *` by contract and the C consumer needs the
   terminator. Changing them is an ABI break for no gain. Keep.

2. **`core/otutable.cpp:164,165`** —
   `assert((h.data() == nullptr) or (h.data()[h.size()] == '\0'))`. This reads
   **one past** the view's extent on purpose: it is the assertion that the
   header is NUL-terminated, which is what lets `regexec` see it. `h[h.size()]`
   would trip `View::operator[]`'s own bounds assert. Keep the `.data()`
   spelling.

   **Comment decision (2026-08-09): this one gains a sentence.** Nothing
   currently says the overread is deliberate, and `h[h.size()]` is exactly what
   a later sweep would reach for. The added comment states the intent and names
   `regexec` as the reason; it does not modify or remove anything already
   there. The other two comment candidates stay **verbatim**:
   `core/showalign.cpp:385` (it becomes the argument for the View if
   `TBD_20260808_msa_cigar.md` converts the cigar walkers, and that plan may
   cite it) and `core/getseq.cpp:200` (the recorded evidence for trap 2).

3. **`utils/maps.cpp:373,378,383,388,393,398,403,408`** — the eight
   `chrmap_*()` accessors return `.data()` of a file-static `std::array`. These
   are the 256-entry lookup tables indexed in the innermost loops of the
   parsers and the masker, and **`TBD_20260714_maps_perf_fix.md` exists because
   the last time these were touched it cost 15-20% on two hot loops** (vsearch
   is built without LTO, so every accessor is a real cross-TU call). Whatever
   happens to them is that document's business, not this one's. Keep. The same
   applies to `core/fastq.cpp:413,524`, which pass those tables onward.

4. **`core/getseq.cpp:200`** already carries the comment *"`label.data()` return
   nullptr, crashing downstream scanners"* — a recorded instance of trap 2
   below. Keep the comment; it is evidence, not clutter.

5. **`core/fasta.cpp:268,270` and `core/fastq.cpp:333,335,337,339`** —
   `buffer.data()[0] = 0` writes the terminator into a `FastxBuffer` whose
   *length* is 0 but whose storage is allocated. `buffer.front()` would be UB
   (or a `_GLIBCXX_DEBUG` abort) precisely when it matters. Keep unless
   `FastxBuffer` grows a `clear()` that expresses the intent.


## A finding this census turned up — the copy-and-terminate idiom

Not a `.data()` or `.c_str()` site, so not group C. It surfaced from the C9
decision at `makeudb_usearch.cpp:249` (write the terminator, do not read it out
of someone else's buffer) and it turns out the tree had already reached that
conclusion three times, independently, without sharing the result.

**The helper exists** — `copy_and_terminate()`, `core/chimera.cpp:296-302`, in
chimera's anonymous namespace, used five times inside that file
(`:2275,2276,2296,2297,2972`):

```cpp
auto copy_and_terminate(View<char> const source,
                        std::vector<char> & destination) -> View<char> {
  assert(destination.size() > source.size());
  std::copy(source.cbegin(), source.cend(), destination.begin());
  destination[source.size()] = '\0';
  return make_view(destination).first(source.size());
}
```

**Nine sites hand-roll it.** All nine target a `std::vector<char>`, so
promoting the helper needs no template:

| file | sites |
|---|---|
| `core/search.cpp` | `:115-118` header, `:123-125` sequence (plus strand only — the minus branch terminates inside `reverse_complement`) |
| `commands/sintax.cpp` | `:619-622` header, `:623-625` sequence |
| `commands/search_exact.cpp` | `:512-515` header, `:516-518` sequence |
| `commands/fastq_mergepairs.cpp` | `:534-535` fwd header, `:541-542` rev header |
| `core/chimera.cpp` | **`:2975-2976`** — the file not using its own helper, one line after using it |

`sintax.cpp` and `search_exact.cpp` are byte-identical modulo `si_plus[t]`
versus `state.si_plus[t]`.

**Three comments in three files already say the same thing**, which is the
actual finding — the convention is settled, only its implementation is not:

- `core/chimera.cpp:290-293` — "the terminator is written here rather than read
  from the source: the byte after a header or sequence view is a `'\0'` in the
  reader's and the Database's buffers alike, but it is not part of the view
  being copied"
- `core/search.cpp:111-114` — "the callers all happen to pass a NUL-terminated
  header, but that byte is outside the view they hand over"
- `core/chimera.cpp:2973-2974` — "rather than reading a
  `query.sequence.size() + 1`'th byte the view does not cover"

**Scope, if this is taken up.** Promote `copy_and_terminate` to `utils/`
(a new `utils/copy_terminate.hpp`; it does not belong in `view.hpp`, which is
about the view abstraction and not about C-string copies) and use it at the
nine sites. The `resize()` that precedes three of them stays at the call site —
it is a headroom decision (`buffer_headroom`), not part of the copy. Chimera's
sibling helper `copy_label()` (`:310-317`, the fixed-`std::array` form with
truncation, seven users) would move with it or stay; decide when the first one
moves.

**Do not fold in:**

- `utils/string_normalize.cpp:78` and `utils/reverse_complement.cpp:81` are
  *producers* (`std::transform`), not copies. Both already carry the same
  `assert(dest.size() > src.size()); // room for the '\0' terminator` prologue
  — a third and fourth independent statement of the convention — but each is
  one well-named function and merging them would gain nothing.
- `commands/sff_convert.cpp:630,631`, `core/showalign.cpp:308-324`, and the
  `buffer[0] = 0` resets in `core/fasta.cpp:268,270`, `core/fastq.cpp:333-339`,
  `core/fastx.cpp:124,155` are **truncation**, not append. Different operation,
  same-looking line.

**One asymmetry to look at, reported as an observation and not a defect.**
`commands/fastq_mergepairs.cpp:534-563` copies six views — fwd/rev ×
header/sequence/quality — and terminates only the two headers. Those two carry
`// fix issue when reusing allocated mem`, which records a past bug where a
stale tail was read as a C string; the four sequence/quality copies are
presumably always consumed by length. That may well be correct. It is also
exactly the kind of 2-of-6 split a shared helper makes visible and forces
someone to justify, which is most of the value here. **Preserve that comment
verbatim** if the site is touched.

**Sequencing.** Independent of every phase in this plan, and of C4 and C9. It
should not be bundled with them: this campaign is about call-site spellings,
and this one adds a shared function. Its own commit, or its own small branch,
whenever the maintainer wants it.


## Traps

Six things that would make a plausible-looking conversion wrong. Numbers 1-3
are specific to this campaign; 4-6 are the standing ones.

1. **`.data()` on an empty `std::vector` is unspecified**, so
   `container.data() == nullptr` is not a portable "absent" test. `core/cluster.cpp`
   currently relies on it and is saved only by the container being a
   `std::array<_, 1>`. See C5.

2. **`nullptr` and "empty" are not interchangeable in this tree, and at two
   sites the difference is output-visible.**

   A `View<char>` has three states, not two: **absent** (`{nullptr, 0}`, what
   the default constructor gives), **present but empty** (`{valid, 0}`), and
   present. `View`'s constructor permits all three — it asserts only
   `(start != nullptr) or (length == 0)`. `FastxBuffer::view()`
   (`core/fastx.hpp:104`) returns `View<char>{storage_.data(), length}`, and
   `storage_` is an allocated `std::vector`, so **a record with an empty header
   (`>` alone on the line) yields a non-null, zero-length view** — state two,
   not state one.

   Six sites branch on the pointer. They are not the same case:

   | site | verdict |
   |---|---|
   | `core/otutable.cpp:169` `has_sample = (query_header.data() != nullptr)` | **load-bearing.** `commands/search_exact.cpp:822` and `commands/usearch_global.cpp:756` pass a literal `View<char>{}` to mean "no sample — this is the db-only pass", while `core/cluster.cpp:443,556` pass a real header. Switching to `.empty()` would reclassify every record with an empty header as having no sample, silently changing `--otutabout`. |
   | `core/otutable.cpp:206` `has_otu` | same shape, same verdict |
   | `core/chimera.cpp:2955-2962` | **already correct, do not "simplify".** It tests `.empty()` and `.data() == nullptr` as two checks with two different `fatal()` messages, and the comment above (`:2947-2954`) says why: the constructor's assert is stripped under `NDEBUG` and this is a library entry point. It is the site that proves the distinction is deliberate. |
   | `core/tax.cpp:84` | **probably safe** — an empty header cannot contain `tax=` either, so both forms return `false`. Convert only with that reasoning written down. |
   | `core/attributes.cpp:139` | **probably safe**, same argument for `size=`. Note the second half of the condition, `attribute.text == nullptr`, is a genuine `char const *` null check and stays regardless. |
   | `core/fastx.hpp:271-275` `get_quality()` | **contract.** A FASTA record reports `nullptr`, and the comment says so. Keep. |

   `core/getseq.cpp:196-201` is the same collision seen from the other side,
   and it is already commented: an empty `std::vector<char>` makes `.data()`
   return `nullptr`, so the loop skips empty labels rather than storing one
   that would read as absent downstream.

   **Rule for this campaign:** before replacing any `x.data() == nullptr` with
   `x.empty()`, find the caller that produces the absent case and check whether
   a present-but-empty value can reach the same branch. Where they differ, the
   `.data()` test is the correct one and stays.

3. **One-past-the-end reads are load-bearing.** `View::operator[]` asserts
   `index < size()`; two current asserts and the `open_output_file` /
   `std::fopen` path deliberately read the byte *after* the extent. Any
   `.data()[i]` → `[i]` conversion needs `i < size()` checked, not assumed.

4. **`.data()` destroys type information** (swarm T2 group B): it drops
   `alignas`, extent, and const-ness distinctions. In `core/align_simd.cpp` the
   `reinterpret_cast<CELL **>(s->qtable.data())` family is the
   container-to-array-base step the cast *needs*; there is no narrower way to
   spell it. In this codebase some `.data()` calls are the bug and some are the
   fix.

5. **Two of the touched files are the hot path.** `core/mask.cpp`'s `wo()` is
   **74%** of a search-command profile (masking runs on every query by
   default), and `core/align_simd.cpp`'s `aligncolumns_*` is ~95% of a
   `usearch_global` profile. Per `CLAUDE.md`, do not attribute a few-percent
   delta to a change here without a `callgrind` instruction count, and measure
   with `--qmask none --dbmask none` when the *target* is anything else.

6. **Run the test suite in the foreground.** swarm's T5 recorded this the hard
   way: launched as a background task the suite stopped after 45 of 1,071
   lines, reported **0 FAIL**, and logged no error at all. Check *how far a run
   got*, not just whether it reported failures. (See also the recorded
   `vsearch-tests` load flakes: `getseqs --labels` and `mergepairs
   --threads 1024` fail spuriously under concurrent build load.)


## Phasing

Branch from `dev` per `CLAUDE.md` (the tree is on `tmp_20260806093624`, so
`git checkout dev` first). One commit per phase; C4 gets its own branch.

| phase | content | sites | risk |
|---|---|---|---|
| 1 | C1 — `fatal()` `.c_str()` deletions | 13 | none |
| 2 | C2 — `set_deferred_error` / `warn` / `fastq_fatal` overloads | 9 | none |
| 3 | C7 — `random_shuffle` → `Span` | 1 | none |
| 4 | C6 — View/Span re-flattening, `decimal::to_text` | ~15 | low |
| 5 | C7 — `fasta_print`, `write_report`, `derep.cpp:1029`, `mergepairs` assigns | 6 | low |
| 6a | C5 — `searchinfo_s` per-thread arrays (index + reference) | 18 | medium — 4 files, threaded code |
| 6b | C5.2 — `si_m` → `Span<searchinfo_s>` + the missing assert | 3 | medium — see trap 1 |
| 7 | C5.1 — `search_exact` / `chimera` state structs hold the vector | 3 | medium |
| 8 | C3 — `Progress` owns its prompt | 4 + 58 | low, wide |
| 9 | C7 — `dust_core` → `Span<char>` | 1 | **callgrind required** |
| — | C4 — `search16` → `Span`, **own branch** | 21 + 7 | high |
| — | C9 — `largeread` / `largewrite` → `Span<Type>`, **own branch** | 16 | **high — binary format** |
| — | D (the rest), E, C8 | — | not doing; recorded above |

Three branches, decided 2026-08-09:

```
tmp_<ts>_references   phases 1-9   (~70 sites, all reversible)
tmp_<ts>_search16     C4           (SIMD entry point, four backends)
tmp_<ts>_udb          C9           (binary format + round-trip verification)
```

Phases 1-5 are independent of each other and of everything else. 6a before 6b
before 7. Phase 9 is separable and can be dropped if the measurement is not
clean. The two split-off branches keep a SIMD change and a binary-format change
from being merged on the strength of a review that was mostly about `.c_str()`
deletions.

**On C4's sequencing:** the `&parameter` plan noted that D1 is what would let
the three remaining `alignstats` call sites drop their locals the way
`chimera.cpp` did in its phase 4b. That is a reason to do C4 *before* the
`TBD_20260808_msa_cigar.md` work, not a reason to bundle them.


## Verification

Per `CLAUDE.md`, plus what the last two campaigns established as necessary:

- **Build** clean at every phase, zero new warnings, with the debug flags.
- **`cppcheck`** and **`run_clang_tidy.sh`** on each touched file, before and
  after — both report current state, not a diff. Expect
  `bugprone-easily-swappable-parameters` to *drop* in phase 6 and in C4; that
  is the point of those two.
- **Test suite in the foreground**: `(cd ~/Documents/src/vsearch-tests/ ;
  bash run_all_tests.sh ../vsearch/bin/vsearch | grep "FAIL")`, and check the
  line count of the full log, per trap 6.
- **Library**: `api_examples` `make test` against a **release** `libvsearch.a`
  (`make -C src libvsearch.a` is a separate explicit target). Required for
  phases 6-7 and C4 — the library session path goes through
  `cluster_query_init` and the chimera core. Expect 40 PASS / 0 FAIL.
- **Differential**, byte-identical after each phase, `--threads 1`, banner line
  stripped: `usearch_global` (alnout, userout with cigar+caln, blast6out,
  samout, uc, fastapairs), `cluster_fast` (msaout, consout, profile, centroids,
  uc), `allpairs_global`, `uchime_denovo`, `uchime_ref`, `search_exact`,
  `sintax` (**with `--randseed`** — it is nondeterministic without one),
  `derep_fulllength`, `sff_convert`, `fastx_syncpairs`, `otutabout`.
  Use `git worktree add --detach` for the reference build; `vsearch-tests`
  already holds `dev`.
- **Cross-compiles**: mingw, POWER, RISC-V/mips64el, each with a warning set
  identical to the pre-campaign baseline. Required for C4 (four SIMD backends)
  and phase 9.
- **Performance**: `callgrind` instruction counts on `usearch_global`,
  `cluster_fast`, `uchime_denovo` for phases 6, 9 and C4; nothing for phases
  1-5, which do not touch a loop. Expect neutrality — this is a typing and
  API-surface change, and `Span` is two words where a pointer was one. Where a
  count moves, disassemble both forms before believing a wall-clock number.


## What swarm learned (carried over)

From `~/Documents/src/swarm/TBD_20260806_misc.md:593-880`, the five items that
apply here directly:

1. **The split is the work.** A blanket sweep breaks code; swarm's group E was
   eight sites where `.data()` was the *fix*, each already carrying a comment
   saying so.
2. **Expect no measurable win.** It is a typing and API-surface change. Bound
   the mechanism first, and when a hot kernel is involved diff the disassembly
   rather than trusting hyperfine.
3. **`fatal()` was a deletion, not a conversion** — the same is true here, and
   it is 13 of vsearch's 38 `.c_str()` sites.
4. **Hoisting `.data()` into a local outside a hot loop is a design, not a
   smell** — measured at two prologue instructions and 16 bytes when inlined.
   vsearch's `core/unique.cpp` and `core/align_simd.cpp` are the same shape.
5. **Before turning a tautological assert into a real one, check whether the
   new condition is reachable from user input.** swarm's `db.cpp:755`
   (`assert(filename.c_str() != nullptr)`, which can never fire) became
   `!empty()` — and `swarm ""` then aborted with SIGABRT on ordinary bad input,
   violating the documented exit status. The assert was not the fix; the
   missing validation was. vsearch has the same shape available at
   `core/attributes.cpp:139` and `core/tax.cpp:84`.

Two structural differences worth noting:

- swarm's `fatal()` is a variadic template forwarding to `std::cerr <<`, and it
  has `utils/view_stream.hpp` so a `View<char>` can be an argument directly.
  vsearch's `fatal()` takes `char const *` or `std::string const &`, so a
  `View` argument still has to become a string. If C6 grows, the cheap
  equivalent is a `fatal(View<char>)` overload — **not** in this campaign.
- vsearch's `Span` → `View` conversion is `explicit`
  (`utils/span.hpp:108-112`, and its comment says "as in swarm"), so
  read-only-ness is marked at the call site. Expect visible
  `static_cast<View<char>>` at some of the phase-6 sites; that is the feature.


## Decisions (settled 2026-08-09, before any code was written)

1. **C3 — `Progress` owns its prompt.** `std::string prompt_`, constructor
   takes `std::string` by value. Deletes the undocumented lifetime
   requirement instead of documenting it. Phase 8.
2. **C7 — the clusters filename becomes a `std::string`**, with `.c_str()` at
   the `std::fopen` boundary as an honest group-B site. Not a `View<char>`
   overload. Phase 5.
3. **C8 — the cigar walkers are out of scope**, and belong to
   `TBD_20260808_msa_cigar.md`. Three plans converging on `alignstats` is how
   comments go stale.
4. **C9 — `largeread`/`largewrite` are converted**, against the census's own
   recommendation, because `Span`'s bounds assert is worth most exactly where
   the risk is highest. The extent is spelled at every call site. **Own
   branch.**
5. **C9's five "no container" sites need no second overload.** Three cases,
   all expressible: a `subspan` (`makeudb:211`), an existing view
   (`makeudb:263`), and a friend-access `Span` over the reserved buffer
   (`udb.cpp:408,455`). See the C9 table.
6. **`makeudb:249` writes the terminator instead of reading it**, and there is
   **no `header_view_z()` accessor**. Reading `len + 1` bytes out of
   `db.getheader(i)` couples the UDB writer to `Database`'s internal layout; an
   accessor would have made that coupling official, and writing the byte at the
   write site deletes it. The single byte goes through `largewrite` so it keeps
   the stream check and the `pos` accounting.
7. **`warn()` moves into `core/fastx.cpp`'s anonymous namespace** in the same
   commit as its overload, after `nm` confirms no other TU references it.
8. **Comments:** `core/otutable.cpp:164` gains a sentence explaining the
   deliberate one-past-the-end read. `core/showalign.cpp:385` and
   `core/getseq.cpp:200` stay **verbatim**.


## Open questions

None. Implementation can start at phase 1.

Two things to settle *during* the work rather than before it, because the
answer depends on what the code looks like by then:

- **C4's extent convention** — `make_span(v).first(count)` at every call site
  (one source of truth, but seven one-element spans in `cluster.cpp`'s hot
  loop) versus keeping the `sequences` parameter with an
  `assert(span.size() >= sequences)`. Decide with the callgrind numbers in
  hand, not from the plan.
- **Where C9's `static_assert(sizeof(unsigned int) == 4)` lives** — next to
  each of `largeread`/`largewrite`, or once in a shared UDB header. There is no
  shared UDB header today, so the answer is probably "twice, with the same
  message".
