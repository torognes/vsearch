# Extending the View/Span comparison to the remaining char comparisons (2026-07-27)

Follow-on to `TBD_20260725_replace_C_functions.md`. That document's `strcmp`
section ended with `View`/`Span` gaining a `char`-aware ordering
(`utils/element_order.hpp`, commit `d1169e95`) plus `compare()` and
`operator!=`. This one inventories the *other* places in the tree that hand-roll
a byte comparison the view operators now express directly.

Nothing here is a bug. These are readability and safety simplifications, with
one exception noted under "The safety point" below, where the current code is
correct only by virtue of a short-circuit.


## The pattern to look for

```cpp
if ((len_a == len_b) and std::equal(ptr_a, ptr_a + len_a, ptr_b))
```

A length test `and` a `std::equal` over a raw pointer and a length — which is
exactly what `View::operator==` does:

```cpp
auto operator==(View<Type> const & other) const -> bool {
  return size() == other.size()
    and std::equal(cbegin(), cend(), other.cbegin());
}
```

Four sites match. Two more are structural: a pair of parallel
pointer/length arrays that `View` exists to replace, and one file that spells
the same comparison two different ways.

Census taken with:

```sh
cd "${HOME}/Documents/src/vsearch/src/"
grep -rn "std::equal(\|std::memcmp" --include=*.cc --include=*.cpp --include=*.hpp . \
  | grep -v "^\./vendored/"
```


## 1. LCA taxonomy levels — `core/results.cpp:615-618` and `:650-653`

The best candidate: two near-identical copies, both inside
`results_show_lcaout`, in a file already converted to views.

```cpp
auto const * const cand_level = db.getheader(static_cast<uint64_t>(cand[k])) + cand_level_start[k][j];
auto const * const query_level = db.getheader(static_cast<uint64_t>(seqno)) + new_level_start[j];
if ((new_level_len[j] != cand_level_len[k][j]) or
    (not std::equal(cand_level, cand_level + new_level_len[j], query_level)))
  {
    match = false;
    break;
  }
```

becomes

```cpp
auto const cand_level = db.header_view(static_cast<uint64_t>(cand[k]))
                          .subspan(cand_level_start[k][j], cand_level_len[k][j]);
auto const query_level = db.header_view(static_cast<uint64_t>(seqno))
                           .subspan(new_level_start[j], new_level_len[j]);
if (cand_level != query_level)
  {
    match = false;
    break;
  }
```

Removes two pointer-arithmetic expressions, the manual length test, and the
unbounded `std::equal` overload (see below). `Database::header_view()` and
`View::subspan()` both already exist (`db.hpp:223`, `view.hpp:171`).

Note the offsets and lengths come from `tax_split()`, which fills
`int *` arrays (`tax.hpp:64`), so the `subspan` arguments widen from `int` to
`std::size_t`. They are offsets into a header and are non-negative by
construction, but the conversion should be explicit rather than implicit.

### The safety point

The current form uses the **three-iterator `std::equal`**, which bounds only
the first range: it reads `new_level_len[j]` bytes starting at `query_level`
without ever being told how long that range is. It is correct today *only*
because the `or` short-circuits when the two lengths differ, so `std::equal` is
never reached with mismatched lengths. That is a correctness argument the
reader has to reconstruct from operator precedence.

`View::operator==` compares the two sizes itself, and `subspan` asserts
`offset <= size()` and `count <= size() - offset` in a debug build, so the
bound becomes structural instead of incidental. This is the one place in this
document where the rewrite buys more than readability.


## 2. Sintax candidate levels — `commands/sintax.cpp:158-159` and `:209-212`

Same comparison, but here the data structure is the real finding:

```cpp
std::array<std::array<char const *, tax_levels>, bootstrap_count> cand_level_name_start {{}};
std::array<std::array<int,          tax_levels>, bootstrap_count> cand_level_name_len   {{}};
```

Two parallel arrays holding the pointer/length pair `View<char>` bundles. The
comparison at `:209` is then the familiar shape:

```cpp
if ((cand_level_name_len[cand_i][level] == cand_level_name_len[cand_j][level]) &&
    std::equal(cand_level_name_start[cand_i][level],
               cand_level_name_start[cand_i][level] + cand_level_name_len[cand_i][level],
               cand_level_name_start[cand_j][level]))
```

Collapsing the two arrays into one
`std::array<std::array<View<char>, tax_levels>, bootstrap_count>` drops the
second array, drops the pointer arithmetic at the fill site (`:181`), and
reduces the test to:

```cpp
if (cand_level_name[cand_i][level] == cand_level_name[cand_j][level])
```

Larger than a one-liner — the fill site and every index of both arrays move —
but it is the same thread as item 1, and the two files solve the same problem
(`results.cpp` keeps offsets and adds them at compare time; `sintax.cpp`
converts to pointers up front). Worth doing them together so the tree ends up
with one idiom.


## 3. FASTQ plus-line check — `core/fastq.cpp:477-484`

```cpp
if (input_handle->header_buffer.length == input_handle->plusline_buffer.length)
  {
    if (not std::equal(input_handle->header_buffer.data(),
               input_handle->header_buffer.data() + input_handle->header_buffer.length,
                input_handle->plusline_buffer.data()))
      {
        plusline_invalid = true;
      }
  }
else
  {
    /* a bare "+" or "+\r" is allowed */
  }
```

The inner test is `header_view() != plusline_view()`. `fastx_s::header_view()`
already exists (`fastx.hpp:254`); `plusline_buffer` would need the matching
accessor, which is two lines beside it.

A partial win only: the `else` branch has real logic (a plus line of length
0, or 2 when the second byte is `'\r'`, is legal), so the `if`/`else`
structure has to stay. Do it when touching `fastq.cpp` for another reason
rather than on its own.


## 4. Userfield name lookup — `utils/userfields.cpp:149` — marginal

```cpp
if ((std::strlen(*valid_userfield) == field_length) and std::equal(ptr, ptr + field_length, *valid_userfield))
```

Textbook `View::operator==`, but `valid_userfield` walks a table of
`char const *` literals, so the `std::strlen` survives the rewrite and the
result is no shorter:

```cpp
if (View<char>{ptr, field_length} == View<char>{*valid_userfield, std::strlen(*valid_userfield)})
```

Only worth doing as part of turning the table itself into `View<char>` (or
`std::string`) entries, which would also retire the `strlen` — a separate
decision about `userfields.cpp`, not a comparison cleanup.


## 5. Not a view candidate: the magic-byte tests — `core/fastx.cpp`

The same gzip/bzip2 sniff is spelled two ways in one file:

```cpp
/* :386, :390 */  if (std::equal(magic.begin(), magic.end(), magic_gzip.begin()))
/* :488, :493 */  if (std::memcmp(first, magic_gzip.data(), 2) == 0)
```

Worth unifying — it also retires the last two `std::memcmp` outside the SIMD
code — but **on `std::equal`, not on `View`**. `magic_gzip` and `magic_bzip`
are `constexpr std::array<unsigned char, 2>` (`fastx.cpp:95-96`) and already
carry `begin()`/`end()`; a view over them would add a cast from
`unsigned char` to `char` and buy nothing. The fix is to use the `:386`
spelling at `:488` as well, with the literal `2` gone from both.


## Explicitly not candidates

These compare with a custom predicate. `View::operator==` compares bytes
exactly, so it is the wrong tool and they stay as they are:

| Location | Predicate |
|---|---|
| `commands/cut.cpp:145` | `is_equivalent_4bit_rhs` — nucleotide equivalence, so `A` matches `N` |
| `utils/compare_strings_nocase.cpp:95`, `:105`, `:117` | `to_upper(lhs) == to_upper(rhs)` (`:75`), through `utils/ascii_case.hpp` |
| `utils/seqcmp.cpp:77-78` | `map_4bit` — nucleotide ordering, plus an early stop at `'\0'` |

The three `std::strcmp` that survive (`cli.cc:4448`, `core/fastx.cpp:342`,
`utils/open_file.cpp:89`) are also not candidates: they compare an `argv`
pointer against the one-character literal `"-"`, where a view buys nothing and
a `std::string` temporary would allocate. See the `strcmp` section of
`TBD_20260725_replace_C_functions.md`.


## Suggested ordering

1. **`core/results.cpp` ×2** — smallest diff, the only one with a safety
   argument, and the file is already view-based.
2. **`commands/sintax.cpp`** — the parallel-array collapse, done in the same
   pass so both taxonomy comparisons end up with one idiom.
3. **`core/fastx.cpp` magic bytes** — a consistency fix that happens to
   retire the last non-SIMD `memcmp`. Independent of the rest.
4. **`core/fastq.cpp` plus line** — bundle with other work in that file.
5. **`utils/userfields.cpp`** — only alongside a decision about the
   `valid_userfield` table.


## Verification

All five are output-affecting: items 1 and 2 decide `--lcaout` and `--sintax`
taxonomy strings, item 3 decides which fatal message a compressed file
produces, item 4 decides a FASTQ parse error. The bar is the same one used for
the `strlen` and `strcmp` passes — see the recipe in
`TBD_20260725_replace_C_functions.md`:

- a `dev` reference binary from a `git worktree add --detach` (plain
  `git worktree add … dev` fails: `~/Documents/src/vsearch-tests` already holds
  `dev`), run with `--threads 1` so output order is deterministic, diffed
  byte-for-byte after normalising the command-line first line;
- datasets with **non-ASCII bytes in headers**, since taxonomy fields are cut
  out of headers and that is where the `char` signedness divergence lives;
- for items 1 and 2 specifically, a taxonomy dataset where two candidates share
  a prefix but differ in level length, so the length half of `operator==` is
  exercised and not just the byte half.

Both have existing suite coverage to fall back on: `scripts/sintax.sh` is part
of the default `run_all_tests.sh` run (unlike `scripts/orient.sh`, which is
commented out at line 69), and `--lcaout` is exercised by
`scripts/usearch_global.sh`, `search_exact.sh`, `allpairs_global.sh` and
`fixed_bugs.sh`. That coverage is ASCII-only, though, which is precisely the
half that cannot catch a signedness or bounds mistake — so the purpose-built
dataset above is the real check, not the suite.
