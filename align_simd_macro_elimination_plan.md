# Plan: eliminate the remaining function-like macros in `align_simd.cc`

## Goal

Replace every remaining function-like macro in `src/align_simd.cc` with
type-safe (`inline`/`constexpr`) functions, following the C++11 refactoring
already applied to the surrounding `v_*` intrinsic wrappers and to the
equivalent code in swarm.

This is a **plan only**; no source change to `align_simd.cc` is part of this
commit.

## Reference: how swarm solved the same problem

swarm went through this exact conversion (`src/search16.cc`,
`src/arch/*/intrinsics_to_functions.{h,cc}`):

- `v_load` / `v_store` macros became `v_load16()` / `v_store16()` functions,
  with a single overloaded `cast_vector16()` helper holding the pointer cast
  that the macro used to hide.
- the `ALIGNCORE` macro became one `onestep_16()` `inline` function taking the
  DP cells **by reference**, with the per-cell temporaries as locals.
- the PPC-only byte-permute (`VECTORBYTEPERMUTE` + `vec_cmpgt` + `perm`) was
  folded **inside** a PPC `v_mask_*16()` function that returns the same
  `uint16_t` bitmask shape as x86_64/aarch64, so the single `onestep_16`
  body serves all three architectures.

We mirror that here.

## Inventory of remaining macros

| Macro | Defs | Location | Notes |
|-------|------|----------|-------|
| `v_init(a..h)` | 3 (PPC / aarch64 / x86) | 157, 210, 290 | aarch64/PPC use a compound literal that already triggers `-Wpedantic` (see comment at line 209) |
| `v_load(a)` | 3 | 158, 211, 291 | wraps a C-style/`reinterpret_cast` pointer cast |
| `v_store(a,b)` | 3 | 159, 212, 292 | same |
| `VECTORBYTEPERMUTE` | 2 (GNU / IBM) | 691, 693 | PPC only |
| `ALIGNCORE(...)` | 2 (PPC packed-`RES` / x86+aarch64 `PATH`) | 702, 737 | the innermost DP kernel; relies on enclosing-scope `E`, `HE`, `HF` (and on PPC `RES*`) |

Call sites: `v_load`/`v_store` in `dprofile_fill16` (≈583-663) and the two
`aligncolumns_*` kernels; `ALIGNCORE` in `aligncolumns_first` (757-922) and
`aligncolumns_rest` (925-…); `v_init` once at line 1619.

## Why convert

- type safety + no token-pasting / double-evaluation pitfalls;
- the macros currently leak and mutate enclosing-scope names (`E`, `HE`, `HF`,
  `RES`, `RES1`, `RES2`), which is fragile and forces the callers to declare
  those temporaries;
- removes C-style casts hidden in the load/store/permute macros
  (CLAUDE.md: prefer `static_cast`, avoid C-style casts);
- converting `v_init` to a function removes the `-Wpedantic` compound-literal
  warning documented at line 209;
- consistency with swarm and with the already-converted `v_add`, `v_max`, … .

## Step-by-step conversion

### 1. `v_load` / `v_store`  → functions + `cast_vector16` helper
- Add one overloaded helper per arch that performs the pointer cast in a single
  place, e.g. `cast_vector16(short const *) -> VECTOR_SHORT const *` and the
  non-const variant (PPC keeps its `vec_ld`/`vec_st` wrapping).
- Define `inline auto v_load(VECTOR_SHORT const * ptr) -> VECTOR_SHORT` and
  `inline auto v_store(VECTOR_SHORT * ptr, VECTOR_SHORT v) -> void`.
- Update call sites to wrap their pointer argument with `cast_vector16(...)`,
  and replace `base + d[k] + i` pointer arithmetic with `std::next(...)`
  (CLAUDE.md: avoid pointer arithmetic). Keep the existing index expressions
  otherwise unchanged.

### 2. `v_init` → function
- `inline auto v_init(short a, … short h) -> VECTOR_SHORT` per arch
  (x86 `_mm_set_epi16(h,…,a)`, aarch64/PPC build the vector without a compound
  literal). One call site to update.
- **Design question for review:** eight `short` parameters are maximally
  "easily swappable" (CLAUDE.md). Options: (a) keep the 8 scalars to match the
  intrinsic 1:1, or (b) take `std::array<short, CHANNELS>`. The single call
  uses `{-1,0,0,0,0,0,0,0}`, so (b) reads cleanly. Recommend (b); confirm.

### 3. `VECTORBYTEPERMUTE` + PPC `ALIGNCORE`  (the only real decision)
Two options:

- **Option A — unify (recommended, mirrors swarm).** Add a PPC
  `inline auto v_mask_gt(VECTOR_SHORT, VECTOR_SHORT) -> unsigned short` that
  folds in `vec_cmpgt` + the byte-permute (currently `VECTORBYTEPERMUTE` +
  `perm`) and returns the same bitmask layout as the existing x86_64/aarch64
  `v_mask_gt`. `ALIGNCORE` then collapses to **one** function for all arches
  (step 4). This deletes `VECTORBYTEPERMUTE`, the packed-`RES` machinery
  (`RES`/`RES1`/`RES2`, `vec_perm`, `perm_merge_long_low/high`) and the
  `#ifdef __PPC__` blocks inside both `aligncolumns_*`.
  - **Risk:** it changes how the PPC `dir` buffer is written (4 words per step
    instead of one packed vector). `backtrack16` already consumes the
    4-words-per-step layout on x86_64/aarch64, so unifying makes PPC use the
    same, already-tested path — but this cannot be **run** locally (no PPC
    hardware; cross-compile only). Must be flagged for human review and, ideally,
    validated on real PPC/CI before merge.

- **Option B — conservative.** Translate each `ALIGNCORE` variant into its own
  `inline` function, preserving the PPC packed-`RES` output byte-for-byte
  (keep the permute as a small `inline` function instead of a macro). Lower
  risk, keeps the `#ifdef` in the kernels, less cleanup.

Recommend **A** for maximal simplification and parity with swarm, but defer to
the maintainer because the PPC path is not locally runnable.

### 4. `ALIGNCORE` → `onestep` inline function
- Signature (Option A, single function):
  `inline auto onestep(VECTOR_SHORT & H, VECTOR_SHORT & N, VECTOR_SHORT & F,
   VECTOR_SHORT V, unsigned short * path, VECTOR_SHORT & E,
   VECTOR_SHORT QR_q, VECTOR_SHORT R_q, VECTOR_SHORT QR_t, VECTOR_SHORT R_t,
   VECTOR_SHORT & h_min, VECTOR_SHORT & h_max) -> void`.
- `HE`/`HF` become **locals** inside `onestep`; drop their declarations from
  both `aligncolumns_*` kernels. `E` is taken by reference (it is updated).
- Body is the current x86/aarch64 `ALIGNCORE` (four `*path` writes via
  `v_mask_gt`); on PPC the same body works once step 3A provides PPC
  `v_mask_gt`.
- **Design question for review:** this carries many same-typed (`VECTOR_SHORT`)
  parameters — swappable-parameter smell. Optionally group the four
  penalty vectors into a small `struct gap_penalties { VECTOR_SHORT QR_q, R_q,
  QR_t, R_t; }` to cut the count and improve readability. Recommend grouping;
  confirm.

### 5. Cleanup
- Remove the now-unused `perm_merge_long_low/high`, `RES*` declarations,
  `VECTORBYTEPERMUTE`, and the `#ifdef __PPC__` direction-packing blocks
  (Option A only).

## Verification

- **x86_64 debug** build (`-O0 -ggdb3`) and **release** build: zero new
  warnings; confirm the `-Wpedantic` compound-literal warning is gone.
- `cppcheck --force --enable=warning,style --language=c++ --std=c++11 align_simd.cc`.
- **Cross-compile** the non-x86 paths to compile-check them: PPC
  (`--host=powerpc64le-linux-gnu`) and the MIPS/RISC-V host from CLAUDE.md;
  confirm whether an aarch64 cross-toolchain is available (else NEON path is
  compile-checked only via simde).
- **Test suite:** `(cd ~/Documents/src/vsearch-tests/ ; bash run_all_tests.sh
  ../vsearch/bin/vsearch | grep FAIL)`. Pay attention to alignment-exercising
  commands that hit `backtrack16` (`--usearch_global … --alnout/--userfields`
  with cigar, `--allpairs_global`, `--cluster_*`).
- **Performance:** `ALIGNCORE` is the innermost DP loop. Run the
  `usearch_global` hyperfine test before/after, and `make -j PROFILE=1`, to
  confirm the `inline` functions are inlined and there is no regression.

## Suggested commit breakdown (small commits, per CLAUDE.md)

1. `v_load`/`v_store` + `cast_vector16` helper (all arches, call sites).
2. `v_init` function (all arches) — also removes the `-Wpedantic` warning.
3. *(Option A)* PPC `v_mask_gt` folding `VECTORBYTEPERMUTE`.
4. `ALIGNCORE` → `onestep` inline function; drop `HE`/`HF` decls and the
   per-arch `#ifdef` in both kernels.
5. Remove dead `perm_merge_long_*` / `RES*` / `VECTORBYTEPERMUTE`.

## Open questions for the maintainer

1. **PPC `ALIGNCORE`: Option A (unify, swarm-style) or B (conservative)?**
   A is cleaner but the PPC `dir` layout change is not runnable locally.
2. `v_init`: eight `short` args vs `std::array<short, CHANNELS>`?
3. `onestep`: keep the four gap-penalty vectors flat, or group them in a
   `struct` to reduce swappable parameters?
4. Is real aarch64/PPC runtime testing available in CI, or do we rely on
   cross-compilation plus the x86_64 test suite?
