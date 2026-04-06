# vsearch Library API

vsearch can be embedded as a static library (`libvsearch.a`) in other
applications. This document provides comprehensive documentation of the
library API: build integration, initialization protocol, available
operations, result structures, error handling, memory ownership, thread
safety, and platform support.

## Table of contents

- [Building](#building)
- [Architecture overview](#architecture-overview)
- [Initialization](#initialization)
- [Database management](#database-management)
- [Chimera detection](#chimera-detection)
- [Global search](#global-search)
- [Clustering](#clustering)
- [Dereplication](#dereplication)
- [Paired-end merging](#paired-end-merging)
- [Sequence masking](#sequence-masking)
- [Configuration reference](#configuration-reference)
- [Error handling](#error-handling)
- [Memory management](#memory-management)
- [Thread safety](#thread-safety)
- [Platform support](#platform-support)
- [Examples](#examples)

---

## Building

### Static library

```bash
# From a fresh clone (autotools not pre-generated):
autoreconf -fi
./configure
make -C src libvsearch.a

# From the release branch (autotools pre-generated):
./configure
make -C src libvsearch.a
```

This produces `src/libvsearch.a`, a position-independent static archive
suitable for linking into shared libraries (e.g., loadable extensions).
All objects are compiled with `-fPIC`. The library excludes `main()` via
the `VSEARCH_NO_MAIN` preprocessor guard.

### CMake integration (ExternalProject)

For projects that build vsearch as a dependency via CMake:

```cmake
include(ExternalProject)
ExternalProject_Add(vsearch_build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/path/to/vsearch
    CONFIGURE_COMMAND ./configure --disable-pdfman --disable-zlib --disable-bzip2
    BUILD_COMMAND make -C src libvsearch.a
    BUILD_IN_SOURCE TRUE
    INSTALL_COMMAND ""
)
```

On macOS cross-compilation (arm64 runner targeting x86_64), pass
`-arch x86_64` via `CXXFLAGS` and `CFLAGS`, and set
`--host=x86_64-apple-darwin` on the configure command.

Pre-generated autotools files are committed to the repository so
`autoreconf` is not required at build time. If timestamps cause
`make` to attempt regeneration, touch the generated files before
building.

### Link dependencies

Required:

- `-lpthread` (threading primitives)
- `-ldl` (dynamic library loading)

Optional (detected at configure time, can be disabled):

- `-lz` (zlib — `--disable-zlib` to omit)
- `-lbz2` (bzip2 — `--disable-bzip2` to omit)

### Include path

Add `src/` to your include path. The single entry point header is:

```cpp
#include "vsearch_api.h"
```

This transitively includes `vsearch.h` (all global declarations) and
the module headers for each API subsystem.

### C++ standard

The library requires C++11 or later. Build with `-std=c++11` (or
higher).

---

## Architecture overview

### Global state model

vsearch was designed as a CLI tool and uses approximately 200 global
`opt_*` variables to control all behavior. The library API inherits
this architecture. There is no context object — configuration is set
by writing to globals, and internal state (database, k-mer index) is
also global.

**Consequences:**

- Only one session can be active per process at a time
- Initialization and configuration are not thread-safe
- Calling `vsearch_init_defaults()` resets ALL options

A session mutex serializes access: `vsearch_init_defaults()` acquires
it, `vsearch_session_end()` releases it. Concurrent callers block until
the mutex is released.

### Per-thread working state

Operations that support multi-threaded execution (chimera detection,
search, clustering) provide opaque per-thread state objects. Each
thread allocates its own instance; the objects must not be shared
across threads. The global database and k-mer index are read-only
after initialization, so no locking is needed during computation.

### Output and I/O

The library defaults to silent operation: `opt_quiet` and
`opt_no_progress` are both `true` after `vsearch_init_defaults()`.
No output is written to stdout or stderr by library functions under
default settings. The caller has full control over I/O.

---

## Initialization

### Required sequence

Every session must follow this protocol:

```cpp
#include "vsearch_api.h"

// 1. Set all ~200 globals to library defaults (quiet, no progress).
//    Acquires the session mutex — blocks if another session is active.
vsearch_init_defaults();

// 2. Override globals for your use case.
opt_wordlength = 8;
opt_id = 0.97;

// 3. Resolve sentinel values to computed defaults.
//    Must be called AFTER setting overrides.
vsearch_apply_defaults_fixups();

// 4-N. Use the library (load DB, run operations, etc.)
// ...

// Final. Release the session mutex.
vsearch_session_end();
```

**Step 3 is mandatory.** `vsearch_apply_defaults_fixups()` resolves
sentinel values that depend on other options:

- `opt_maxhits`: 0 → `INT64_MAX`
- `opt_minwordmatches`: -1 → wordlength-based default from
  `searchcore.h` lookup table
- Gap penalties: raw values are adjusted by subtracting extension
  penalties (e.g., `opt_gap_open_query_interior` 20 → 18 after
  subtracting `opt_gap_extension_query_interior` 2)

Skipping this step causes silent failures — candidates are rejected
due to wrong gap penalties, or no k-mer matches are found.

### Re-initialization

Multiple sequential sessions in the same process are supported.
Repeat the full sequence (init → configure → fixups → work → cleanup
→ session end) for each session. `vsearch_init_defaults()` resets all
globals to fresh defaults, and `vsearch_apply_defaults_fixups()`
correctly re-applies gap penalty adjustments each time.

See `examples/example_reinit.cc` for a tested multi-session example.

### Initialization functions

| Function | Description |
|----------|-------------|
| `vsearch_init_defaults()` | Set all opt_* globals to library defaults. Acquires session mutex. |
| `vsearch_apply_defaults_fixups()` | Resolve sentinel values. Call after setting overrides. |
| `vsearch_session_end()` | Release session mutex. Call after all cleanup. |

---

## Database management

vsearch uses a global in-memory sequence database. Sequences are loaded
via `db_init()` + `db_add()`, then indexed for k-mer search.

### Loading sequences

```cpp
// Reset database state (also frees any previous data)
db_init();

// Add sequences one at a time
for (int i = 0; i < n_seqs; i++) {
    db_add(false,             // is_fastq: false for FASTA
           headers[i],        // null-terminated header string
           sequences[i],      // null-terminated DNA sequence (ACGT, uppercase)
           nullptr,           // quality string (nullptr for FASTA)
           strlen(headers[i]),
           strlen(sequences[i]),
           abundance);        // size annotation (1 for uniform abundance)
}
```

**Requirements:**

- `db_init()` must be called before `db_add()`. Without it, internal
  statistics (shortest sequence length) are corrupted.
- Sequences must be DNA (ACGT). Convert U to T before calling.
- `db_init()` calls `db_free()` internally, so it is safe to call on
  an already-populated database.

### Reading from files

For file-based workflows, `db_read()` loads sequences from FASTA/FASTQ
files directly:

```cpp
db_read("sequences.fasta", 0);  // 0 = do not upcase
```

### Masking and indexing

After loading sequences, apply masking and build the k-mer index:

```cpp
dust_all();                              // DUST low-complexity masking
dbindex_prepare(1, opt_dbmask);          // allocate index (1 = use bitmap)
dbindex_addallsequences(opt_dbmask);     // index all sequences
```

### Sorting

Some operations require the database to be pre-sorted:

```cpp
db_sortbylength();               // for cluster_fast
db_sortbylength_shortest_first();
db_sortbyabundance();            // for cluster_size
```

### Cleanup

```cpp
dbindex_free();
db_free();
```

### Database functions

| Function | Description |
|----------|-------------|
| `db_init()` | Reset database state. Call before `db_add()`. |
| `db_add(is_fastq, header, seq, qual, hlen, slen, abund)` | Add one sequence. |
| `db_read(filename, upcase)` | Load sequences from file. |
| `db_free()` | Free all database memory. |
| `db_getsequencecount()` | Number of loaded sequences. |
| `db_getnucleotidecount()` | Total nucleotides across all sequences. |
| `db_getlongestsequence()` | Length of longest sequence. |
| `db_getshortestsequence()` | Length of shortest sequence. |
| `db_getlongestheader()` | Length of longest header. |
| `db_getheader(seqno)` | Pointer to header string (database-owned). |
| `db_getsequence(seqno)` | Pointer to sequence string (database-owned). |
| `db_getsequencelen(seqno)` | Sequence length. |
| `db_getheaderlen(seqno)` | Header length. |
| `db_getabundance(seqno)` | Abundance annotation. |
| `db_getquality(seqno)` | Quality string (FASTQ only). |
| `db_is_fastq()` | Whether database contains FASTQ data. |
| `db_sortbylength()` | Sort by length, longest first. |
| `db_sortbylength_shortest_first()` | Sort by length, shortest first. |
| `db_sortbyabundance()` | Sort by abundance, highest first. |

### K-mer index functions

| Function | Description |
|----------|-------------|
| `dbindex_prepare(bitmap, mask)` | Allocate k-mer index. `bitmap=1` enables bitmap mode (required for clustering). |
| `dbindex_addallsequences(mask)` | Index all loaded sequences. |
| `dbindex_addsequence(seqno, mask)` | Index a single sequence (for incremental indexing). |
| `dbindex_free()` | Free k-mer index memory. |

---

## Chimera detection

Reference-based chimera detection (UCHIME algorithm: Edgar et al. 2011,
Bioinformatics 27:2194-2200) with vsearch's SIMD-optimized
Needleman-Wunsch alignment.

### Lifecycle

**Single-threaded** (convenience wrappers):

```cpp
struct chimera_info_s * ci = chimera_info_alloc();
chimera_detect_init(ci);       // session + per-thread init

struct chimera_result_s result;
chimera_detect_single(ci, query_seq, query_head,
                      query_len, query_abundance, &result);

chimera_detect_cleanup(ci);    // per-thread + session cleanup
chimera_info_free(ci);
```

**Multi-threaded** (split API):

```cpp
// Session init: once, single-threaded, after DB indexed
chimera_session_init();

// Per-thread init: one handle per thread
struct chimera_info_s * ci1 = chimera_info_alloc();
chimera_detect_thread_init(ci1);
struct chimera_info_s * ci2 = chimera_info_alloc();
chimera_detect_thread_init(ci2);

// Detection: thread-safe with separate handles
// thread 1: chimera_detect_single(ci1, ...)
// thread 2: chimera_detect_single(ci2, ...)

// Per-thread cleanup
chimera_detect_thread_cleanup(ci1);
chimera_info_free(ci1);
chimera_detect_thread_cleanup(ci2);
chimera_info_free(ci2);

// Session cleanup: once, single-threaded
chimera_session_cleanup();
```

### Result structure

`chimera_result_s` matches vsearch's 18-column `--uchimeout` format:

| Field | Type | Description |
|-------|------|-------------|
| `score` | `double` | h-score |
| `query_label` | `char[1024]` | Query header (truncated to 1023 chars) |
| `parent_a_label` | `char[1024]` | Parent A header |
| `parent_b_label` | `char[1024]` | Parent B header |
| `closest_parent_label` | `char[1024]` | Closest parent header |
| `id_query_model` | `double` | Query-to-model identity % |
| `id_query_a` | `double` | Query-to-parent-A identity % |
| `id_query_b` | `double` | Query-to-parent-B identity % |
| `id_a_b` | `double` | Parent-A-to-parent-B identity % |
| `id_query_top` | `double` | Query-to-closest-parent identity % |
| `left_yes` | `int` | Left segment YES votes |
| `left_no` | `int` | Left segment NO votes |
| `left_abstain` | `int` | Left segment ABSTAIN votes |
| `right_yes` | `int` | Right segment YES votes |
| `right_no` | `int` | Right segment NO votes |
| `right_abstain` | `int` | Right segment ABSTAIN votes |
| `divergence` | `double` | Divergence between parents |
| `flag` | `char` | Classification: `'Y'`, `'N'`, or `'?'` |

For non-chimeric results (`flag='N'`), only `query_label` and `flag`
are populated; all other fields are zero/empty.

### De novo mode

For de novo chimera detection (`--uchime_denovo` equivalent), set
`opt_chimeras_denovo` to a non-null value before calling
`chimera_detect_init()`. In this mode:

- Process queries in decreasing abundance order
- After classifying a query as non-chimeric, add it to the reference
  database with `db_add()` and index it with `dbindex_addsequence()`
- The `query_abundance` parameter to `chimera_detect_single()` controls
  abundance skew filtering

De novo mode is inherently sequential (single-threaded).

### Key options

| Global | Type | Default | Description |
|--------|------|---------|-------------|
| `opt_wordlength` | `int64_t` | 8 | K-mer length. Set before `vsearch_apply_defaults_fixups()`. |
| `opt_minh` | `double` | 0.28 | Minimum h-score for chimera classification. |
| `opt_xn` | `double` | 8.0 | Weight of no votes. |
| `opt_dn` | `double` | 1.4 | Pseudo-count prior for votes. |
| `opt_mindiv` | `double` | 0.8 | Minimum divergence. |
| `opt_mindiffs` | `int64_t` | 3 | Minimum differences from parents. |
| `opt_abskew` | `double` | 2.0 | Abundance skew (de novo mode). |

### Functions

| Function | Description |
|----------|-------------|
| `chimera_info_alloc()` | Allocate opaque per-thread state. Returns heap pointer. |
| `chimera_info_free(ci)` | Free per-thread state. Does NOT call cleanup. Null-safe. |
| `chimera_session_init()` | Session-level init: set search globals, init mutexes. Call once after DB indexed. |
| `chimera_session_cleanup()` | Session-level teardown: destroy mutexes. Call after all per-thread cleanup. |
| `chimera_detect_thread_init(ci)` | Per-thread init: allocate SIMD aligners, k-mer finders. |
| `chimera_detect_thread_cleanup(ci)` | Per-thread teardown: free resources. |
| `chimera_detect_single(ci, seq, head, len, abund, result)` | Detect chimera for one query. Returns 0 on success. |
| `chimera_detect_init(ci)` | Convenience: `session_init` + `thread_init`. Single-threaded only. |
| `chimera_detect_cleanup(ci)` | Convenience: `thread_cleanup` + `session_cleanup`. Single-threaded only. |

---

## Global search

Search query sequences against the loaded database, reporting hits
above a minimum identity threshold.

### Lifecycle

```cpp
// Configure search parameters
opt_id = 0.97;           // minimum 97% identity
opt_maxaccepts = 3;      // return up to 3 hits
opt_maxrejects = 16;     // give up after 16 rejects

// Allocate and initialize per-thread state
struct searchinfo_s * si = search_info_alloc();
search_init(si);         // call after DB indexed

// Search queries
struct search_result_s results[10];
int result_count;
search_single(si, query_seq, query_head, query_len,
              query_abundance, results, 10, &result_count);

// Cleanup
search_cleanup(si);
search_info_free(si);
```

### Result structure

`search_result_s` contains per-hit alignment details:

| Field | Type | Description |
|-------|------|-------------|
| `target` | `int` | Database sequence index. |
| `target_label` | `char[1024]` | Target header (truncated to 1023 chars). |
| `id` | `double` | Percent identity (per `opt_iddef`). |
| `matches` | `int` | Matching columns. |
| `mismatches` | `int` | Mismatching columns. |
| `gaps` | `int` | Gap columns. |
| `alignment_length` | `int` | Total alignment length. |
| `query_length` | `int` | Query sequence length. |
| `target_length` | `int` | Target sequence length. |
| `accepted` | `bool` | Whether hit passed the identity threshold. |

Results are ordered by identity (descending). The caller provides the
result array and its capacity; `result_count` is set to the number of
hits returned (up to `max_results` and `opt_maxaccepts`).

### Key options

| Global | Type | Default | Description |
|--------|------|---------|-------------|
| `opt_id` | `double` | 0.0 | Minimum identity threshold (0.0–1.0). |
| `opt_maxaccepts` | `int64_t` | 1 | Stop after N accepted hits. |
| `opt_maxrejects` | `int64_t` | 32 | Stop after N rejected candidates. |
| `opt_iddef` | `int64_t` | 2 | Identity definition (0–4). |
| `opt_wordlength` | `int64_t` | 8 | K-mer length for candidate selection. |

### Functions

| Function | Description |
|----------|-------------|
| `search_info_alloc()` | Allocate opaque per-thread state. Returns heap pointer. |
| `search_info_free(si)` | Free per-thread state. Null-safe. |
| `search_init(si)` | Initialize SIMD aligners and k-mer finders. Call after DB indexed. |
| `search_single(si, seq, head, len, size, results, max, count)` | Search one query. Thread-safe with per-thread si. |
| `search_cleanup(si)` | Free search resources. Call before `search_info_free`. |

---

## Clustering

Greedy centroid-based clustering (equivalent to `--cluster_fast` or
`--cluster_size`). Sequences are processed sequentially; each is either
assigned to an existing cluster or becomes a new centroid.

### Lifecycle

```cpp
// Sort database (required)
db_sortbylength();          // for cluster_fast behavior
// or: db_sortbyabundance() for cluster_size behavior

// Prepare index — do NOT call dbindex_addallsequences().
// Centroids are indexed incrementally during clustering.
dbindex_prepare(1, opt_qmask);

// Allocate and initialize session
struct cluster_session_s * cs = cluster_session_alloc();
cluster_session_init(cs);

// Assign sequences sequentially (0, 1, 2, ...)
struct cluster_result_s result;
for (int i = 0; i < db_getsequencecount(); i++) {
    cluster_assign_single(cs, i, &result);
    // result.is_centroid, result.cluster_id, result.identity, ...
}

// Cleanup
cluster_session_cleanup(cs);
cluster_session_free(cs);
```

### Result structure

| Field | Type | Description |
|-------|------|-------------|
| `is_centroid` | `bool` | `true` if sequence starts a new cluster. |
| `cluster_id` | `int` | Cluster number (0-based). |
| `centroid_seqno` | `int` | Database seqno of the cluster centroid. |
| `centroid_label` | `char[1024]` | Centroid header (truncated to 1023 chars). |
| `identity` | `double` | Identity to centroid (100.0 if centroid). |
| `cigar` | `char[4096]` | CIGAR alignment string (empty if centroid). |
| `cigar_truncated` | `bool` | `true` if CIGAR was truncated to fit buffer. |

### Key options

| Global | Type | Default | Description |
|--------|------|---------|-------------|
| `opt_id` | `double` | 0.0 | Identity threshold for cluster membership. |
| `opt_strand` | `int64_t` | 1 | 1 = plus strand only, 2 = both strands. |
| `opt_maxaccepts` | `int64_t` | 1 | Accepts before stopping search. |
| `opt_maxrejects` | `int64_t` | 32 | Rejects before stopping search. |

### Functions

| Function | Description |
|----------|-------------|
| `cluster_session_alloc()` | Allocate opaque session state. |
| `cluster_session_free(cs)` | Free session state. Null-safe. |
| `cluster_session_init(cs)` | Initialize session. DB must be sorted; `dbindex_prepare()` called but NOT `dbindex_addallsequences`. |
| `cluster_assign_single(cs, seqno, result)` | Assign one sequence. Must be called in seqno order (0, 1, 2, ...). |
| `cluster_session_cleanup(cs)` | Free session resources. |

---

## Dereplication

Full-length dereplication: collapse identical sequences and sum their
abundances. This operation does NOT use the global database — it
maintains its own internal state.

### Lifecycle

```cpp
struct derep_session_s * ds = derep_session_alloc();
derep_session_init(ds);

// Add sequences (normalized internally: uppercase, U→T)
for (int i = 0; i < n_seqs; i++) {
    derep_add_sequence(ds, headers[i], sequences[i],
                       strlen(sequences[i]), 1);
}

// Retrieve results (sorted by abundance, descending)
struct derep_result_s results[1000];
int result_count;
derep_get_results(ds, results, 1000, &result_count);

derep_session_cleanup(ds);
derep_session_free(ds);
```

### Result structure

| Field | Type | Description |
|-------|------|-------------|
| `header` | `const char *` | Representative header. Session-owned; valid until cleanup. |
| `sequence` | `const char *` | Normalized sequence. Session-owned; valid until cleanup. |
| `abundance` | `uint64_t` | Total abundance (sum of identical sequences). |
| `seqlen` | `uint64_t` | Sequence length. |
| `count` | `int` | Number of input sequences collapsed. |

**Pointer lifetime:** The `header` and `sequence` pointers are owned by
the `derep_session_s` and are valid until `derep_session_cleanup()` is
called. Copy them if you need them afterward.

### Functions

| Function | Description |
|----------|-------------|
| `derep_session_alloc()` | Allocate opaque session state. |
| `derep_session_free(ds)` | Free session state. Null-safe. |
| `derep_session_init(ds)` | Initialize session. No database required. |
| `derep_add_sequence(ds, header, seq, len, abund)` | Add one sequence. Normalized internally. |
| `derep_get_results(ds, results, max, count)` | Retrieve results sorted by abundance. Caller-provided array. |
| `derep_session_cleanup(ds)` | Free session resources. Invalidates result pointers. |

---

## Paired-end merging

Merge overlapping forward and reverse FASTQ reads into a single
consensus sequence with combined quality scores.

### Lifecycle

```cpp
// Initialize quality score lookup table (once)
mergepairs_init();

// Merge one pair
struct merge_result_s result;
int rc = mergepairs_single(
    fwd_seq, fwd_qual, fwd_len,
    rev_seq, rev_qual, rev_len,
    fwd_header, rev_header, &result);

if (rc == 0 && result.merged) {
    // result.merged_sequence, result.merged_quality, etc.
}
```

No per-thread state is needed. Each call to `mergepairs_single()` is
fully independent and thread-safe after `mergepairs_init()`.

### Result structure

| Field | Type | Description |
|-------|------|-------------|
| `merged` | `bool` | `true` if merge succeeded. |
| `merged_length` | `int` | Length of merged sequence. |
| `merged_sequence` | `char[10000]` | Merged DNA sequence (null-terminated). |
| `merged_quality` | `char[10000]` | Merged quality string (null-terminated). |
| `ee_merged` | `double` | Expected errors in merged sequence. |
| `ee_fwd` | `double` | Expected errors from forward read. |
| `ee_rev` | `double` | Expected errors from reverse read. |
| `fwd_errors` | `int` | Mismatches attributed to forward read. |
| `rev_errors` | `int` | Mismatches attributed to reverse read. |
| `overlap_length` | `int` | Length of overlap region. |

### Key options

| Global | Type | Default | Description |
|--------|------|---------|-------------|
| `opt_fastq_ascii` | `int64_t` | 33 | ASCII offset for quality scores. |
| `opt_fastq_minovlen` | `int64_t` | 10 | Minimum overlap length. |
| `opt_fastq_maxdiffs` | `int64_t` | 10 | Maximum mismatches in overlap. |
| `opt_fastq_maxdiffpct` | `double` | 100.0 | Maximum mismatch percentage in overlap. |
| `opt_fastq_maxee` | `double` | `DBL_MAX` | Maximum expected errors in merged read. |
| `opt_fastq_minmergelen` | `int64_t` | 0 | Minimum merged length. |
| `opt_fastq_maxmergelen` | `int64_t` | `INT_MAX` | Maximum merged length. |
| `opt_fastq_allowmergestagger` | `bool` | false | Allow staggered merges (read-through). |

### Functions

| Function | Description |
|----------|-------------|
| `mergepairs_init()` | Initialize quality lookup table. Call once before any merging. |
| `mergepairs_single(fwd_s, fwd_q, fwd_l, rev_s, rev_q, rev_l, fwd_h, rev_h, result)` | Merge one pair. Returns 0 on success, -1 on failure. Thread-safe. |

---

## Sequence masking

DUST low-complexity masking, available both as a batch operation on the
global database and as a standalone per-sequence function.

```cpp
// Batch: mask all sequences in the global database
dust_all();

// Standalone: mask a single sequence in-place (thread-safe)
char seq[] = "AAAAAAAAAAAACCGTACGT";
dust_single(seq, strlen(seq), false);  // false = soft-mask (lowercase)
```

`dust_single()` is fully thread-safe with no global state dependencies.
It can be called without any initialization.

### Functions

| Function | Description |
|----------|-------------|
| `dust_all()` | Apply DUST masking to all database sequences. |
| `dust_single(seq, len, hardmask)` | Mask one sequence in-place. `hardmask=true` replaces with N; `false` lowercases. Thread-safe. |

### Masking constants

| Constant | Value | Meaning |
|----------|-------|---------|
| `MASK_NONE` | 0 | No masking. |
| `MASK_DUST` | 1 | DUST low-complexity masking. |
| `MASK_SOFT` | 2 | Soft masking (lowercase). |

`opt_dbmask` and `opt_qmask` control which masking is applied to
database and query sequences, respectively. Both default to `MASK_DUST`.

---

## Configuration reference

All configuration is done by setting global `opt_*` variables after
`vsearch_init_defaults()` and before `vsearch_apply_defaults_fixups()`.
The full list of ~200 globals is declared as `extern` in `vsearch.h`.
Below are the most commonly used options, grouped by subsystem.

### Alignment scoring

| Global | Type | Default | Description |
|--------|------|---------|-------------|
| `opt_match` | `int64_t` | 2 | Match reward. |
| `opt_mismatch` | `int64_t` | -4 | Mismatch penalty. |
| `opt_gap_open_query_interior` | `int64_t` | 20 | Interior gap open (raw; adjusted by fixups). |
| `opt_gap_extension_query_interior` | `int64_t` | 2 | Interior gap extension. |
| `opt_gap_open_query_left` | `int64_t` | 2 | Left terminal gap open (raw). |
| `opt_gap_extension_query_left` | `int64_t` | 1 | Left terminal gap extension. |
| `opt_gap_open_query_right` | `int64_t` | 2 | Right terminal gap open (raw). |
| `opt_gap_extension_query_right` | `int64_t` | 1 | Right terminal gap extension. |

Target-side gap penalties mirror query-side with `_target_` in the name.
Defaults are identical.

Gap open penalties are stored as raw values. `vsearch_apply_defaults_fixups()`
subtracts the corresponding extension penalty (e.g., 20 - 2 = 18 for
interior). This matches the internal scoring convention.

### Search control

| Global | Type | Default | Description |
|--------|------|---------|-------------|
| `opt_id` | `double` | 0.0 | Minimum identity (0.0–1.0). |
| `opt_iddef` | `int64_t` | 2 | Identity definition (0=CD-HIT, 1=edit distance, 2=default, 3=marine bio, 4=BLAST). |
| `opt_maxaccepts` | `int64_t` | 1 | Stop after N accepted hits. |
| `opt_maxrejects` | `int64_t` | 32 | Stop after N rejected candidates. |
| `opt_maxhits` | `int64_t` | 0 | Maximum total hits (0 = unlimited, resolved by fixups). |
| `opt_wordlength` | `int64_t` | 8 | K-mer length for candidate selection. |
| `opt_minwordmatches` | `int64_t` | -1 | Minimum k-mer matches (-1 = auto, resolved by fixups). |
| `opt_strand` | `int64_t` | 1 | Strand: 1 = plus only, 2 = both. |

### Sequence filtering

| Global | Type | Default | Description |
|--------|------|---------|-------------|
| `opt_minseqlength` | `int64_t` | 0 | Minimum sequence length. |
| `opt_maxseqlength` | `int64_t` | `INT_MAX` | Maximum sequence length. |
| `opt_minsize` | `int64_t` | 0 | Minimum abundance. Set explicitly for your use case. |
| `opt_maxsize` | `int64_t` | `INT_MAX` | Maximum abundance. |

### Abundance annotation

| Global | Type | Default | Description |
|--------|------|---------|-------------|
| `opt_sizein` | `bool` | false | Parse `;size=N` from input headers. |
| `opt_sizeout` | `bool` | false | Write `;size=N` to output headers. |

### Output control

| Global | Type | Default | Description |
|--------|------|---------|-------------|
| `opt_quiet` | `bool` | true (library) | Suppress all stderr output. |
| `opt_no_progress` | `bool` | true (library) | Suppress progress bars. |
| `opt_threads` | `int64_t` | 0 | Number of threads (0 = auto-detect). |

---

## Error handling

vsearch was designed as a CLI tool and handles errors by printing a
message to stderr and calling `exit()` via the internal `fatal()`
function. **The library inherits this behavior.** On invalid input or
unrecoverable errors, the library will terminate the process.

Known fatal error triggers include:

- Out of memory (xmalloc/xrealloc failure)
- Invalid FASTA/FASTQ format in `db_read()`
- File I/O failures

**Library API return codes:**

| Function | Success | Failure |
|----------|---------|---------|
| `chimera_detect_single()` | 0 | (fatal on internal error) |
| `mergepairs_single()` | 0 | -1 (merge failed; result.merged = false) |
| All other functions | void | (fatal on internal error) |

**Recommendations for embedding:**

- Validate input data before passing to vsearch (DNA alphabet,
  non-empty sequences, reasonable lengths).
- If process termination is unacceptable, run vsearch operations in a
  child process or use signal handlers / `setjmp`/`longjmp` to catch
  `exit()` calls (fragile — may leak resources).
- With `opt_quiet = true` (the library default), no messages are
  written to stderr during normal operation. Fatal errors still write
  to stderr before exiting.

---

## Memory management

### Allocation and ownership

All result structures are **caller-allocated, library-populated**. The
library writes into caller-provided structs; no heap allocation is
performed for results.

```cpp
// Caller allocates on stack or heap
struct chimera_result_s result;
chimera_detect_single(ci, seq, head, len, abund, &result);
// result is fully populated — no cleanup needed
```

**Exception:** `derep_result_s` contains `const char *` pointers that
point into session-owned memory. These pointers are valid until
`derep_session_cleanup()` is called. Copy the strings if you need
them afterward.

### Opaque state objects

Per-thread and session objects are allocated by the library and must be
freed by the matching free function:

| Allocate | Free | Notes |
|----------|------|-------|
| `chimera_info_alloc()` | `chimera_info_free(ci)` | Call cleanup first. |
| `search_info_alloc()` | `search_info_free(si)` | Call cleanup first. |
| `cluster_session_alloc()` | `cluster_session_free(cs)` | Call cleanup first. |
| `derep_session_alloc()` | `derep_session_free(ds)` | Call cleanup first. |

All free functions are null-safe.

### Fixed-size buffers in result structs

Several result structs use fixed-size character buffers:

| Buffer | Size | Truncation |
|--------|------|------------|
| `chimera_result_s` labels | 1024 chars | Silent truncation to 1023 + null. |
| `cluster_result_s.centroid_label` | 1024 chars | Silent truncation. |
| `cluster_result_s.cigar` | 4096 chars | Truncation flagged via `cigar_truncated`. |
| `search_result_s.target_label` | 1024 chars | Silent truncation. |
| `merge_result_s` sequences | 10000 chars | Silent truncation. |

### Database pointers

Pointers returned by `db_getheader()`, `db_getsequence()`, and
`db_getquality()` point into global database memory. They are valid
only while the database is loaded (until `db_free()` is called).

### Internal allocations

`vsearch_init_defaults()` allocates `opt_ee_cutoffs_values` via
`xmalloc`. This is freed and reallocated on the next
`vsearch_init_defaults()` call. No action is required from the caller.

---

## Thread safety

### Summary

| Phase | Thread safety |
|-------|---------------|
| `vsearch_init_defaults()` | Single-threaded. Acquires session mutex. |
| Option overrides (`opt_* = ...`) | Single-threaded. |
| `vsearch_apply_defaults_fixups()` | Single-threaded. |
| Database loading (`db_init`, `db_add`, `dust_all`) | Single-threaded. |
| Index building (`dbindex_prepare`, `dbindex_addallsequences`) | Single-threaded. |
| Session init (`chimera_session_init`, etc.) | Single-threaded. |
| Per-thread init (`chimera_detect_thread_init`, etc.) | Safe for different instances. |
| Computation (`chimera_detect_single`, `search_single`, etc.) | Thread-safe with per-thread state. |
| Cleanup | Single-threaded. Join all threads first. |

### Rules

1. **One session at a time.** The session mutex prevents concurrent
   initialization but does not protect against misuse within a session.

2. **Each thread gets its own state object.** Never share
   `chimera_info_s`, `searchinfo_s`, or `cluster_session_s` across
   threads.

3. **Database is read-only during computation.** After indexing
   completes, the global database and k-mer index are safe to read
   from any thread.

4. **Dereplication is self-contained.** `derep_session_s` has no global
   state dependencies and can be used concurrently across sessions
   (each with its own instance).

5. **Merging is stateless.** After `mergepairs_init()`, calls to
   `mergepairs_single()` are fully independent and thread-safe.

6. **Masking:** `dust_single()` is thread-safe. `dust_all()` operates
   on the global database and is single-threaded.

---

## Platform support

### Supported targets

| Platform | Architecture | SIMD | Status |
|----------|-------------|------|--------|
| Linux | x86_64 | SSE2/SSSE3 (native) | Fully supported |
| Linux | aarch64 | NEON | Fully supported |
| Linux | ppc64le | AltiVec | Supported (upstream) |
| macOS | arm64 (native) | NEON | Supported |
| macOS | x86_64 (cross-compile) | SIMDe | Supported |
| Windows | x86_64 | — | Not supported (autotools) |
| WebAssembly | wasm32 | — | Not supported |

### SIMD architecture

vsearch uses architecture-specific SIMD for alignment:

- **x86_64:** Native SSE2 and SSSE3 intrinsics via separate compilation
  units (`libcpu_sse2.a`, `libcpu_ssse3.a`).
- **aarch64:** Native ARM NEON intrinsics.
- **ppc64le:** AltiVec intrinsics.
- **Other architectures:** Falls back to
  [SIMDe](https://github.com/simd-everywhere/simde) (vendored as
  `third_party/simde`), which provides portable SIMD emulation.

On macOS, the `-march=x86-64` flag is omitted (unsupported by Apple
Clang). SIMD instruction flags (`-msse2`, `-mssse3`) are applied when
`TARGET_X86_64` is true, regardless of Apple platform.

### Build configuration flags

| Flag | Effect |
|------|--------|
| `--disable-zlib` | Omit gzip support |
| `--disable-bzip2` | Omit bzip2 support |
| `--disable-pdfman` | Skip PDF manual generation |
| `--enable-debug` | Debug build with sanitizers and extra warnings |
| `--enable-profiling` | Profiling build (`-pg -O1`) |

---

## Examples

Working examples are in the `examples/` directory. Each demonstrates a
complete lifecycle: initialization, database loading, computation, and
cleanup.

| Example | Operation | Key features |
|---------|-----------|--------------|
| `example_chimera.cc` | Chimera detection (reference) | Single-threaded convenience API, 18-column output |
| `example_search.cc` | Global search | Identity filtering, multi-hit results |
| `example_cluster.cc` | Greedy clustering | Length-sorted DB, incremental centroid indexing |
| `example_cluster_strand.cc` | Clustering with strands | `opt_strand` behavior, reverse complement |
| `example_derep.cc` | Dereplication | Self-contained (no global DB), abundance output |
| `example_merge.cc` | Paired-end merging | FASTQ quality handling, stateless API |
| `example_dust.cc` | DUST masking | Standalone per-sequence, no init needed |
| `example_reinit.cc` | Re-initialization | Multi-session, multi-handle, gap penalty regression |

### Building examples

```bash
cd examples
g++ -std=c++11 -O2 -I../src -o example_chimera example_chimera.cc \
    ../src/libvsearch.a -lpthread -ldl
./example_chimera
```

### Minimal complete example

```cpp
#include "vsearch_api.h"
#include <cstdio>
#include <cstring>

int main() {
    // Initialize
    vsearch_init_defaults();
    opt_wordlength = 8;
    opt_id = 0.97;
    opt_maxaccepts = 1;
    vsearch_apply_defaults_fixups();

    // Load database
    db_init();
    const char * h = "ref1";
    const char * s = "ACGTACGTACGTACGTACGT";
    db_add(false, h, s, nullptr, strlen(h), strlen(s), 1);
    dust_all();
    dbindex_prepare(1, opt_dbmask);
    dbindex_addallsequences(opt_dbmask);

    // Search
    struct searchinfo_s * si = search_info_alloc();
    search_init(si);

    const char * qh = "query1";
    const char * qs = "ACGTACGTACGTACGTACGT";
    struct search_result_s results[5];
    int count;
    search_single(si, qs, qh, strlen(qs), 1, results, 5, &count);

    for (int i = 0; i < count; i++) {
        printf("%s\t%s\t%.1f%%\n", qh, results[i].target_label, results[i].id);
    }

    // Cleanup
    search_cleanup(si);
    search_info_free(si);
    dbindex_free();
    db_free();
    vsearch_session_end();
}
```
