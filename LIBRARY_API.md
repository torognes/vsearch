# vsearch Library API

vsearch can be embedded as a static library (`libvsearch.a`) in other
applications. This document provides comprehensive documentation of the
library API: build integration, initialization protocol, available
operations, result structures, error handling, memory ownership, thread
safety, and platform support.

## API version

The API version is exposed both at compile time and at runtime:

```cpp
#include "vsearch_api.h"

// Compile-time:
#if VSEARCH_API_VERSION < 1002000  // less than 1.2.0
#error "vsearch 1.2.0 or later required"
#endif

// Runtime:
int v = vsearch_api_version();          // MAJOR*1000000 + MINOR*1000 + PATCH
const char * s = vsearch_api_version_string();  // "MAJOR.MINOR.PATCH"
```

The numeric version (`VSEARCH_API_VERSION`) is encoded as
`MAJOR * 1000000 + MINOR * 1000 + PATCH` (same convention as OpenSSL
and libcurl), so each component may range `0..999` without collision.
The three components are also available individually as
`VSEARCH_API_VERSION_MAJOR`, `_MINOR`, and `_PATCH`.

### Version bump rules (semver)

- **MAJOR** — incompatible changes: removing, renaming, or retyping a
  public function; removing or reordering a field in a result struct;
  changing the semantics of an existing function. Pre-1.0
  (`MAJOR == 0`) the API is considered unstable and MINOR bumps may
  also break compatibility.
- **MINOR** — backward-compatible additions: new public functions,
  new fields appended to a result struct or to the `Parameters` struct.
- **PATCH** — backward-compatible fixes that don't change the API
  surface (bug fixes, doc updates, internal refactors).

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
# The generated build files are shipped and authoritative; no autoreconf needed.
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

The generated autotools files are committed to the repository and are
authoritative: maintainer mode is disabled (`AM_MAINTAINER_MODE`), so
`autoreconf` is not required at build time and `make` will not attempt to
regenerate them (no matching autoconf/automake version is needed, and
timestamps no longer matter). Run `./autogen.sh` (or `autoreconf -fi`)
only if you modify `configure.ac` or a `Makefile.am`.

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

vsearch was designed as a CLI tool with roughly 200 options. Configuration
now lives entirely in a `Parameters` struct threaded through the API — there
are no `opt_*` configuration globals. Internal compute state (the sequence
database and k-mer index) is still process-global, so there is no per-session
context object for it.

**Consequences:**

- Only one session can be active per process at a time
- Initialization and configuration are not thread-safe
- Each session is configured by a fresh `Parameters` struct

A session mutex serializes access: `vsearch_session_begin()` acquires
it, `vsearch_session_end()` releases it. A second `vsearch_session_begin()`
while a session is active fails with a fatal diagnostic (it does not block).

### Per-thread working state

Operations that support multi-threaded execution (chimera detection,
search, clustering) provide opaque per-thread state objects. Each
thread allocates its own instance; the objects must not be shared
across threads. The global database and k-mer index are read-only
after initialization, so no locking is needed during computation.

### Output and I/O

The library defaults to silent operation: `parameters.opt_quiet` and
`parameters.opt_no_progress` are both `true` in a default-constructed
`Parameters`. No output is written to stdout or stderr by library functions
under default settings. The caller has full control over I/O.

---

## Initialization

### Required sequence

Every session must follow this protocol:

```cpp
#include "vsearch_api.h"

// 1. A fresh, self-defaulting configuration (quiet, no progress).
struct Parameters parameters;

// 2. Override options for your use case.
parameters.opt_wordlength = 8;
parameters.opt_id = 0.97;

// 3. Begin the session: acquires the session mutex (fatal if one is already
//    active), resolves sentinel values, and applies the configuration.
vsearch_session_begin(parameters);

// 4-N. Use the library (load DB, run operations, etc.)
// ...

// Final. Release the session mutex.
vsearch_session_end();
```

`vsearch_session_begin()` resolves sentinel values that depend on other
options (via `vsearch_apply_defaults_fixups(parameters)`):

- `opt_maxhits`: 0 → `INT64_MAX`
- `opt_minwordmatches`: -1 → wordlength-based default from
  `searchcore.h` lookup table
- Gap penalties: raw values are adjusted by subtracting extension
  penalties (e.g., `opt_gap_open_query_interior` 20 → 18 after
  subtracting `opt_gap_extension_query_interior` 2)

You configure only the fields you care about; the rest keep their library
defaults. The `Parameters` struct is the single configuration source — the
compute engines read it directly, so there are no `opt_*` globals to set.

### Re-initialization

Multiple sequential sessions in the same process are supported.
Repeat the full sequence (fresh `Parameters` → configure →
`vsearch_session_begin` → work → cleanup → `vsearch_session_end`) for each
session. Each fresh struct re-applies the gap penalty adjustments from raw
defaults.

See `api_examples/example_reinit.cc` for a tested multi-session example.

### Session functions

| Function | Description |
|----------|-------------|
| `vsearch_session_begin(Parameters &)` | Acquire session mutex, resolve sentinels, apply config. Call once after configuring. |
| `vsearch_apply_defaults_fixups(Parameters &)` | Resolve a struct's sentinel values (called by session_begin; exposed for inspection). |
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
db_read("sequences.fasta", 0, parameters);  // 0 = do not upcase
```

### Masking and indexing

After loading sequences, apply masking and build the k-mer index:

```cpp
dust_all();                                       // DUST low-complexity masking
Dbindex dbindex;                                  // the k-mer index (owns its buffers, RAII)
dbindex.prepare(1, opt_dbmask, parameters);       // allocate index (1 = use bitmap)
dbindex.add_all_sequences(opt_dbmask, parameters);// index all sequences
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
// the Dbindex frees its buffers automatically when it goes out of scope (RAII);
// call clear() explicitly only to release it early or to reuse it for another build.
db_free();
```

### Database functions

| Function | Description |
|----------|-------------|
| `db_init()` | Reset database state. Call before `db_add()`. |
| `db_add(is_fastq, header, seq, qual, hlen, slen, abund)` | Add one sequence. |
| `db_read(filename, upcase, parameters)` | Load sequences from file. |
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
| `Dbindex dbindex;` | The k-mer index, an owned object (RAII, non-copyable). Pass it by reference to the session/batch entry points below. |
| `dbindex.prepare(bitmap, mask, parameters)` | Allocate k-mer index. `bitmap=1` enables bitmap mode (required for clustering). |
| `dbindex.add_all_sequences(mask, parameters)` | Index all loaded sequences. |
| `dbindex.add_sequence(seqno, mask)` | Index a single sequence (for incremental indexing). |
| `dbindex.clear()` | Free k-mer index memory (also done automatically by the destructor). |

---

## Chimera detection

Reference-based chimera detection (UCHIME algorithm: Edgar et al. 2011,
Bioinformatics 27:2194-2200) with vsearch's SIMD-optimized
Needleman-Wunsch alignment.

### Lifecycle

**Single-threaded** (convenience wrappers):

```cpp
struct chimera_info_s * ci = chimera_info_alloc();
chimera_detect_init(ci, parameters, dbindex);   // session + per-thread init

struct chimera_result_s result;
chimera_detect_single(ci, query_seq, query_head,
                      query_len, query_abundance, &result);

chimera_detect_cleanup(ci);    // per-thread + session cleanup
chimera_info_free(ci);
```

**Multi-threaded** (split API):

```cpp
// Session init: once, single-threaded, after DB indexed
chimera_session_init(parameters);

// Per-thread init: one handle per thread
struct chimera_info_s * ci1 = chimera_info_alloc();
chimera_detect_thread_init(ci1, parameters, dbindex);
struct chimera_info_s * ci2 = chimera_info_alloc();
chimera_detect_thread_init(ci2, parameters, dbindex);

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

For de novo chimera detection, set `opt_chimeras_denovo` to a non-null
value before calling `chimera_detect_init()`. This selects the
`--chimeras_denovo` algorithm (detection in long exact sequences), which
is distinct from `--uchime_denovo` — the two produce different
classifications and output. In this mode:

- Process queries in decreasing abundance order
- After classifying a query as non-chimeric, add it to the reference
  database with `db_add()` and index it with `dbindex.add_sequence()`
- The `query_abundance` parameter to `chimera_detect_single()` controls
  abundance skew filtering

De novo mode is inherently sequential (single-threaded).

### Key options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `opt_wordlength` | `int64_t` | 8 | K-mer length. Set before `vsearch_session_begin()`. |
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
| `chimera_session_init(parameters)` | Session-level init hook. Currently a no-op — detection config is built per-thread in `chimera_detect_thread_init` — but kept as a stable API symbol; call once after DB indexed. |
| `chimera_session_cleanup()` | Session-level teardown: destroy mutexes. Call after all per-thread cleanup. |
| `chimera_detect_thread_init(ci, parameters, dbindex)` | Per-thread init: allocate SIMD aligners, k-mer finders; searches `dbindex`. |
| `chimera_detect_thread_cleanup(ci)` | Per-thread teardown: free resources. |
| `chimera_detect_single(ci, seq, head, len, abund, result)` | Detect chimera for one query. Returns 0 on success. |
| `chimera_detect_init(ci, parameters, dbindex)` | Convenience: `session_init` + `thread_init`. Single-threaded only. |
| `chimera_detect_cleanup(ci)` | Convenience: `thread_cleanup` + `session_cleanup`. Single-threaded only. |

---

## Global search

Search query sequences against the loaded database, reporting hits
above a minimum identity threshold.

### Lifecycle

```cpp
// Configure search parameters on the struct, before vsearch_session_begin()
parameters.opt_id = 0.97;           // minimum 97% identity
parameters.opt_maxaccepts = 3;      // return up to 3 hits
parameters.opt_maxrejects = 16;     // give up after 16 rejects
parameters.opt_strand = true;       // optional: search both strands
// ... vsearch_session_begin(parameters), load and index the DB ...

// Allocate and initialize the session
struct search_session_s * ss = search_session_alloc();
search_session_init(ss, parameters, dbindex);  // call after DB indexed

// Search queries
struct search_result_s results[10];
int result_count;
search_session_single(ss, query_seq, query_head, query_len,
                      query_abundance, results, 10, &result_count);

// Target headers are looked up from the database by index
for (int i = 0; i < result_count; i++) {
    printf("%s\t%s\t%.1f\n",
           query_head,
           db_getheader(results[i].target),
           results[i].id);
}

// Cleanup
search_session_cleanup(ss);
search_session_free(ss);
```

For bulk workloads, `search_batch()` parallelizes across `opt_threads`:

```cpp
search_batch(parameters, dbindex, q_seqs, q_heads, q_lens, q_sizes, query_count,
             results, max_results_per_query, result_counts);
```

### Result structure

`search_result_s` contains per-hit alignment details:

| Field | Type | Description |
|-------|------|-------------|
| `target` | `int` | Database sequence index. Use `db_getheader(target)` for the header string and `db_getsequence(target)` for the sequence. |
| `id` | `double` | Percent identity (per `opt_iddef`). |
| `matches` | `int` | Matching columns. |
| `mismatches` | `int` | Mismatching columns. |
| `gaps` | `int` | Gap columns. |
| `alignment_length` | `int` | Total alignment length. |
| `query_length` | `int` | Query sequence length. |
| `target_length` | `int` | Target sequence length. |
| `accepted` | `bool` | Whether hit passed the identity threshold. |
| `strand` | `int` | `0` = plus strand, `1` = minus strand (when `opt_strand` is true). |

Results are ordered by identity (descending). The caller provides the
result array and its capacity; `result_count` is set to the number of
hits returned (up to `max_results` and `opt_maxaccepts`).

### Key options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `opt_id` | `double` | 0.0 | Minimum identity threshold (0.0–1.0). |
| `opt_maxaccepts` | `int64_t` | 1 | Stop after N accepted hits. |
| `opt_maxrejects` | `int64_t` | 32 | Stop after N rejected candidates. |
| `opt_iddef` | `int64_t` | 2 | Identity definition (0–4). |
| `opt_wordlength` | `int64_t` | 8 | K-mer length for candidate selection. |
| `opt_strand` | `bool` | `false` | `false` = plus strand only, `true` = both strands. |

### Functions

| Function | Description |
|----------|-------------|
| `search_session_alloc()` | Allocate opaque session state. |
| `search_session_free(ss)` | Free session state. Null-safe (cleanup is implicit). |
| `search_session_init(ss, parameters, dbindex)` | Initialize session. Call after DB indexed. Respects `opt_strand`. Stores a reference to `dbindex`, which must outlive the session. |
| `search_session_single(ss, seq, head, len, size, results, max, count)` | Search one query (both strands when `opt_strand` is true). One session per process; do not share across threads. |
| `search_session_cleanup(ss)` | Free per-session resources. Call before `search_session_free`. |
| `search_batch(parameters, dbindex, seqs, heads, lens, sizes, n, results, max_per, counts)` | Bulk-parallel search of `dbindex`. Internally uses `opt_threads`. |

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

// Prepare index — do NOT call dbindex.add_all_sequences().
// Centroids are indexed incrementally during clustering.
Dbindex dbindex;
dbindex.prepare(1, opt_qmask, parameters);

// Allocate and initialize session
struct cluster_session_s * cs = cluster_session_alloc();
cluster_session_init(cs, parameters, dbindex);

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

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `opt_id` | `double` | 0.0 | Identity threshold for cluster membership. |
| `opt_strand` | `bool` | `false` | `false` = plus strand only, `true` = both strands. |
| `opt_maxaccepts` | `int64_t` | 1 | Accepts before stopping search. |
| `opt_maxrejects` | `int64_t` | 32 | Rejects before stopping search. |

### Functions

| Function | Description |
|----------|-------------|
| `cluster_session_alloc()` | Allocate opaque session state. |
| `cluster_session_free(cs)` | Free session state. Null-safe. |
| `cluster_session_init(cs, parameters, dbindex)` | Initialize session. DB must be sorted; `dbindex.prepare()` called but NOT `add_all_sequences`. Stores a reference to `dbindex` (mutated as centroids are added); it must outlive the session. |
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
mergepairs_init(parameters);

// Merge one pair. Zero-initialize so the pointers start at nullptr.
struct merge_result_s result = {};
int rc = mergepairs_single(
    parameters,
    fwd_seq, fwd_qual, fwd_len,
    rev_seq, rev_qual, rev_len,
    fwd_header, rev_header, &result);

if (rc == 0 && result.merged) {
    // use result.merged_sequence / result.merged_quality
    //   (null-terminated, length == result.merged_length)
}

// Release caller-owned buffers before reusing or discarding result
merge_result_free(&result);
```

No per-thread state is needed. Each call to `mergepairs_single()` is
fully independent and thread-safe after `mergepairs_init()`.

### Result structure

| Field | Type | Description |
|-------|------|-------------|
| `merged` | `bool` | `true` if merge succeeded. |
| `merged_length` | `int` | Length of merged sequence. |
| `merged_sequence` | `char *` | `xmalloc`'d merged DNA (null-terminated), sized to `merged_length + 1`. `nullptr` on failure or before first call. Caller must release via `merge_result_free()`. |
| `merged_quality` | `char *` | `xmalloc`'d quality (null-terminated). Ownership matches `merged_sequence`. |
| `ee_merged` | `double` | Expected errors in merged sequence. |
| `ee_fwd` | `double` | Expected errors from forward read. |
| `ee_rev` | `double` | Expected errors from reverse read. |
| `fwd_errors` | `int` | Mismatches attributed to forward read. |
| `rev_errors` | `int` | Mismatches attributed to reverse read. |
| `overlap_length` | `int` | Length of overlap region. |

Buffers grow as needed — there is no fixed upper bound on merged
length. Each call to `mergepairs_single()` overwrites the pointer
fields without freeing them, so reuse the same struct only after
calling `merge_result_free()` between calls.

### Key options

| Option | Type | Default | Description |
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
| `mergepairs_init(parameters)` | Initialize quality lookup table. Call once before any merging. |
| `mergepairs_single(parameters, fwd_s, fwd_q, fwd_l, rev_s, rev_q, rev_l, fwd_h, rev_h, result)` | Merge one pair. Allocates `result->merged_sequence` / `merged_quality` via xmalloc. Returns 0 on success, -1 on failure. Thread-safe. |
| `merge_result_free(result)` | Free the merged sequence/quality buffers and null the pointers. Null-safe on either field. |

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

### Masking modes

`opt_dbmask` and `opt_qmask` are of type `Masking`, an `enum struct`
declared in `mask.h`, and control which masking is applied to database
and query sequences, respectively.

| Enumerator | Value | Meaning |
|------------|-------|---------|
| `Masking::none` | 0 | No masking. |
| `Masking::dust` | 1 | DUST low-complexity masking. |
| `Masking::soft` | 2 | Soft masking (lowercase). |

Both default to `Masking::dust`. (`Masking::error`, value -1, is an
internal sentinel produced when an unrecognised `--qmask` / `--dbmask`
argument is parsed; it is not a valid masking mode.)

---

## Configuration reference

All configuration is done by setting `opt_*` fields on a `Parameters` struct
before `vsearch_session_begin()`. The full set of ~200 options is defined in
the `Parameters` struct in `vsearch.h`. The tables below list them by option
name (accessed as `parameters.opt_name`), grouped by subsystem.

### Alignment scoring

| Option | Type | Default | Description |
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

Gap open penalties are stored as raw values. `vsearch_session_begin()`
subtracts the corresponding extension penalty (e.g., 20 - 2 = 18 for
interior). This matches the internal scoring convention.

### Search control

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `opt_id` | `double` | 0.0 | Minimum identity (0.0–1.0). |
| `opt_iddef` | `int64_t` | 2 | Identity definition (0=CD-HIT, 1=edit distance, 2=default, 3=marine bio, 4=BLAST). |
| `opt_maxaccepts` | `int64_t` | 1 | Stop after N accepted hits. |
| `opt_maxrejects` | `int64_t` | 32 | Stop after N rejected candidates. |
| `opt_maxhits` | `int64_t` | 0 | Maximum total hits (0 = unlimited, resolved by fixups). |
| `opt_wordlength` | `int64_t` | 8 | K-mer length for candidate selection. |
| `opt_minwordmatches` | `int64_t` | -1 | Minimum k-mer matches (-1 = auto, resolved by fixups). |
| `opt_strand` | `bool` | `false` | Strand: `false` = plus only, `true` = both. |

### Sequence filtering

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `opt_minseqlength` | `int64_t` | -1 | Minimum sequence length. The default is the -1 "unset" sentinel; `db_read()` treats any non-positive value as no lower bound (the CLI resolves it to a command-specific 1 or 32). |
| `opt_maxseqlength` | `int64_t` | `INT_MAX` | Maximum sequence length. |
| `opt_minsize` | `int64_t` | 0 | Minimum abundance. Set explicitly for your use case. |
| `opt_maxsize` | `int64_t` | `INT_MAX` | Maximum abundance. |

### Abundance annotation

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `opt_sizein` | `bool` | false | Parse `;size=N` from input headers. |
| `opt_sizeout` | `bool` | false | Write `;size=N` to output headers. |

### Output control

| Option | Type | Default | Description |
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
| `search_session_alloc()` | `search_session_free(ss)` | Call cleanup first. |
| `cluster_session_alloc()` | `cluster_session_free(cs)` | Call cleanup first. |
| `derep_session_alloc()` | `derep_session_free(ds)` | Call cleanup first. |

All free functions are null-safe.

### Fixed-size buffers in result structs

Several result structs still use fixed-size character buffers:

| Buffer | Size | Truncation |
|--------|------|------------|
| `chimera_result_s` labels | 1024 chars | Silent truncation to 1023 + null. |
| `cluster_result_s.centroid_label` | 1024 chars | Silent truncation. |
| `cluster_result_s.cigar` | 4096 chars | Truncation flagged via `cigar_truncated`. |

`merge_result_s.merged_sequence` / `merged_quality` are dynamically
allocated by `mergepairs_single()` and owned by the caller — release
with `merge_result_free()`. `search_result_s` no longer carries a
target header copy; look it up with `db_getheader(result.target)`.

### Database pointers

Pointers returned by `db_getheader()`, `db_getsequence()`, and
`db_getquality()` point into global database memory. They are valid
only while the database is loaded (until `db_free()` is called).

### Internal allocations

`vsearch_session_begin()` derives the internal `opt_ee_cutoffs_values` array
from `parameters.opt_ee_cutoffs`, allocating it via `xmalloc` and freeing any
previous allocation. No action is required from the caller.

---

## Thread safety

### Summary

| Phase | Thread safety |
|-------|---------------|
| Configuring `Parameters` (`parameters.opt_* = ...`) | Single-threaded. |
| `vsearch_session_begin()` | Single-threaded. Acquires session mutex. |
| Database loading (`db_init`, `db_add`, `dust_all`) | Single-threaded. |
| Index building (`dbindex.prepare`, `dbindex.add_all_sequences`) | Single-threaded. |
| Session init (`chimera_session_init`, etc.) | Single-threaded. |
| Per-thread init (`chimera_detect_thread_init`, etc.) | Safe for different instances. |
| Computation (`chimera_detect_single`, `search_session_single`, etc.) | Thread-safe with per-thread state. |
| Cleanup | Single-threaded. Join all threads first. |

### Rules

1. **One session at a time.** The session mutex prevents concurrent
   initialization but does not protect against misuse within a session.

2. **Each thread gets its own state object.** Never share
   `chimera_info_s`, `search_session_s`, or `cluster_session_s` across
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

Working examples are in the `api_examples/` directory. Each demonstrates a
complete lifecycle: initialization, database loading, computation, and
cleanup.

| Example | Operation | Key features |
|---------|-----------|--------------|
| `example_chimera.cc` | Chimera detection (reference) | Single-threaded convenience API, 18-column output |
| `example_search.cc` | Global search | Identity filtering, multi-hit results |
| `example_cluster.cc` | Greedy clustering | Length-sorted DB, incremental centroid indexing |
| `example_derep.cc` | Dereplication | Self-contained (no global DB), abundance output |
| `example_merge.cc` | Paired-end merging | FASTQ quality handling, stateless API |
| `example_dust.cc` | DUST masking | Standalone per-sequence, no init needed |
| `example_reinit.cc` | Re-initialization | Multi-session, multi-handle, gap penalty regression |
| `example_lifecycle.cc` | API memory/error contracts | Free null-safety, `merge_result_free` idempotency, `mergepairs_single` -1 return, result reuse |
| `example_dbinfo.cc` | Database query/indexing surface | `db_read`, statistical accessors, quality, sort ordering, incremental `dbindex.add_sequence` |

### Building examples

```bash
cd api_examples
g++ -std=c++11 -O3 -I../src -o example_chimera example_chimera.cc \
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
    struct Parameters parameters;
    parameters.opt_wordlength = 8;
    parameters.opt_id = 0.97;
    parameters.opt_maxaccepts = 1;
    vsearch_session_begin(parameters);

    // Load database
    db_init();
    const char * h = "ref1";
    const char * s = "ACGTACGTACGTACGTACGT";
    db_add(false, h, s, nullptr, strlen(h), strlen(s), 1);
    dust_all();
    Dbindex dbindex;
    dbindex.prepare(1, opt_dbmask, parameters);
    dbindex.add_all_sequences(opt_dbmask, parameters);

    // Search
    struct search_session_s * ss = search_session_alloc();
    search_session_init(ss, parameters, dbindex);

    const char * qh = "query1";
    const char * qs = "ACGTACGTACGTACGTACGT";
    struct search_result_s results[5];
    int count;
    search_session_single(ss, qs, qh, strlen(qs), 1, results, 5, &count);

    for (int i = 0; i < count; i++) {
        printf("%s\t%s\t%.1f%%\n",
               qh, db_getheader(results[i].target), results[i].id);
    }

    // Cleanup
    search_session_cleanup(ss);
    search_session_free(ss);
    dbindex.clear();  // (also freed automatically at scope exit)
    db_free();
    vsearch_session_end();
}
```
