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

This transitively includes `vsearch.hpp` (all global declarations) and
the module headers for each API subsystem.

### C++ standard

The library requires C++11 or later. Build with `-std=c++11` (or
higher).

---

## Architecture overview

### Global state model

vsearch was designed as a CLI tool with roughly 200 options. Configuration
now lives entirely in a `Parameters` struct threaded through the API — there
are no `opt_*` configuration globals. The sequence database and k-mer index
are now caller-owned objects (`Database` and `Dbindex`) rather than process
globals. Some lower-level compute state is still process-global (e.g. the
`showalign` alignment row buffers and the `derep_prefix` sort comparator
bridge), so there is not yet a full per-session context object and only
one session may be active at a time.

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
across threads. The database and k-mer index are read-only
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

// 4. An owned, empty database (RAII; frees itself at scope exit). Populate it
//    with db.read(file, ...) or db.init() + db.add(), then mask and index it.
Database db;

// 5-N. Use the library (load DB, run operations, etc.)
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

### RAII session guard (C++)

C++ callers can replace the explicit `vsearch_session_begin()` /
`vsearch_session_end()` pair with the `VsearchSession` guard declared in
`vsearch_api.h`: its constructor begins the session and its destructor ends
it, so the session is released automatically at scope exit — including on an
early `return`. It is non-copyable and non-movable (a session is a single
process-wide resource).

```cpp
#include "vsearch_api.h"

struct Parameters parameters;
parameters.opt_wordlength = 8;
parameters.opt_id = 0.97;

{
  VsearchSession const session(parameters);   // begins the session here

  Database db;
  // ... configure, load, and use the library ...

}  // session ends here: vsearch_session_end() runs in the destructor
```

The guard is a thin convenience wrapper; the two functions remain the primary
interface, so code that needs an unscoped begin/end pair (or a C-style
interface) can keep calling them directly. Most programs in `api_examples/`
use the guard; `example_reinit.cc` uses the explicit calls.

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
| `VsearchSession session(Parameters &)` | RAII guard (C++): begins the session in its constructor, ends it in its destructor at scope exit. Optional convenience over the two functions above; non-copyable and non-movable. |

---

## Database management

vsearch stores sequences in a caller-owned in-memory database, `Database db;`
(a non-copyable RAII object declared in `core/db.hpp` that frees itself when it
goes out of scope). Sequences are loaded via `db.init()` + `db.add()`, then
indexed for k-mer search.

### Loading sequences

```cpp
// The owned database (RAII; frees itself at scope exit)
Database db;

// Reset database state (also frees any previous data)
db.init();

// Add sequences one at a time. The header/sequence/quality are passed as a
// SeqRecord of non-owning View<char> windows (pointer + length); the views need
// not be NUL-terminated. For FASTA the quality view is empty ({}).
for (int i = 0; i < n_seqs; i++) {
    db.add(false,             // is_fastq: false for FASTA
           SeqRecord{View<char>{headers[i], strlen(headers[i])},
                     View<char>{sequences[i], strlen(sequences[i])},
                     View<char>{}},   // quality view empty for FASTA
           abundance);        // size annotation (1 for uniform abundance)
}
// A FASTQ record passes a non-empty quality view (same length as the sequence)
// and is_fastq = true:
//     db.add(true, SeqRecord{View<char>{header, hlen},
//                            View<char>{seq, slen},
//                            View<char>{qual, slen}}, abundance);
```

**Requirements:**

- `db.init()` must be called before `db.add()`. Without it, internal
  statistics (shortest sequence length) are corrupted.
- Sequences must be DNA (ACGT). Convert U to T before calling.
- `db.init()` clears any previous data internally (like `db.clear()`), so
  it is safe to call on an already-populated database.

### Reading from files

For file-based workflows, `db.read()` loads sequences from FASTA/FASTQ
files directly:

```cpp
db.read("sequences.fasta", 0, parameters);  // 0 = do not upcase
```

### Masking and indexing

After loading sequences, apply masking and build the k-mer index:

```cpp
dust_all(db, parameters);                              // DUST low-complexity masking
Dbindex dbindex;                                       // the k-mer index (owns its buffers, RAII)
dbindex.prepare(1, opt_dbmask, db, parameters);        // allocate index (1 = use bitmap)
dbindex.add_all_sequences(opt_dbmask, db, parameters); // index all sequences
```

### Sorting

Some operations require the database to be pre-sorted:

```cpp
db.sortbylength(parameters);               // for cluster_fast
db.sortbylength_shortest_first(parameters);
db.sortbyabundance(parameters);            // for cluster_size
```

### Cleanup

```cpp
// Both the Database and the Dbindex free their buffers automatically when they
// go out of scope (RAII); call clear() explicitly only to release memory early
// or to reuse the object for another build.
db.clear();
```

### Database functions

| Function | Description |
|----------|-------------|
| `db.init()` | Reset database state. Call before `db.add()`. |
| `db.add(is_fastq, SeqRecord{header, seq, qual}, abund)` | Add one sequence (header/seq/qual as View<char> windows; empty quality view for FASTA). |
| `db.read(filename, upcase, parameters)` | Load sequences from file. |
| `db.clear()` | Free all database memory (also done automatically by the destructor). |
| `db.getsequencecount()` | Number of loaded sequences. |
| `db.getnucleotidecount()` | Total nucleotides across all sequences. |
| `db.getlongestsequence()` | Length of longest sequence. |
| `db.getshortestsequence()` | Length of shortest sequence. |
| `db.getlongestheader()` | Length of longest header. |
| `db.getheader(seqno)` | Read-only pointer to header string (database-owned). |
| `db.getsequence(seqno)` | Read-only pointer to sequence string (database-owned). |
| `db.getsequencelen(seqno)` | Sequence length. |
| `db.getheaderlen(seqno)` | Header length. |
| `db.getabundance(seqno)` | Abundance annotation. |
| `db.getquality(seqno)` | Quality string (FASTQ only). |
| `db.is_fastq()` | Whether database contains FASTQ data (accessor). |
| `db.sortbylength(parameters)` | Sort by length, longest first. |
| `db.sortbylength_shortest_first(parameters)` | Sort by length, shortest first. |
| `db.sortbyabundance(parameters)` | Sort by abundance, highest first. |

### K-mer index functions

| Function | Description |
|----------|-------------|
| `Dbindex dbindex;` | The k-mer index, an owned object (RAII, non-copyable). Pass it by reference to the session/batch entry points below. |
| `dbindex.prepare(bitmap, mask, db, parameters)` | Allocate k-mer index. `bitmap=1` enables bitmap mode (required for clustering). |
| `dbindex.add_all_sequences(mask, db, parameters)` | Index all loaded sequences. |
| `dbindex.add_sequence(seqno, mask, db)` | Index a single sequence (for incremental indexing). |
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
chimera_detect_init(ci, parameters, dbindex, db);   // session + per-thread init

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
chimera_detect_thread_init(ci1, parameters, dbindex, db);
struct chimera_info_s * ci2 = chimera_info_alloc();
chimera_detect_thread_init(ci2, parameters, dbindex, db);

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
  database with `db.add()` and index it with `dbindex.add_sequence()`
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
| `chimera_detect_thread_init(ci, parameters, dbindex, db)` | Per-thread init: allocate SIMD aligners, k-mer finders; searches `dbindex`. |
| `chimera_detect_thread_cleanup(ci)` | Per-thread teardown: free resources. |
| `chimera_detect_single(ci, seq, head, len, abund, result)` | Detect chimera for one query. Returns 0 on success. |
| `chimera_detect_init(ci, parameters, dbindex, db)` | Convenience: `session_init` + `thread_init`. Single-threaded only. |
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
search_session_init(ss, parameters, dbindex, db);  // call after DB indexed

// Search queries
struct search_result_s results[10];
int result_count;
search_session_single(ss, query_seq, query_head, query_len,
                      query_abundance, results, 10, &result_count);

// Target headers are looked up from the database by index
for (int i = 0; i < result_count; i++) {
    printf("%s\t%s\t%.1f\n",
           query_head,
           db.getheader(results[i].target),
           results[i].id);
}

// Cleanup
search_session_cleanup(ss);
search_session_free(ss);
```

For bulk workloads, `search_batch()` parallelizes across `opt_threads`:

```cpp
search_batch(parameters, dbindex, db, q_seqs, q_heads, q_lens, q_sizes, query_count,
             results, max_results_per_query, result_counts);
```

### Result structure

`search_result_s` contains per-hit alignment details:

| Field | Type | Description |
|-------|------|-------------|
| `target` | `int` | Database sequence index. Use `db.getheader(target)` for the header string and `db.getsequence(target)` for the sequence. |
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
| `search_session_init(ss, parameters, dbindex, db)` | Initialize session. Call after DB indexed. Respects `opt_strand`. Stores a reference to `dbindex`, which must outlive the session. |
| `search_session_single(ss, seq, head, len, size, results, max, count)` | Search one query (both strands when `opt_strand` is true). One session per process; do not share across threads. |
| `search_session_cleanup(ss)` | Free per-session resources. Call before `search_session_free`. |
| `search_batch(parameters, dbindex, db, seqs, heads, lens, sizes, n, results, max_per, counts)` | Bulk-parallel search of `dbindex`. Internally uses `opt_threads`. |

---

## Clustering

Greedy centroid-based clustering (equivalent to `--cluster_fast` or
`--cluster_size`). Sequences are processed sequentially; each is either
assigned to an existing cluster or becomes a new centroid.

### Lifecycle

```cpp
// Sort database (required)
db.sortbylength(parameters);          // for cluster_fast behavior
// or: db.sortbyabundance(parameters) for cluster_size behavior

// Prepare index — do NOT call dbindex.add_all_sequences().
// Centroids are indexed incrementally during clustering.
Dbindex dbindex;
dbindex.prepare(1, opt_qmask, db, parameters);

// Allocate and initialize session
struct cluster_session_s * cs = cluster_session_alloc();
cluster_session_init(cs, parameters, dbindex, db);

// Assign sequences sequentially (0, 1, 2, ...)
struct cluster_result_s result;
for (int i = 0; i < db.getsequencecount(); i++) {
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
| `cluster_session_init(cs, parameters, dbindex, db)` | Initialize session. DB must be sorted; `dbindex.prepare()` called but NOT `add_all_sequences`. Stores a reference to `dbindex` (mutated as centroids are added); it must outlive the session. |
| `cluster_assign_single(cs, seqno, result)` | Assign one sequence. Must be called in seqno order (0, 1, 2, ...). |
| `cluster_session_cleanup(cs)` | Free session resources. |

---

## Dereplication

Full-length dereplication: collapse identical sequences and sum their
abundances. This operation does NOT use the database — it
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
// Create a merge session (once); it builds and privately holds the quality
// score lookup tables and is reused for every merge.
MergePairs const merger(parameters);

// Merge one pair. sequence and quality are read-only views (same length).
MergeResult const result = merger.merge(
    parameters,
    MergeInput{View<char>{fwd_seq, fwd_len}, View<char>{fwd_qual, fwd_len}},
    MergeInput{View<char>{rev_seq, rev_len}, View<char>{rev_qual, rev_len}});

if (result.merged) {
    // use result.sequence / result.quality (std::string,
    //   length == result.merged_length); freed with the result (RAII)
}
```

No per-thread state is needed. A `MergePairs` session is read-only once
constructed, so a single instance can be shared across threads; `merge()`
is `const` and every call is fully independent and thread-safe. The returned
`MergeResult` owns its buffers and frees them automatically when it goes out
of scope — there is nothing for the caller to release.

### Result structure

| Field | Type | Description |
|-------|------|-------------|
| `merged` | `bool` | `true` if merge succeeded. |
| `merged_length` | `int` | Length of merged sequence. |
| `sequence` | `std::string` | Merged DNA, `merged_length` characters. Empty on failure. |
| `quality` | `std::string` | Merged quality (ASCII). Empty on failure. Same length as `sequence`. |
| `ee_merged` | `double` | Expected errors in merged sequence. |
| `ee_fwd` | `double` | Expected errors from forward read. |
| `ee_rev` | `double` | Expected errors from reverse read. |
| `fwd_errors` | `int` | Mismatches attributed to forward read. |
| `rev_errors` | `int` | Mismatches attributed to reverse read. |
| `overlap_length` | `int` | Length of overlap region. |
| `error` | `MergeError` | Hard input error, if any: `none` (success or ordinary non-merge), `quality_below_qmin`, or `quality_above_qmax`. Added in API 0.16.0. |
| `error_value` | `int` | The offending FASTQ quality value when `error != none`. Added in API 0.16.0. |

The strings grow as needed — there is no fixed upper bound on merged
length — and are freed with the `MergeResult` (RAII). `merge()` returns a
fresh result by value each call, so results are always independent.

A non-merge (`merged == false`) is ordinary — poor overlap, too many
differences, etc. — unless `error != MergeError::none`, which flags a hard
input error: a FASTQ quality value outside `[fastq_qmin, fastq_qmax]`
(`error_value` carries it). The CLI treats that as fatal; the library reports
it on the result instead of throwing, so a caller can skip the pair or stop the
batch as it chooses.

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
| `MergePairs(parameters)` | Construct a merge session. Builds the quality lookup tables once and holds them privately; reuse the session for every merge. |
| `MergePairs::merge(parameters, fwd, rev)` | Merge one pair (`fwd`/`rev` are `MergeInput` = sequence + equal-length quality `View<char>`). Returns a `MergeResult` by value; on failure `result.merged` is `false` and the strings are empty (`result.error` distinguishes an out-of-range quality from an ordinary non-merge). `const`, thread-safe. |
| `MergeInput{sequence, quality}` | One read: two read-only `View<char>` of equal length (one quality symbol per base). |

---

## Sequence masking

DUST low-complexity masking, available both as a batch operation on the
database and as a standalone per-sequence function.

```cpp
// Batch: mask all sequences in the database
dust_all(db, parameters);

// Standalone: mask a single sequence in-place (thread-safe)
char seq[] = "AAAAAAAAAAAACCGTACGT";
dust_single(seq, strlen(seq), false);  // false = soft-mask (lowercase)
```

`dust_single()` is fully thread-safe with no global state dependencies.
It can be called without any initialization.

### Functions

| Function | Description |
|----------|-------------|
| `dust_all(db, parameters)` | Apply DUST masking to all database sequences. |
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
the `Parameters` struct in `vsearch.hpp`. The tables below list them by option
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
| `opt_minseqlength` | `int64_t` | -1 | Minimum sequence length. The default is the -1 "unset" sentinel; `db.read()` treats any non-positive value as no lower bound (the CLI resolves it to a command-specific 1 or 32). |
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

vsearch was designed as a CLI tool: internally, an unrecoverable error
prints a message to stderr and ends the program through the `fatal()`
function. In the **CLI** that means `std::exit()`. In the **library**,
`fatal()` instead **throws** a `VsearchError` for the duration of a
session, so a consumer can catch it, skip the bad input and carry on
rather than have the whole host process die.

```cpp
struct VsearchError {   // declared in utils/fatal.hpp, included by vsearch_api.h
  std::string message;  // the same text fatal() writes to stderr
};
```

Throwing mode is enabled by `vsearch_session_begin()` and cleared by
`vsearch_session_end()` (so the `VsearchSession` guard turns it on for
its scope). It is a thread-local flag set only on the thread that opened
the session: worker threads spawned by the compute engines never throw
(they stop cooperatively and the main thread reports the error after
joining), because a C++ exception must not escape a `std::thread`.

Known error triggers include:

- Invalid FASTA/FASTQ format, or an unreadable/missing file, in
  `db.read()` / `udb_read()`
- FASTQ quality values outside `[opt_fastq_qmin, opt_fastq_qmax]`
- File I/O failures
- Out of memory (xmalloc/xrealloc failure)

**Catching and recovering:**

```cpp
#include "vsearch_api.h"

struct Parameters parameters;
try {
  VsearchSession const session(parameters);   // enables throwing mode
  Database db;
  db.read("maybe_bad.fasta", 0, parameters);   // any fatal() now unwinds
  // ... use the library ...
}                                              // ~VsearchSession ends the session
catch (VsearchError const & error) {
  std::fprintf(stderr, "vsearch: %s\n", error.message.c_str());
  // recover: the session was ended by the guard's destructor during the
  // unwind (files closed, database/index freed), so a fresh session and a
  // new run in the same process are safe.
}
```

Recovery relies on stack unwinding running destructors, so **own every
resource with RAII** (`VsearchSession`, `Database`, `Dbindex`, the
`OutputFileHandle` openers, the `*_session`/`*_info` handles via their
free functions or a `unique_ptr` wrapper). A `fatal()` caught this way
runs those destructors on the way out. If you drive the session with the
bare `vsearch_session_begin()` / `vsearch_session_end()` pair instead of
the guard, you must call `vsearch_session_end()` yourself in the catch
block before starting another session.

> **Resource cleanup on recovery.** A caught `VsearchError` unwinds
> cleanly: the engines' owning resources are RAII-managed (open files,
> the k-mer index, per-thread aligner/heap state, and each hit's
> alignment string, which `struct hit` now owns as a `std::string`), so
> they are released whether the fatal was a deterministic input error
> (bad or missing file, malformed record, out-of-range quality) or an
> out-of-memory / internal error deep in the search/cluster engines.
> Recovery leaves no leaks and never corrupts state — a subsequent run
> in the same process is safe. (This assumes your own objects are RAII
> too; see the try/catch example above.)

**Library API return codes** (independent of the throwing channel — these
are ordinary, non-fatal outcomes):

| Function | Success | Failure |
|----------|---------|---------|
| `chimera_detect_single()` | 0 | (throws `VsearchError` on internal error) |
| `MergePairs::merge()` | `result.merged == true` | `result.merged == false` (merge failed; strings empty) |
| All other functions | void | (throws `VsearchError` on error) |

**Recommendations for embedding:**

- Validate input data before passing to vsearch (DNA alphabet,
  non-empty sequences, reasonable lengths) to avoid the throw path
  entirely on the hot loop.
- Wrap each run in `try { VsearchSession ... } catch (VsearchError ...)`
  and keep all owned objects RAII, so a caught error cleans up.
- With `opt_quiet = true` (the library default), no messages are
  written to stderr during normal operation. A `fatal()` still writes
  its message to stderr before throwing (and to the `--log` file if one
  is open); the same text is available as `VsearchError::message`.

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

`MergeResult.sequence` / `quality` are `std::string`s owned by the result
(RAII) — no caller-side release. `search_result_s` no longer carries a
target header copy; look it up with `db.getheader(result.target)`.

### Database pointers

Pointers returned by `db.getheader()`, `db.getsequence()`, and
`db.getquality()` point into database-owned memory. They are valid
only while the database is loaded (until `db.clear()` is called or the
`Database` goes out of scope).

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
| Database loading (`db.init`, `db.add`, `dust_all`) | Single-threaded. |
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
   completes, the database and k-mer index are safe to read
   from any thread.

4. **Dereplication is self-contained.** `derep_session_s` has no global
   state dependencies and can be used concurrently across sessions
   (each with its own instance).

5. **Merging is stateless.** Once a `MergePairs` session is constructed,
   `merge()` is `const`; calls are fully independent and thread-safe.

6. **Masking:** `dust_single()` is thread-safe. `dust_all()` operates
   on the database and is single-threaded.

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
| `example_derep.cc` | Dereplication | Self-contained (no sequence DB), abundance output |
| `example_merge.cc` | Paired-end merging | FASTQ quality handling, stateless API |
| `example_dust.cc` | DUST masking | Standalone per-sequence, no init needed |
| `example_reinit.cc` | Re-initialization | Multi-session, multi-handle, gap penalty regression |
| `example_lifecycle.cc` | API memory/error contracts | Free null-safety, `MergeResult` RAII success/failure contracts, `MergePairs` session statelessness |
| `example_dbinfo.cc` | Database query/indexing surface | `db.read`, statistical accessors, quality, sort ordering, incremental `dbindex.add_sequence` |

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
    Database db;
    db.init();
    const char * h = "ref1";
    const char * s = "ACGTACGTACGTACGTACGT";
    db.add(false, h, s, nullptr, strlen(h), strlen(s), 1);
    dust_all(db, parameters);
    Dbindex dbindex;
    dbindex.prepare(1, opt_dbmask, db, parameters);
    dbindex.add_all_sequences(opt_dbmask, db, parameters);

    // Search
    struct search_session_s * ss = search_session_alloc();
    search_session_init(ss, parameters, dbindex, db);

    const char * qh = "query1";
    const char * qs = "ACGTACGTACGTACGTACGT";
    struct search_result_s results[5];
    int count;
    search_session_single(ss, qs, qh, strlen(qs), 1, results, 5, &count);

    for (int i = 0; i < count; i++) {
        printf("%s\t%s\t%.1f%%\n",
               qh, db.getheader(results[i].target), results[i].id);
    }

    // Cleanup
    search_session_cleanup(ss);
    search_session_free(ss);
    dbindex.clear();  // (also freed automatically at scope exit)
    db.clear();
    vsearch_session_end();
}
```
