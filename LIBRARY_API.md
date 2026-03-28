# vsearch Library API

vsearch can be embedded as a static library (`libvsearch.a`) in other
applications. This document describes the library API, build process,
and usage patterns.

## Building

### Static library

```bash
autoreconf -fi
./configure
make -C src libvsearch.a
```

This produces `src/libvsearch.a`, a position-independent static archive
suitable for linking into shared libraries (e.g., loadable extensions).
All objects are compiled with `-fPIC`. The library excludes `main()`
via the `VSEARCH_NO_MAIN` preprocessor guard.

### Link dependencies

The library requires the following system libraries at link time:

- `-lpthread` (threading)
- `-ldl` (dynamic library loading for CPU feature detection)

Optional (detected at configure time):

- `-lz` (zlib, for gzip-compressed I/O)
- `-lbz2` (bzip2, for bzip2-compressed I/O)

### Include path

Add `src/` to your include path. The primary header for library
consumers is:

```c
#include "vsearch_api.h"
```

This includes `vsearch.h` (the aggregate header with all global
declarations) plus the module headers needed for the API functions.

## Initialization

vsearch uses approximately 200 global `opt_*` variables to control all
behavior. These **must** be initialized before calling any library
function. Missing or incorrect globals cause silent failures (e.g.,
wrong gap penalties cause all alignment candidates to be rejected).

### Required initialization sequence

```c
#include "vsearch_api.h"

// 1. Set all ~200 globals to CLI-equivalent defaults
vsearch_init_defaults();

// 2. Override globals for your use case
opt_wordlength = 8;          // k-mer word length
opt_minh = 0.28;             // chimera h-score threshold
// ... (see "Chimera detection options" below)

// 3. Resolve sentinel values (-1, 0) to computed defaults
//    Must be called AFTER setting overrides like opt_wordlength,
//    since opt_minwordmatches depends on it.
vsearch_apply_defaults_fixups();
```

### Global state warning

All vsearch state is global. There is no context object. This means:

- Only one "session" can be active at a time per process
- Initialization is not thread-safe
- Calling `vsearch_init_defaults()` resets ALL options, including any
  you may have set from a previous session

## Database management

vsearch uses a global sequence database. Sequences are loaded via
`db_init()` + `db_add()`, then indexed for k-mer search.

```c
// Reset database state (frees any previous data)
db_init();

// Load sequences one at a time
for (each sequence) {
    db_add(false,           // is_fastq (false for FASTA)
           header,          // null-terminated header string
           sequence,        // null-terminated sequence (ACGT, uppercase)
           nullptr,         // quality (nullptr for FASTA)
           header_length,   // strlen(header)
           sequence_length, // strlen(sequence)
           abundance);      // size annotation (1 for uchime_ref)
}

// Apply DUST low-complexity masking
dust_all();

// Build k-mer index
dbindex_prepare(1, opt_dbmask);   // 1 = use bitmap
dbindex_addallsequences(opt_dbmask);
```

### Cleanup

```c
dbindex_free();
db_free();
```

### Notes

- `db_init()` must be called before `db_add()` when not using `db_read()`.
  Without it, the `shortest` sequence length counter stays at 0, corrupting
  internal statistics.
- `db_init()` calls `db_free()` internally, so it is safe to call on an
  already-populated database.
- Sequences must be DNA (ACGT). RNA sequences (with U) must be converted
  to T before calling `db_add()`.

## Chimera detection

The chimera detection API provides per-query detection against a
pre-loaded reference database. It implements the UCHIME algorithm
(Edgar et al. 2011, Bioinformatics 27:2194-2200) with vsearch's
SIMD-optimized Needleman-Wunsch alignment.

### Per-thread working state

Each thread needs its own `chimera_info_s` instance. This is an opaque
type — allocate and free via the API functions:

```c
// Allocate opaque working state (one per thread)
struct chimera_info_s * ci = chimera_info_alloc();

// Initialize search infrastructure (SIMD aligners, k-mer finders, heaps)
// Must be called AFTER database is loaded and indexed.
chimera_detect_init(ci);

// Run detection for each query
struct chimera_result_s result;
chimera_detect_single(ci,
                      query_sequence,  // null-terminated, uppercase DNA
                      query_header,    // null-terminated
                      query_length,    // strlen(query_sequence)
                      query_abundance, // 1 for uchime_ref
                      &result);

// result.flag is 'Y' (chimera), 'N' (non-chimera), or '?' (borderline)
// result.score is the h-score
// result.parent_a_label, result.parent_b_label are populated for Y/?

// Cleanup (reverse order of init)
chimera_detect_cleanup(ci);
chimera_info_free(ci);
```

### Result struct

`chimera_result_s` contains all 18 fields from vsearch's `--uchimeout`
format:

| Field | Type | Description |
|-------|------|-------------|
| `score` | `double` | h-score |
| `query_label` | `char[1024]` | Query header (may truncate) |
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
| `flag` | `char` | Classification: 'Y', 'N', or '?' |

For non-chimeric results (`flag='N'`), only `query_label` and `flag`
are populated; all other fields are zero.

Labels may be silently truncated to 1023 characters.

### Chimera detection options

After calling `vsearch_init_defaults()`, override these for chimera
detection:

```c
// Required: set word length before calling vsearch_apply_defaults_fixups()
opt_wordlength = 8;

// UCHIME algorithm parameters (defaults shown)
opt_minh = 0.28;       // minimum h-score for chimera classification
opt_xn = 8.0;          // weight of no votes
opt_dn = 1.4;          // pseudo-count prior for votes
opt_mindiv = 0.8;      // minimum divergence
opt_mindiffs = 3;      // minimum differences from parents
opt_abskew = 2.0;      // abundance skew (de novo mode only)
```

Alignment scoring and masking options use vsearch defaults, which are
appropriate for chimera detection. Override only if you know what you
are doing.

### Thread safety

- **Database loading** (steps 1-7): single-threaded only.
- **Detection** (`chimera_detect_single`): thread-safe if each thread
  has its own `chimera_info_s`. The global database and k-mer index
  are read-only after indexing.
- **Cleanup**: single-threaded only.

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

## Complete example

```c
#include "vsearch_api.h"

void detect_chimeras(const char **headers, const char **sequences,
                     int n_refs, const char **query_headers,
                     const char **query_sequences, int n_queries) {
    // Initialize
    vsearch_init_defaults();
    opt_wordlength = 8;
    vsearch_apply_defaults_fixups();

    // Load reference database
    db_init();
    for (int i = 0; i < n_refs; i++) {
        db_add(false, headers[i], sequences[i], nullptr,
               strlen(headers[i]), strlen(sequences[i]), 1);
    }
    dust_all();
    dbindex_prepare(1, opt_dbmask);
    dbindex_addallsequences(opt_dbmask);

    // Set up per-thread state
    struct chimera_info_s *ci = chimera_info_alloc();
    chimera_detect_init(ci);

    // Detect chimeras
    for (int i = 0; i < n_queries; i++) {
        struct chimera_result_s result;
        chimera_detect_single(ci, query_sequences[i], query_headers[i],
                              strlen(query_sequences[i]), 1, &result);
        if (result.flag == 'Y') {
            printf("Chimera: %s (score=%.1f, parents=%s + %s)\n",
                   result.query_label, result.score,
                   result.parent_a_label, result.parent_b_label);
        }
    }

    // Cleanup
    chimera_detect_cleanup(ci);
    chimera_info_free(ci);
    dbindex_free();
    db_free();
}
```

## API reference

### Global initialization

| Function | Description |
|----------|-------------|
| `vsearch_init_defaults()` | Set all opt_* globals to CLI defaults. Call once before any other function. |
| `vsearch_apply_defaults_fixups()` | Resolve sentinel values. Call after setting overrides. |

### Database

| Function | Description |
|----------|-------------|
| `db_init()` | Reset database state. Call before `db_add()`. |
| `db_add(...)` | Add a sequence to the global database. |
| `db_free()` | Free database memory. |
| `db_getsequencecount()` | Return number of loaded sequences. |
| `db_getlongestsequence()` | Return length of longest sequence. |
| `dust_all()` | Apply DUST masking to all loaded sequences. |

### K-mer index

| Function | Description |
|----------|-------------|
| `dbindex_prepare(bitmap, mask)` | Allocate k-mer index structures. |
| `dbindex_addallsequences(mask)` | Index all loaded sequences. |
| `dbindex_addsequence(seqno, mask)` | Index a single sequence (de novo mode). |
| `dbindex_free()` | Free k-mer index memory. |

### Chimera detection

| Function | Description |
|----------|-------------|
| `chimera_info_alloc()` | Allocate opaque per-thread working state. |
| `chimera_info_free(ci)` | Free per-thread state (null-safe). |
| `chimera_detect_init(ci)` | Initialize search infrastructure. Call after DB is indexed. |
| `chimera_detect_single(ci, ...)` | Detect chimera for one query. Thread-safe with per-thread ci. |
| `chimera_detect_cleanup(ci)` | Free search infrastructure. Call before `chimera_info_free`. |

### Global search

| Function | Description |
|----------|-------------|
| `search_info_alloc()` | Allocate opaque per-thread search state. |
| `search_info_free(si)` | Free per-thread state (null-safe). |
| `search_init(si)` | Initialize SIMD aligners and k-mer finders. Call after DB is indexed. |
| `search_single(si, ...)` | Search one query against the database. Thread-safe with per-thread si. |
| `search_cleanup(si)` | Free search infrastructure. Call before `search_info_free`. |

Set `opt_id` to the minimum identity threshold (e.g., 0.97 for 97%)
and `opt_maxaccepts`/`opt_maxrejects` to control search depth before
calling `search_init()`.

### Paired-end merging

| Function | Description |
|----------|-------------|
| `mergepairs_init()` | Initialize quality score lookup table. Call once. |
| `mergepairs_single(...)` | Merge one forward/reverse read pair. Thread-safe. |

No per-thread state needed — each call is fully independent.
Override `opt_fastq_minovlen`, `opt_fastq_maxdiffs`, etc. for
merge parameters.

### Clustering

| Function | Description |
|----------|-------------|
| `cluster_session_alloc()` | Allocate opaque session state. |
| `cluster_session_free(cs)` | Free session state (null-safe). |
| `cluster_session_init(cs)` | Initialize for incremental clustering. Call after DB is loaded and `dbindex_prepare()` is called (but NOT `dbindex_addallsequences`). |
| `cluster_assign_single(cs, seqno, ...)` | Assign one sequence to a cluster. Must be called sequentially (seqno 0, 1, 2, ...). |
| `cluster_session_cleanup(cs)` | Free session resources. Call before `cluster_session_free`. |

Clustering is inherently sequential. Database must be pre-sorted by
length (cluster_fast) or abundance (cluster_size). Set `opt_id` to
the identity threshold. New centroids are indexed incrementally.
