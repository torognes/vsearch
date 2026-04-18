/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2026, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
  All rights reserved.

  This software is dual-licensed and available under a choice
  of one of two licenses, either under the terms of the GNU
  General Public License version 3 or the BSD 2-Clause License.

  See the file LICENSE for details.

*/

#pragma once

/*
 * vsearch library API
 *
 * Public header for embedding vsearch as a static library. Exposes the
 * same core operations as the CLI: sequence database management, DUST
 * masking, dereplication, paired-end merging, global search, greedy
 * clustering, and chimera detection. Callers drive these operations
 * by setting the opt_* globals that the CLI parses from command-line
 * flags and calling the per-subsystem functions declared here and in
 * the module headers included below.
 *
 * See LIBRARY_API.md for a full walk-through of each subsystem,
 * result structures, thread-safety rules, and worked examples.
 *
 * === API versioning ===
 *
 * The API follows semantic versioning. The version is exposed via
 * the VSEARCH_API_VERSION_* macros and the vsearch_api_version() and
 * vsearch_api_version_string() functions below.
 *
 *   MAJOR — incremented for incompatible changes:
 *             * removing, renaming, or retyping any public function
 *             * removing or reordering fields in a result struct
 *             * changing the semantics of an existing function
 *           Pre-1.0 (MAJOR == 0), the API is considered unstable
 *           and MINOR bumps may also break compatibility.
 *
 *   MINOR — incremented for backward-compatible additions:
 *             * new public functions
 *             * new fields appended to the end of a result struct
 *             * new opt_* globals (extern declarations)
 *
 *   PATCH — incremented for backward-compatible fixes that do not
 *           change the API surface (bug fixes, doc updates, internal
 *           refactors).
 *
 * Callers can compile-time check against VSEARCH_API_VERSION
 * (encoded as MAJOR*1000000 + MINOR*1000 + PATCH, matching the
 * OpenSSL / libcurl convention) or query the runtime version via
 * vsearch_api_version() to gate features that depend on a specific
 * version.
 *
 * === Initialization sequence ===
 *
 * Library users MUST initialize vsearch's global state before calling
 * any other function. The required sequence is:
 *
 *   1. vsearch_init_defaults()       — set all ~200 opt_* globals to
 *                                      safe library defaults (quiet,
 *                                      no progress output)
 *   2. Override any opt_* globals    — e.g., opt_wordlength = 8 for
 *      needed for your use case        chimera detection
 *   3. vsearch_apply_defaults_fixups() — resolve sentinel values
 *                                        (opt_minwordmatches, opt_maxhits)
 *   4. db_init()                     — reset database state
 *   5. db_add() ...                  — load sequences
 *   6. dust_all()                    — apply DUST masking
 *   7. dbindex_prepare() +           — build k-mer index
 *      dbindex_addallsequences()
 *   8. Per-subsystem session setup   — e.g., chimera_session_init(),
 *                                      search_session_init(), etc.
 *   9. Per-thread working state      — e.g., chimera_info_alloc() +
 *                                      chimera_detect_thread_init()
 *                                      (repeat for each thread)
 *  10. Per-query calls               — e.g., chimera_detect_single(),
 *                                      search_session_single()
 *  11. Per-thread teardown + per-subsystem cleanup (reverse of 8-9)
 *  12. dbindex_free() + db_free()    — release database
 *  13. vsearch_session_end()         — release session lock
 *
 * === Thread safety ===
 *
 * - Initialization (steps 1-8): single-threaded only
 * - Per-thread init (step 9): safe for different thread-local instances
 * - Computation (step 10): thread-safe IF each thread uses its own
 *   per-thread state. The global database and k-mer index are
 *   read-only after step 7.
 * - Cleanup (steps 11-13): single-threaded only
 *
 * vsearch_init_defaults() acquires a session mutex that blocks
 * concurrent callers until vsearch_session_end() releases it.
 * This prevents two threads from corrupting shared global state.
 * Forgetting to call vsearch_session_end() will cause the next
 * vsearch_init_defaults() call to block indefinitely.
 *
 * === Global state warning ===
 *
 * vsearch uses ~200 global opt_* variables. vsearch_init_defaults()
 * sets ALL of them. If you override any, do so BEFORE calling
 * vsearch_apply_defaults_fixups(). Missing or wrong globals cause
 * silent search failures (e.g., wrong gap penalties → rejected
 * candidates, wrong opt_minwordmatches → no k-mer matches).
 *
 * === Re-initialization ===
 *
 * Multiple sequential sessions in the same process are supported.
 * Repeat the full initialization sequence for each session.
 * vsearch_apply_defaults_fixups() correctly re-applies gap penalty
 * adjustments from the freshly-reset defaults each time.
 */

/* === API version === */

#define VSEARCH_API_VERSION_MAJOR 0
#define VSEARCH_API_VERSION_MINOR 1
#define VSEARCH_API_VERSION_PATCH 0

/* Encoded as MAJOR*1000000 + MINOR*1000 + PATCH (OpenSSL/libcurl
   convention), so each component may range 0..999 without collision:
   #if VSEARCH_API_VERSION >= 1002000  // >= 1.2.0 */
#define VSEARCH_API_VERSION \
  ((VSEARCH_API_VERSION_MAJOR * 1000000) + \
   (VSEARCH_API_VERSION_MINOR * 1000) + \
   VSEARCH_API_VERSION_PATCH)

/* Build "MAJOR.MINOR.PATCH" from the numeric macros at preprocess time
   so the string literal can never drift out of sync. */
#define VSEARCH_API_STRINGIFY_(x) #x
#define VSEARCH_API_STRINGIFY(x) VSEARCH_API_STRINGIFY_(x)
#define VSEARCH_API_VERSION_STRING                      \
  VSEARCH_API_STRINGIFY(VSEARCH_API_VERSION_MAJOR) "."  \
  VSEARCH_API_STRINGIFY(VSEARCH_API_VERSION_MINOR) "."  \
  VSEARCH_API_STRINGIFY(VSEARCH_API_VERSION_PATCH)

/* Aggregate header — includes config.h, system headers, and all
   module headers. Also declares all global opt_* extern variables. */
#include "vsearch.h"

/* Module headers for API functions not declared in vsearch.h.
   searchcore.h is intentionally NOT included here — it defines
   internal types (searchinfo_s, struct hit, minwordmatches_defaults)
   that are implementation details of the search subsystem and not
   part of the public ABI. */
#include "chimera.h"
#include "cluster.h"
#include "db.h"
#include "dbindex.h"
#include "derep.h"
#include "fastq_mergepairs.h"
#include "mask.h"
#include "search.h"

/* === API version queries === */

/* Runtime API version, encoded identically to VSEARCH_API_VERSION.
   Lets callers loaded against a different libvsearch.a verify that the
   header they built against matches the runtime version. */
auto vsearch_api_version() -> int;

/* Runtime API version as "MAJOR.MINOR.PATCH". */
auto vsearch_api_version_string() -> const char *;

/* === Global initialization === */

/* Set all opt_* globals to their CLI-equivalent default values.
   Must be called once before any other vsearch function.
   Allocates opt_ee_cutoffs_values via xmalloc.
   Acquires the session mutex — blocks if another session is active.
   Caller MUST call vsearch_session_end() when the session is done. */
auto vsearch_init_defaults() -> void;

/* Release the session mutex acquired by vsearch_init_defaults().
   Call after all cleanup (dbindex_free, db_free, etc.) is complete.
   Omitting this call will cause the next vsearch_init_defaults() to
   block indefinitely. */
auto vsearch_session_end() -> void;

/* Resolve sentinel values in opt_* globals to computed defaults.
   Call after vsearch_init_defaults() and after setting any overrides.
   Handles: opt_maxhits (0 → INT64_MAX), opt_minwordmatches (-1 →
   wordlength-based default from searchcore.h minwordmatches_defaults).

   Note: opt_minsize is NOT resolved here — it has command-specific
   defaults (1 for most commands, 8 for cluster_unoise). Library users
   should set opt_minsize explicitly if needed (default from
   vsearch_init_defaults is 0). */
auto vsearch_apply_defaults_fixups() -> void;
