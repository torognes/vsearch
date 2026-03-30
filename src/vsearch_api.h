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
 * Public header for embedding vsearch as a static library.
 * Provides chimera detection and database management functions.
 *
 * === Initialization sequence ===
 *
 * Library users MUST initialize vsearch's global state before calling
 * any other function. The required sequence is:
 *
 *   1. vsearch_init_defaults()       — set all ~200 opt_* globals to
 *                                      safe CLI-equivalent defaults
 *   2. Override any opt_* globals    — e.g., opt_wordlength = 8 for
 *      needed for your use case        chimera detection
 *   3. vsearch_apply_defaults_fixups() — resolve sentinel values
 *                                        (opt_minwordmatches, opt_maxhits)
 *   4. db_init()                     — reset database state
 *   5. db_add() ...                  — load sequences
 *   6. dust_all()                    — apply DUST masking
 *   7. dbindex_prepare() +           — build k-mer index
 *      dbindex_addallsequences()
 *   8. chimera_info_alloc() +        — per-thread working state
 *      chimera_detect_init()
 *   9. chimera_detect_single()       — per-query detection
 *  10. chimera_detect_cleanup() +    — teardown
 *      chimera_info_free()
 *  11. dbindex_free() + db_free()    — release database
 *
 * === Thread safety ===
 *
 * - Steps 1-7 (initialization): single-threaded only
 * - Step 9 (detection): thread-safe IF each thread has its own
 *   chimera_info_s (from chimera_info_alloc + chimera_detect_init).
 *   The global database and k-mer index are read-only after step 7.
 * - Steps 10-11 (cleanup): single-threaded only
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
 * Repeat the full initialization sequence (steps 1-11) for each
 * session. vsearch_apply_defaults_fixups() correctly re-applies gap
 * penalty adjustments from the freshly-reset defaults each time.
 */

/* Aggregate header — includes config.h, system headers, and all
   module headers. Also declares all global opt_* extern variables. */
#include "vsearch.h"

/* Module headers for API functions not declared in vsearch.h */
#include "chimera.h"
#include "cluster.h"
#include "db.h"
#include "dbindex.h"
#include "derep.h"
#include "fastq_mergepairs.h"
#include "mask.h"
#include "search.h"
#include "searchcore.h"

/* === Global initialization === */

/* Set all opt_* globals to their CLI-equivalent default values.
   Must be called once before any other vsearch function.
   Allocates opt_ee_cutoffs_values via xmalloc. */
auto vsearch_init_defaults() -> void;

/* Resolve sentinel values in opt_* globals to computed defaults.
   Call after vsearch_init_defaults() and after setting any overrides.
   Handles: opt_maxhits (0 → INT64_MAX), opt_minwordmatches (-1 →
   wordlength-based default from searchcore.h minwordmatches_defaults).

   Note: opt_minsize is NOT resolved here — it has command-specific
   defaults (1 for most commands, 8 for cluster_unoise). Library users
   should set opt_minsize explicitly if needed (default from
   vsearch_init_defaults is 0). */
auto vsearch_apply_defaults_fixups() -> void;
