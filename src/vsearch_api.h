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
 * by configuring a Parameters struct (the same options the CLI parses
 * from command-line flags) and calling the per-subsystem functions
 * declared here and in the module headers included below.
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
 *             * new fields appended to the end of the Parameters struct
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
 * Library users configure vsearch through a Parameters struct and open a
 * session with it before calling any other function. The required sequence is:
 *
 *   1. struct Parameters parameters; — a fresh, self-defaulting configuration
 *                                      (library defaults: quiet, no progress)
 *   2. parameters.opt_x = ...        — override any options for your use case
 *                                      (e.g., parameters.opt_wordlength = 8)
 *   3. VsearchSession session(parameters); — open the session: resolve the
 *                                      struct's sentinels/config and enable the
 *                                      recoverable-error mode on this thread
 *                                      (RAII; the session ends at scope exit)
 *   4. Database db;                  — an owned, empty database (RAII; passed
 *                                      by reference to steps 6-8)
 *   5. db.add() ...                  — load sequences (or db.read(file, ...));
 *      (db.init() first if           call db.init() before the first db.add()
 *       building with add())         when assembling a database programmatically
 *   6. dust_all(db, parameters)      — apply DUST masking
 *   7. Dbindex dbindex;              — build k-mer index (an owned object,
 *      dbindex.prepare(.., db, ..) +   passed by reference to step 8/9)
 *      dbindex.add_all_sequences(db)
 *   8. Per-subsystem session setup   — e.g., chimera_session_init(),
 *                                      search_session_init(ss, params, dbindex, db), etc.
 *   9. Per-thread working state      — e.g., chimera_info_alloc() +
 *                                      chimera_detect_thread_init()
 *                                      (repeat for each thread)
 *  10. Per-query calls               — e.g., chimera_detect_single(),
 *                                      search_session_single()
 *  11. Per-thread teardown + per-subsystem cleanup (reverse of 8-9)
 *  12. (everything released)         — db, dbindex and the session free
 *                                      themselves when they go out of scope (RAII)
 *
 * === Thread safety ===
 *
 * - Initialization (steps 1-8): single-threaded only
 * - Per-thread init (step 9): safe for different thread-local instances
 * - Computation (step 10): thread-safe IF each thread uses its own
 *   per-thread state. The database and k-mer index are
 *   read-only after step 7.
 * - Cleanup (step 11): single-threaded only
 *
 * There is no process-wide session lock: vsearch keeps no shared mutable
 * state, so independent sessions may run concurrently in different threads,
 * each owning its own Parameters, Database, Dbindex and VsearchSession. The
 * recoverable-error mode is per-thread, so construct the VsearchSession on the
 * thread that will catch VsearchError (worker threads spawned internally by a
 * compute engine keep the default, non-throwing behaviour).
 *
 * === Configuration ===
 *
 * All ~200 options live in the Parameters struct with correct library
 * defaults, so a default-constructed Parameters is already a valid, quiet
 * configuration; you only set the fields you care about. The struct is the
 * single configuration source: the compute engines read it directly (threaded
 * through the session/batch handles and the per-query state), so there are no
 * opt_* globals to set or keep in sync.
 *
 * === Re-initialization ===
 *
 * Multiple sequential sessions in the same process are supported: build a
 * fresh Parameters and open a new VsearchSession for each. Each fresh struct
 * re-applies the gap penalty adjustments from raw defaults.
 */

/* === API version === */

#define VSEARCH_API_VERSION_MAJOR 0
#define VSEARCH_API_VERSION_MINOR 17
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
   module headers, and defines the Parameters configuration struct. */
#include "vsearch.hpp"

/* Recoverable-error payload thrown by the engines when a fatal condition
   is hit inside a library session (see the error-handling section of
   LIBRARY_API.md). Catch VsearchError to skip a bad input and continue
   instead of the process being killed by std::exit(). */
#include "utils/fatal.hpp"

/* Module headers for API functions not declared in vsearch.h.
   searchcore.h is intentionally NOT included here — it defines
   internal types (searchinfo_s, struct hit, minwordmatches_defaults)
   that are implementation details of the search subsystem and not
   part of the public ABI. */
#include "core/chimera.hpp"
#include "core/cluster.hpp"
#include "core/db.hpp"
#include "core/dbindex.hpp"
#include "core/derep.hpp"
#include "core/mergepairs.hpp"
#include "core/mask.hpp"
#include "core/search.hpp"

/* === API version queries === */

/* Runtime API version, encoded identically to VSEARCH_API_VERSION.
   Lets callers loaded against a different libvsearch.a verify that the
   header they built against matches the runtime version. */
auto vsearch_api_version() -> int;

/* Runtime API version as "MAJOR.MINOR.PATCH". */
auto vsearch_api_version_string() -> const char *;

/* === Session lifecycle === */

/* A vsearch library session, owned by the caller. Construct one from a
   configured Parameters after setting any overrides on the struct and before
   any other vsearch function: the constructor resolves the struct's sentinel
   values via vsearch_apply_defaults_fixups() and puts fatal() into throwing
   mode on the calling thread, so a fatal condition unwinds as a catchable
   VsearchError instead of terminating the process. The destructor restores the
   thread's previous mode, so the session ends automatically at scope exit
   (including on an early return or while the stack unwinds). Non-copyable and
   non-movable.

   There is no process-wide lock and no begin/end pair to keep in sync: vsearch
   holds no shared mutable state, so independent sessions may run concurrently
   in different threads, each owning its own Parameters/Database/Dbindex/session.
   Construct the session on the thread that will catch VsearchError; nested
   sessions on one thread compose (the previous mode is saved and restored).

   Note: parameters.opt_minsize is NOT resolved here — it has command-specific
   defaults (1 for most commands, 8 for cluster_unoise). Set it explicitly if
   needed (the struct default is 0).

     {
       VsearchSession const session(parameters);
       // ... use the library ...
     }  // the session ends here
*/
class VsearchSession {
public:
  explicit VsearchSession(struct Parameters & parameters);
  ~VsearchSession();
  VsearchSession(VsearchSession const &) = delete;
  auto operator=(VsearchSession const &) -> VsearchSession & = delete;
private:
  bool previous_throw_mode;
};
