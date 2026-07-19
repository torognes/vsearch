/*
 * example_recover.cc — Library recoverable-error channel.
 *
 * This behaviour is unique to *embedding* vsearch as a library: while a
 * session is active, a fatal condition throws a VsearchError instead of
 * calling std::exit(), so the host process survives a bad input and can carry
 * on. The CLI never reaches this — it runs a single operation and, on a fatal
 * error, exits. Every check here is self-validating (no comparison against
 * native vsearch output): they verify the recover-and-continue contract from
 * the error-handling section of LIBRARY_API.md. The fatals are triggered at
 * increasing depth — a missing file (fails at open), a malformed file (fails
 * mid-read, with a partially-loaded database and an open handle to unwind),
 * and a compute-engine entry point (chimera_detect_single) — plus the session
 * contracts (message content, nested begin, and the explicit begin/end
 * recovery path).
 *
 * Build:  g++ -std=c++11 -O3 -I../src -o example_recover example_recover.cc ../src/libvsearch.a -lpthread -ldl
 * Run:    ./example_recover
 */

#include "vsearch_api.h"

#include <cstdio>
#include <cstring>
#include <string>


/* --- Test 1: a fatal inside a session is catchable, and the process
   continues ---
   A malformed input (here, a path that cannot be opened) makes db.read() hit
   fatal(), which — because a session is active — unwinds as a VsearchError
   rather than killing the process. After catching it, a second *good* read in
   the same process must still succeed: that is the recover-and-continue
   contract a library consumer relies on to skip a bad file and move on. */
static int test_recover_and_continue(struct Parameters const & parameters)
{
  int failures = 0;

  bool caught = false;
  std::string message;
  try
    {
      Database bad;
      bad.read("data/this_file_does_not_exist.fasta", 0, parameters);
      std::fprintf(stderr, "FAIL: reading a missing file did not raise VsearchError\n");
      ++failures;
    }
  catch (VsearchError const & error)
    {
      caught = true;
      message = error.message;   /* same text fatal() printed to stderr */
    }

  if (caught)
    {
      /* The payload must carry the same text fatal() printed — the only check
         that the throwing path populated VsearchError::message (an empty
         message would still be "caught" but useless to a consumer trying to
         report or classify the failure). */
      if (message.empty())
        {
          std::fprintf(stderr, "FAIL: caught VsearchError carried an empty message\n");
          ++failures;
        }
      else if (message.find("Unable to open file for reading") == std::string::npos)
        {
          std::fprintf(stderr, "FAIL: message did not describe the open failure: \"%s\"\n",
                       message.c_str());
          ++failures;
        }
      else
        {
          std::fprintf(stderr, "PASS: recovered from a fatal error: \"%s\"\n", message.c_str());
        }
    }

  /* The process is still alive, so a valid read in the same session works. */
  Database good;
  good.read("data/chimera_ref.fasta", 0, parameters);
  if (good.getsequencecount() != 6)
    {
      std::fprintf(stderr, "FAIL: post-recovery read got %lu sequences, expected 6\n",
                   (unsigned long) good.getsequencecount());
      ++failures;
    }
  else
    {
      std::fprintf(stderr, "PASS: a good read after recovery loaded all 6 sequences\n");
    }
  good.clear();

  return failures;
}


/* --- Test 2: recover from a fatal raised MID-READ, deep in db.read() ---
   test_recover_and_continue triggers the shallowest fatal: a missing file,
   which fails at open() before any read buffer or sequence storage exists. A
   *malformed* file is the deeper path the RAII unwind audit was written for:
   the fatal fires after the fastx handle is open and at least one record has
   been stored, so the unwind must free both the open handle and the
   partially-populated Database. data/malformed.fastq holds one good record
   followed by one whose quality line is shorter than its sequence. */
static int test_recover_from_midread(struct Parameters const & parameters)
{
  int failures = 0;

  bool caught = false;
  std::string message;
  try
    {
      Database bad;
      bad.read("data/malformed.fastq", 0, parameters);
      std::fprintf(stderr, "FAIL: reading a malformed FASTQ did not raise VsearchError\n");
      ++failures;
    }
  catch (VsearchError const & error)
    {
      caught = true;
      message = error.message;
    }

  if (caught)
    {
      if (message.find("Sequence and quality lines must be equally long") == std::string::npos)
        {
          std::fprintf(stderr, "FAIL: mid-read message was not the expected format error: \"%s\"\n",
                       message.c_str());
          ++failures;
        }
      else
        {
          std::fprintf(stderr, "PASS: recovered from a mid-read fatal: \"%s\"\n", message.c_str());
        }
    }

  /* Recover-and-continue: a good read in the same session still works. */
  Database good;
  good.read("data/chimera_ref.fasta", 0, parameters);
  if (good.getsequencecount() != 6)
    {
      std::fprintf(stderr, "FAIL: post-recovery read got %lu sequences, expected 6\n",
                   (unsigned long) good.getsequencecount());
      ++failures;
    }
  else
    {
      std::fprintf(stderr, "PASS: a good read after a mid-read fatal loaded all 6 sequences\n");
    }
  good.clear();

  return failures;
}


/* --- Test 3: recover from a fatal raised inside a compute-engine entry point
   ---
   The fatals above all originate in db.read(). This one fires inside
   chimera_detect_single(): its input validation (a query_len that disagrees
   with the query string) is a library-only contract — the CLI never calls this
   entry point. The throw unwinds while the per-thread chimera_info_s (SIMD
   aligners and k-mer finders allocated by chimera_detect_init) is live, so this
   checks that (a) an engine-internal fatal is catchable, (b) the per-thread
   handle survives it, and (c) a subsequent valid query on the SAME handle still
   works — the recover-and-continue contract at the engine level. Under the
   sanitizer CI job the caught unwind and the later cleanup must leak nothing. */
static int test_recover_from_engine_fatal(struct Parameters const & parameters)
{
  int failures = 0;

  Database db;
  db.init();
  char const * const header = "ref1";
  char const * const sequence =
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
  db.add(false, SeqRecord{View<char>{header, std::strlen(header)},
                          View<char>{sequence, std::strlen(sequence)},
                          View<char>{}}, 1);
  dust_all(db, parameters);
  Dbindex dbindex;
  dbindex.prepare(1, parameters.opt_dbmask, db, parameters);
  dbindex.add_all_sequences(parameters.opt_dbmask, db, parameters);

  struct chimera_info_s * const info = chimera_info_alloc();
  chimera_detect_init(info, parameters, dbindex, db);

  /* A query_len that does not match strlen(query_seq) is rejected by
     chimera_detect_single before any per-query work, raising a fatal that
     unwinds out of the engine. */
  bool caught = false;
  std::string message;
  try
    {
      struct chimera_result_s result;
      chimera_detect_single(info, sequence, "query1",
                            static_cast<int>(std::strlen(sequence)) + 5, 1, &result);
      std::fprintf(stderr, "FAIL: an inconsistent query length did not raise VsearchError\n");
      ++failures;
    }
  catch (VsearchError const & error)
    {
      caught = true;
      message = error.message;
    }

  if (caught && message.empty())
    {
      std::fprintf(stderr, "FAIL: engine fatal carried an empty message\n");
      ++failures;
    }
  else if (caught)
    {
      std::fprintf(stderr, "PASS: recovered from an engine-internal fatal: \"%s\"\n", message.c_str());
    }

  /* The handle must still be usable: a valid query on the same per-thread state
     succeeds (a corrupted handle would crash or misbehave here). */
  struct chimera_result_s result;
  int const status = chimera_detect_single(info, sequence, "query1",
                                           static_cast<int>(std::strlen(sequence)), 1, &result);
  if (status != 0)
    {
      std::fprintf(stderr, "FAIL: valid query after recovery returned %d (expected 0)\n", status);
      ++failures;
    }
  else
    {
      std::fprintf(stderr, "PASS: a valid query works on the same handle after a caught fatal\n");
    }

  chimera_detect_cleanup(info);
  chimera_info_free(info);
  dbindex.clear();
  db.clear();

  return failures;
}


/* --- Test 4: sessions stay reusable after a caught fatal ---
   A fatal that unwinds must run the RAII cleanup on its way out — in
   particular the VsearchSession guard, whose destructor calls
   vsearch_session_end() to clear throwing mode and release the session lock.
   Verify that by ending the first session (with the guard leaving scope during
   the unwind), then beginning a fresh one and doing real work: if the lock had
   been left held, vsearch_session_begin() below would itself fatally fail. */
static int test_session_reuse_after_recover()
{
  int failures = 0;

  /* First session: trigger and recover from a fatal, then let the guard end
     the session as it leaves this scope. */
  {
    struct Parameters first;
    VsearchSession const session(first);
    try
      {
        Database bad;
        bad.read("data/this_file_does_not_exist.fasta", 0, first);
      }
    catch (VsearchError const &)
      {
        /* recovered */
      }
  }

  /* Second, independent session in the same process must start cleanly. */
  struct Parameters second;
  VsearchSession const session(second);
  Database db;
  db.read("data/chimera_ref.fasta", 0, second);
  if (db.getsequencecount() != 6)
    {
      std::fprintf(stderr, "FAIL: fresh session after recovery got %lu sequences, expected 6\n",
                   (unsigned long) db.getsequencecount());
      ++failures;
    }
  else
    {
      std::fprintf(stderr, "PASS: a fresh session works after a caught fatal\n");
    }
  db.clear();

  return failures;
}


/* --- Test 5: a nested vsearch_session_begin() is a catchable fatal ---
   The session lock is taken with try_lock(), so a second begin while a session
   is active does not block: it raises a fatal diagnostic. On the
   session-owning thread throwing mode is already on, so that fatal surfaces as
   a catchable VsearchError rather than killing the process. The failed begin
   never acquired the lock, so the original session is intact and a single
   vsearch_session_end() releases it. */
static int test_double_begin_is_fatal()
{
  int failures = 0;

  struct Parameters first;
  vsearch_session_begin(first);

  bool caught = false;
  try
    {
      struct Parameters second;
      vsearch_session_begin(second);
      std::fprintf(stderr, "FAIL: a nested vsearch_session_begin() did not raise VsearchError\n");
      ++failures;
    }
  catch (VsearchError const &)
    {
      caught = true;
    }

  /* The first session still holds the (once-taken) lock; release it once. */
  vsearch_session_end();

  if (caught)
    {
      std::fprintf(stderr, "PASS: a nested vsearch_session_begin() is a catchable fatal\n");
    }
  else
    {
      ++failures;
    }

  /* A fresh session after the caught nested-begin must start cleanly, proving
     the lock was released exactly once (not double-unlocked). */
  struct Parameters again;
  vsearch_session_begin(again);
  vsearch_session_end();

  return failures;
}


/* --- Test 6: recovery via the explicit vsearch_session_begin/end pair ---
   Every recover test above uses the VsearchSession RAII guard, whose destructor
   ends the session during the unwind. LIBRARY_API.md documents an alternative:
   a caller driving the bare begin/end pair must call vsearch_session_end()
   itself in the catch block before starting another session. This exercises
   that manual path — a fatal is caught, the session is ended by hand, and a
   fresh explicit session then does real work. */
static int test_recover_with_explicit_session()
{
  int failures = 0;

  struct Parameters parameters;
  vsearch_session_begin(parameters);
  bool caught = false;
  try
    {
      Database bad;
      bad.read("data/this_file_does_not_exist.fasta", 0, parameters);
    }
  catch (VsearchError const &)
    {
      caught = true;
    }
  vsearch_session_end();   /* manual: the RAII guard is not in play here */

  if (not caught)
    {
      std::fprintf(stderr, "FAIL: explicit-session read of a missing file did not raise VsearchError\n");
      ++failures;
    }

  /* A fresh explicit session in the same process must work after the manual end
     above released the lock. */
  struct Parameters again;
  vsearch_session_begin(again);
  Database db;
  db.read("data/chimera_ref.fasta", 0, again);
  if (db.getsequencecount() != 6)
    {
      std::fprintf(stderr, "FAIL: explicit fresh session got %lu sequences, expected 6\n",
                   (unsigned long) db.getsequencecount());
      ++failures;
    }
  else
    {
      std::fprintf(stderr, "PASS: recovery via the explicit begin/end pair works\n");
    }
  db.clear();
  vsearch_session_end();

  return failures;
}


int main()
{
  int failures = 0;

  /* Tests 1-3 each run inside their own single session. */
  {
    struct Parameters parameters;
    VsearchSession const session(parameters);
    failures += test_recover_and_continue(parameters);
  }

  {
    struct Parameters parameters;
    VsearchSession const session(parameters);
    failures += test_recover_from_midread(parameters);
  }

  {
    struct Parameters parameters;
    parameters.opt_wordlength = 8;
    VsearchSession const session(parameters);
    failures += test_recover_from_engine_fatal(parameters);
  }

  /* Tests 4-6 manage their own (sequential) sessions. */
  failures += test_session_reuse_after_recover();
  failures += test_double_begin_is_fatal();
  failures += test_recover_with_explicit_session();

  return failures == 0 ? 0 : 1;
}
