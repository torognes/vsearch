/*
 * example_recover.cc — Library recoverable-error channel.
 *
 * This behaviour is unique to *embedding* vsearch as a library: while a
 * session is active, a fatal condition throws a VsearchError instead of
 * calling std::exit(), so the host process survives a bad input and can carry
 * on. The CLI never reaches this — it runs a single operation and, on a fatal
 * error, exits. Both checks here are self-validating (no comparison against
 * native vsearch output): they verify the recover-and-continue contract from
 * the error-handling section of LIBRARY_API.md.
 *
 * Build:  g++ -std=c++11 -O3 -I../src -o example_recover example_recover.cc ../src/libvsearch.a -lpthread -ldl
 * Run:    ./example_recover
 */

#include "vsearch_api.h"

#include <cstdio>
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
      std::fprintf(stderr, "PASS: recovered from a fatal error: \"%s\"\n", message.c_str());
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


/* --- Test 2: sessions stay reusable after a caught fatal ---
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


int main()
{
  int failures = 0;

  /* Test 1 runs inside a single session created here. */
  {
    struct Parameters parameters;
    VsearchSession const session(parameters);
    failures += test_recover_and_continue(parameters);
  }

  /* Test 2 manages its own (sequential) sessions. */
  failures += test_session_reuse_after_recover();

  return failures == 0 ? 0 : 1;
}
