/*
 * example_lifecycle.cc — Library API contract tests.
 *
 * These checks exercise behaviour that is specific to *embedding* vsearch as a
 * library and that the CLI test-suite never reaches: the CLI runs a single
 * operation and exits, so it never frees a null handle, never inspects a
 * return code, never reuses a result struct, and never relies on the
 * zero-initialisation contract of a result. Every check here is
 * self-validating (no comparison against native vsearch output): it verifies
 * the documented memory-ownership, null-safety and error-return contracts from
 * LIBRARY_API.md.
 *
 * Build:  g++ -std=c++11 -O3 -I../src -o example_lifecycle example_lifecycle.cc ../src/libvsearch.a -lpthread -ldl
 * Run:    ./example_lifecycle
 */

#include "vsearch_api.h"

#include <cstdio>
#include <cstring>
#include <string>


/* Read one FASTQ record (4 lines). Returns false on EOF/error. */
static bool read_fastq_record(char const * path,
                              std::string & header,
                              std::string & sequence,
                              std::string & quality)
{
  std::FILE * const input = std::fopen(path, "r");
  if (input == nullptr) { return false; }
  char line[65536];

  if (std::fgets(line, sizeof(line), input) == nullptr) { std::fclose(input); return false; }
  char * newline = std::strchr(line, '\n');
  if (newline != nullptr) { *newline = '\0'; }
  header = line + 1;   /* skip '@' */

  if (std::fgets(line, sizeof(line), input) == nullptr) { std::fclose(input); return false; }
  newline = std::strchr(line, '\n');
  if (newline != nullptr) { *newline = '\0'; }
  sequence = line;

  if (std::fgets(line, sizeof(line), input) == nullptr) { std::fclose(input); return false; }  /* '+' */

  if (std::fgets(line, sizeof(line), input) == nullptr) { std::fclose(input); return false; }
  newline = std::strchr(line, '\n');
  if (newline != nullptr) { *newline = '\0'; }
  quality = line;

  std::fclose(input);
  return true;
}


/* --- Test 1: every free function is null-safe (documented "null-safe") ---
   A crash here (segfault) fails the test via a non-zero exit code. */
static int test_free_null_safety()
{
  chimera_info_free(nullptr);
  search_session_free(nullptr);
  cluster_session_free(nullptr);
  derep_session_free(nullptr);
  merge_result_free(nullptr);
  std::fprintf(stderr, "PASS: all free functions accept nullptr\n");
  return 0;
}


/* --- Test 2: merge_result_free contract ---
   - a zero-initialised result is a no-op (nothing to free)
   - after freeing a populated result the pointers are nulled
   - a second free of the same result is safe (idempotent)         */
static int test_merge_result_free_contract(struct Parameters const & parameters)
{
  int failures = 0;

  /* zero-initialised: no buffers, free must be a harmless no-op */
  struct merge_result_s empty = {};
  merge_result_free(&empty);
  if ((empty.merged_sequence != nullptr) or (empty.merged_quality != nullptr))
    {
      std::fprintf(stderr, "FAIL: merge_result_free changed a zero-init result\n");
      ++failures;
    }

  /* A known-overlapping pair (same reads example_merge merges successfully). */
  std::string fwd_head, fwd_seq, fwd_qual;
  std::string rev_head, rev_seq, rev_qual;
  if ((not read_fastq_record("data/merge_fwd.fastq", fwd_head, fwd_seq, fwd_qual)) or
      (not read_fastq_record("data/merge_rev.fastq", rev_head, rev_seq, rev_qual)))
    {
      std::fprintf(stderr, "FAIL: could not read merge test reads\n");
      return failures + 1;
    }

  struct merge_result_s result = {};
  int const status = mergepairs_single(parameters,
                                       fwd_seq.c_str(), fwd_qual.c_str(),
                                       static_cast<int>(fwd_seq.size()),
                                       rev_seq.c_str(), rev_qual.c_str(),
                                       static_cast<int>(rev_seq.size()),
                                       fwd_head.c_str(), rev_head.c_str(), &result);
  if ((status != 0) or (not result.merged) or
      (result.merged_sequence == nullptr) or (result.merged_quality == nullptr))
    {
      std::fprintf(stderr, "FAIL: expected overlap pair did not merge (rc=%d)\n", status);
      ++failures;
    }

  merge_result_free(&result);
  if ((result.merged_sequence != nullptr) or (result.merged_quality != nullptr))
    {
      std::fprintf(stderr, "FAIL: merge_result_free did not null the buffers\n");
      ++failures;
    }

  /* double free must not crash and must leave the pointers null */
  merge_result_free(&result);
  if ((result.merged_sequence != nullptr) or (result.merged_quality != nullptr))
    {
      std::fprintf(stderr, "FAIL: second merge_result_free disturbed the null buffers\n");
      ++failures;
    }

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: merge_result_free is a no-op on empty, nulls on free, idempotent\n");
    }
  return failures;
}


/* --- Test 3: mergepairs_single failure return contract ---
   Two sequences that share no overlap must fail to merge: the documented
   contract is rc == -1, result.merged == false, and both output pointers
   left at nullptr so the caller can safely skip merge_result_free (though it
   remains a no-op). This return channel is unique to the library — the CLI
   simply counts a "not merged" read and moves on. */
static int test_merge_failure_contract(struct Parameters const & parameters)
{
  int failures = 0;

  std::string const fwd_seq(60, 'A');
  std::string const rev_seq(60, 'C');   /* rev-comp is 60x G — cannot overlap 60x A */
  std::string const qual(60, 'I');

  struct merge_result_s result = {};
  int const status = mergepairs_single(parameters,
                                       fwd_seq.c_str(), qual.c_str(),
                                       static_cast<int>(fwd_seq.size()),
                                       rev_seq.c_str(), qual.c_str(),
                                       static_cast<int>(rev_seq.size()),
                                       "fwd", "rev", &result);
  if (status != -1)
    {
      std::fprintf(stderr, "FAIL: non-overlapping merge returned %d (expected -1)\n", status);
      ++failures;
    }
  if (result.merged)
    {
      std::fprintf(stderr, "FAIL: failed merge left result.merged == true\n");
      ++failures;
    }
  if ((result.merged_sequence != nullptr) or (result.merged_quality != nullptr))
    {
      std::fprintf(stderr, "FAIL: failed merge left non-null output buffers\n");
      ++failures;
    }
  merge_result_free(&result);   /* must be safe even though nothing was allocated */

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: non-overlapping merge returns -1 with null buffers\n");
    }
  return failures;
}


/* --- Test 4: result-struct reuse contract ---
   The docs allow reusing one merge_result_s for successive merges provided
   merge_result_free() is called between calls. Verify a second merge into the
   same (freed) struct produces valid buffers again. */
static int test_result_reuse(struct Parameters const & parameters)
{
  int failures = 0;

  std::string fwd_head, fwd_seq, fwd_qual;
  std::string rev_head, rev_seq, rev_qual;
  if ((not read_fastq_record("data/merge_fwd.fastq", fwd_head, fwd_seq, fwd_qual)) or
      (not read_fastq_record("data/merge_rev.fastq", rev_head, rev_seq, rev_qual)))
    {
      std::fprintf(stderr, "FAIL: could not read merge test reads\n");
      return failures + 1;
    }

  struct merge_result_s result = {};
  for (int iteration = 0; iteration < 2; ++iteration)
    {
      int const status = mergepairs_single(parameters,
                                           fwd_seq.c_str(), fwd_qual.c_str(),
                                           static_cast<int>(fwd_seq.size()),
                                           rev_seq.c_str(), rev_qual.c_str(),
                                           static_cast<int>(rev_seq.size()),
                                           fwd_head.c_str(), rev_head.c_str(), &result);
      if ((status != 0) or (not result.merged) or (result.merged_sequence == nullptr))
        {
          std::fprintf(stderr, "FAIL: reuse iteration %d did not merge (rc=%d)\n",
                       iteration, status);
          ++failures;
        }
      merge_result_free(&result);
    }

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: a freed merge_result_s can be reused for another merge\n");
    }
  return failures;
}


/* --- Test 5: non-chimeric result zeroing contract ---
   LIBRARY_API.md: "For non-chimeric results (flag='N'), only query_label and
   flag are populated; all other fields are zero/empty." A library caller reads
   these struct fields directly, so the zeroing is a contract (the CLI only
   ever emits formatted text). Detect a clearly non-chimeric query against a
   tiny reference DB and assert the numeric fields are all zero. */
static int test_nonchimera_result_zeroed(struct Parameters const & parameters)
{
  int failures = 0;

  db_init();
  char const * const header = "ref1";
  char const * const sequence =
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
  db_add(false, header, sequence, nullptr,
         std::strlen(header), std::strlen(sequence), 1);
  dust_all(parameters);
  dbindex_prepare(1, parameters.opt_dbmask, parameters);
  dbindex_addallsequences(parameters.opt_dbmask, parameters);

  struct chimera_info_s * const info = chimera_info_alloc();
  chimera_detect_init(info, parameters);

  /* A query identical to the single reference cannot be a chimera. */
  struct chimera_result_s result;
  chimera_detect_single(info, sequence, "query1",
                        static_cast<int>(std::strlen(sequence)), 1, &result);

  if (result.flag != 'N')
    {
      std::fprintf(stderr, "FAIL: single-reference query classified as '%c', expected 'N'\n",
                   result.flag);
      ++failures;
    }
  bool const numbers_zeroed =
    (result.score == 0.0) and (result.id_query_model == 0.0) and
    (result.id_query_a == 0.0) and (result.id_query_b == 0.0) and
    (result.id_a_b == 0.0) and (result.id_query_top == 0.0) and
    (result.left_yes == 0) and (result.left_no == 0) and (result.left_abstain == 0) and
    (result.right_yes == 0) and (result.right_no == 0) and (result.right_abstain == 0) and
    (result.divergence == 0.0);
  if (not numbers_zeroed)
    {
      std::fprintf(stderr, "FAIL: non-chimeric result has non-zero numeric fields\n");
      ++failures;
    }
  if ((result.parent_a_label[0] != '\0') or (result.parent_b_label[0] != '\0'))
    {
      std::fprintf(stderr, "FAIL: non-chimeric result has non-empty parent labels\n");
      ++failures;
    }
  if (std::strcmp(result.query_label, "query1") != 0)
    {
      std::fprintf(stderr, "FAIL: non-chimeric result query_label = '%s', expected 'query1'\n",
                   result.query_label);
      ++failures;
    }

  chimera_detect_cleanup(info);
  chimera_info_free(info);
  dbindex_free();
  db_free();

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: non-chimeric result populates only query_label and flag\n");
    }
  return failures;
}


/* --- Test 6: dust_single honours the hardmask parameter ---
   The soft-mask path is covered by example_dust; this pins the other value of
   the boolean library parameter. Hard-masking replaces a low-complexity run
   with 'N'; soft-masking lowercases the same run. */
static int test_dust_hardmask()
{
  int failures = 0;

  std::string const original =
    "ACGTAGCTAGCTGATCGATCGATTTTTTTTTTTTTTTTTTTTTTTTTTTGCATGCATGCAT";

  std::string soft = original;
  dust_single(&soft[0], static_cast<int>(soft.size()), false);

  std::string hard = original;
  dust_single(&hard[0], static_cast<int>(hard.size()), true);

  bool const soft_lowered = (soft != original) and (soft.find('N') == std::string::npos);
  bool const hard_has_n = hard.find('N') != std::string::npos;
  if (not soft_lowered)
    {
      std::fprintf(stderr, "FAIL: soft dust_single did not lowercase a low-complexity run\n");
      ++failures;
    }
  if (not hard_has_n)
    {
      std::fprintf(stderr, "FAIL: hard dust_single did not insert 'N' in a low-complexity run\n");
      ++failures;
    }

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: dust_single hardmask=true inserts N, hardmask=false lowercases\n");
    }
  return failures;
}


int main()
{
  int failures = 0;

  struct Parameters parameters;
  parameters.opt_wordlength = 8;
  vsearch_session_begin(parameters);
  mergepairs_init(parameters);

  failures += test_free_null_safety();
  failures += test_merge_result_free_contract(parameters);
  failures += test_merge_failure_contract(parameters);
  failures += test_result_reuse(parameters);
  failures += test_nonchimera_result_zeroed(parameters);
  failures += test_dust_hardmask();

  vsearch_session_end();

  return failures == 0 ? 0 : 1;
}
