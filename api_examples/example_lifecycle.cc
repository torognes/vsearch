/*
 * example_lifecycle.cc — Library API contract tests.
 *
 * These checks exercise behaviour that is specific to *embedding* vsearch as a
 * library and that the CLI test-suite never reaches: the CLI runs a single
 * operation and exits, so it never frees a null handle, never inspects a
 * per-call result value, never reuses a merge session across pairs, and never
 * relies on the empty-on-failure contract of a result. Every check here is
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


/* Build a MergeInput (sequence + equal-length quality) from two strings. */
static MergeInput as_input(std::string const & sequence, std::string const & quality)
{
  return MergeInput{View<char>{sequence.data(), sequence.size()},
                    View<char>{quality.data(), quality.size()}};
}


/* --- Test 1: every free function is null-safe (documented "null-safe") ---
   A crash here (segfault) fails the test via a non-zero exit code. The merge
   API no longer has a free function: MergeResult owns its buffers (RAII). */
static int test_free_null_safety()
{
  chimera_info_free(nullptr);
  search_session_free(nullptr);
  cluster_session_free(nullptr);
  derep_session_free(nullptr);
  std::fprintf(stderr, "PASS: all free functions accept nullptr\n");
  return 0;
}


/* --- Test 2: MergePairs::merge success contract ---
   A known-overlapping pair merges; the returned MergeResult carries
   merged == true and sequence/quality strings whose length matches
   merged_length. The strings own their storage, freed when the result goes
   out of scope (no caller action). */
static int test_merge_success_contract(MergePairs const & merger,
                                       struct Parameters const & parameters)
{
  int failures = 0;

  /* A known-overlapping pair (same reads example_merge merges successfully). */
  std::string fwd_head, fwd_seq, fwd_qual;
  std::string rev_head, rev_seq, rev_qual;
  if ((not read_fastq_record("data/merge_fwd.fastq", fwd_head, fwd_seq, fwd_qual)) or
      (not read_fastq_record("data/merge_rev.fastq", rev_head, rev_seq, rev_qual)))
    {
      std::fprintf(stderr, "FAIL: could not read merge test reads\n");
      return failures + 1;
    }

  MergeResult const result = merger.merge(parameters,
                                          as_input(fwd_seq, fwd_qual),
                                          as_input(rev_seq, rev_qual));
  if (not result.merged)
    {
      std::fprintf(stderr, "FAIL: expected overlap pair did not merge\n");
      ++failures;
    }
  if ((result.sequence.size() != static_cast<std::size_t>(result.merged_length)) or
      (result.quality.size() != static_cast<std::size_t>(result.merged_length)))
    {
      std::fprintf(stderr, "FAIL: merged sequence/quality length disagrees with merged_length\n");
      ++failures;
    }
  if (result.merged and result.sequence.empty())
    {
      std::fprintf(stderr, "FAIL: successful merge produced an empty sequence\n");
      ++failures;
    }

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: a successful merge returns populated, self-owned buffers\n");
    }
  return failures;
}


/* --- Test 3: MergePairs::merge failure contract ---
   Two sequences that share no overlap must fail to merge: the documented
   contract is result.merged == false with both output strings left empty.
   This channel is unique to the library — the CLI simply counts a "not
   merged" read and moves on. */
static int test_merge_failure_contract(MergePairs const & merger,
                                       struct Parameters const & parameters)
{
  int failures = 0;

  std::string const fwd_seq(60, 'A');
  std::string const rev_seq(60, 'C');   /* rev-comp is 60x G — cannot overlap 60x A */
  std::string const qual(60, 'I');

  MergeResult const result = merger.merge(parameters,
                                          as_input(fwd_seq, qual),
                                          as_input(rev_seq, qual));
  if (result.merged)
    {
      std::fprintf(stderr, "FAIL: failed merge left result.merged == true\n");
      ++failures;
    }
  if ((not result.sequence.empty()) or (not result.quality.empty()))
    {
      std::fprintf(stderr, "FAIL: failed merge left non-empty output buffers\n");
      ++failures;
    }

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: non-overlapping merge returns merged == false with empty buffers\n");
    }
  return failures;
}


/* --- Test 4: session statelessness ---
   A single MergePairs session is documented as stateless and reusable across
   pairs (and shareable across threads). Verify two successive merges of the
   same pair on one session produce identical results. */
static int test_session_stateless(MergePairs const & merger,
                                  struct Parameters const & parameters)
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

  MergeResult const first = merger.merge(parameters,
                                         as_input(fwd_seq, fwd_qual),
                                         as_input(rev_seq, rev_qual));
  MergeResult const second = merger.merge(parameters,
                                          as_input(fwd_seq, fwd_qual),
                                          as_input(rev_seq, rev_qual));
  if ((not first.merged) or (not second.merged))
    {
      std::fprintf(stderr, "FAIL: repeated merge on one session did not merge\n");
      ++failures;
    }
  if ((first.sequence != second.sequence) or (first.quality != second.quality))
    {
      std::fprintf(stderr, "FAIL: repeated merge on one session gave different results\n");
      ++failures;
    }

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: a MergePairs session is stateless across successive merges\n");
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

  Database db;
  db.init();
  char const * const header = "ref1";
  char const * const sequence =
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
  db.add(false, header, sequence, nullptr,
         std::strlen(header), std::strlen(sequence), 1);
  dust_all(db, parameters);
  Dbindex dbindex;
  dbindex.prepare(1, parameters.opt_dbmask, db, parameters);
  dbindex.add_all_sequences(parameters.opt_dbmask, db, parameters);

  struct chimera_info_s * const info = chimera_info_alloc();
  chimera_detect_init(info, parameters, dbindex, db);

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
  dbindex.clear();
  db.clear();

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
  VsearchSession const session(parameters);
  MergePairs const merger(parameters);

  failures += test_free_null_safety();
  failures += test_merge_success_contract(merger, parameters);
  failures += test_merge_failure_contract(merger, parameters);
  failures += test_session_stateless(merger, parameters);
  failures += test_nonchimera_result_zeroed(parameters);
  failures += test_dust_hardmask();

  return failures == 0 ? 0 : 1;
}
