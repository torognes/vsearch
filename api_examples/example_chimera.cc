/*
 * example_chimera.cc — Chimera detection using the vsearch library API.
 *
 * Part 1: Loads reference sequences from a FASTA file, then classifies each
 * query sequence as chimeric (Y), non-chimeric (N), or borderline (?).
 * Output matches vsearch --uchime_ref --uchimeout format (18 tab-separated
 * columns). Run with no arguments; paths are hardcoded to data/ directory.
 *
 * Part 2: Self-validating batch vs sequential comparison.
 *
 * Build:  g++ -std=c++11 -O2 -I../src -o example_chimera example_chimera.cc ../src/libvsearch.a -lpthread -ldl
 * Run:    ./example_chimera
 * Verify: diff <(sort output.tsv) <(sort data/expected_chimera.tsv)
 */

#include "vsearch_api.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

/* Simple FASTA reader — reads all sequences into memory. */
static void read_fasta(const char * path,
                       std::vector<std::string> & labels,
                       std::vector<std::string> & sequences) {
    std::FILE * fp = std::fopen(path, "r");
    if (fp == nullptr) {
        std::fprintf(stderr, "Error: cannot open %s\n", path);
        return;
    }
    char line[65536];
    std::string label;
    std::string seq;
    while (std::fgets(line, sizeof(line), fp) != nullptr) {
        /* strip newline */
        char * nl = std::strchr(line, '\n');
        if (nl != nullptr) { *nl = '\0'; }
        nl = std::strchr(line, '\r');
        if (nl != nullptr) { *nl = '\0'; }
        if (line[0] == '>') {
            if (!label.empty()) {
                labels.push_back(label);
                sequences.push_back(seq);
            }
            label = line + 1;
            seq.clear();
        } else {
            seq += line;
        }
    }
    if (!label.empty()) {
        labels.push_back(label);
        sequences.push_back(seq);
    }
    std::fclose(fp);
}


/* --- Part 1: TSV output for diff comparison --- */
static int run_chimera_tsv() {
    vsearch_init_defaults();
    opt_wordlength = 8;
    vsearch_apply_defaults_fixups();

    std::vector<std::string> ref_labels, ref_seqs;
    read_fasta("data/chimera_ref.fasta", ref_labels, ref_seqs);

    db_init();
    for (size_t i = 0; i < ref_labels.size(); i++) {
        db_add(false, ref_labels[i].c_str(), ref_seqs[i].c_str(),
               nullptr, ref_labels[i].size(), ref_seqs[i].size(), 1);
    }
    dust_all();
    dbindex_prepare(1, opt_dbmask);
    dbindex_addallsequences(opt_dbmask);

    struct chimera_info_s * ci = chimera_info_alloc();
    chimera_detect_init(ci);

    std::vector<std::string> query_labels, query_seqs;
    read_fasta("data/chimera_queries.fasta", query_labels, query_seqs);

    for (size_t i = 0; i < query_labels.size(); i++) {
        struct chimera_result_s result;
        chimera_detect_single(ci,
                              query_seqs[i].c_str(),
                              query_labels[i].c_str(),
                              static_cast<int>(query_seqs[i].size()),
                              1,
                              &result);

        if (result.flag == 'N') {
            std::printf("%.4f\t%s\t*\t*\t*\t*\t*\t*\t*\t*\t"
                        "0\t0\t0\t0\t0\t0\t*\t%c\n",
                        result.score, result.query_label, result.flag);
        } else {
            std::printf("%.4f\t%s\t%s\t%s\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t"
                        "%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%c\n",
                        result.score,
                        result.query_label,
                        result.parent_a_label,
                        result.parent_b_label,
                        result.closest_parent_label,
                        result.id_query_model,
                        result.id_query_a,
                        result.id_query_b,
                        result.id_a_b,
                        result.id_query_top,
                        result.left_yes, result.left_no, result.left_abstain,
                        result.right_yes, result.right_no, result.right_abstain,
                        result.divergence,
                        result.flag);
        }
    }

    chimera_detect_cleanup(ci);
    chimera_info_free(ci);
    dbindex_free();
    db_free();
    vsearch_session_end();

    return 0;
}


/* --- Part 2: Self-validating batch vs sequential comparison --- */
static int run_batch_tests()
{
  int failures = 0;

  vsearch_init_defaults();
  opt_wordlength = 8;
  opt_threads = 2;
  vsearch_apply_defaults_fixups();

  std::vector<std::string> ref_labels, ref_seqs;
  read_fasta("data/chimera_ref.fasta", ref_labels, ref_seqs);

  db_init();
  for (size_t i = 0; i < ref_labels.size(); i++)
    {
      db_add(false, ref_labels[i].c_str(), ref_seqs[i].c_str(),
             nullptr, ref_labels[i].size(), ref_seqs[i].size(), 1);
    }
  dust_all();
  dbindex_prepare(1, opt_dbmask);
  dbindex_addallsequences(opt_dbmask);

  std::vector<std::string> query_labels, query_seqs;
  read_fasta("data/chimera_queries.fasta", query_labels, query_seqs);
  int const nq = static_cast<int>(query_labels.size());

  /* Sequential */
  std::vector<struct chimera_result_s> seq_results(nq);
  {
    struct chimera_info_s * ci = chimera_info_alloc();
    chimera_detect_init(ci);

    for (int i = 0; i < nq; i++)
      {
        chimera_detect_single(ci,
                              query_seqs[i].c_str(),
                              query_labels[i].c_str(),
                              static_cast<int>(query_seqs[i].size()),
                              1,
                              &seq_results[i]);
      }

    chimera_detect_cleanup(ci);
    chimera_info_free(ci);
  }

  /* Batch */
  std::vector<const char *> q_seqs(nq);
  std::vector<const char *> q_heads(nq);
  std::vector<int> q_lens(nq);
  std::vector<int> q_sizes(nq);
  for (int i = 0; i < nq; i++)
    {
      q_seqs[i] = query_seqs[i].c_str();
      q_heads[i] = query_labels[i].c_str();
      q_lens[i] = static_cast<int>(query_seqs[i].size());
      q_sizes[i] = 1;
    }

  std::vector<struct chimera_result_s> batch_results(nq);
  chimera_detect_batch(q_seqs.data(), q_heads.data(), q_lens.data(),
                       q_sizes.data(), nq, batch_results.data());

  /* Compare */
  for (int i = 0; i < nq; i++)
    {
      auto const & sr = seq_results[i];
      auto const & br = batch_results[i];

      bool mismatch = (sr.flag != br.flag) ||
                      (std::fabs(sr.score - br.score) > 0.001);
      if (sr.flag != 'N')
        {
          mismatch = mismatch ||
            (std::strcmp(sr.parent_a_label, br.parent_a_label) != 0) ||
            (std::strcmp(sr.parent_b_label, br.parent_b_label) != 0) ||
            (std::fabs(sr.id_query_a - br.id_query_a) > 0.1) ||
            (std::fabs(sr.id_query_b - br.id_query_b) > 0.1) ||
            (sr.left_yes != br.left_yes) ||
            (sr.right_yes != br.right_yes);
        }
      if (mismatch)
        {
          std::fprintf(stderr,
                       "FAIL: batch chimera query %d: flag=%c/score=%.4f "
                       "!= flag=%c/score=%.4f\n",
                       i, br.flag, br.score, sr.flag, sr.score);
          ++failures;
        }
    }

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: batch chimera matches sequential "
                   "(%d queries, %ld threads)\n", nq, (long) opt_threads);
    }

  dbindex_free();
  db_free();
  vsearch_session_end();

  return failures;
}


int main() {
    int rc = run_chimera_tsv();
    if (rc != 0)
      {
        return rc;
      }

    int failures = run_batch_tests();
    return failures == 0 ? 0 : 1;
}
