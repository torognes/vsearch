/*
 * example_cluster.cc — Sequence clustering using the vsearch library API.
 *
 * Part 1: Loads sequences from a FASTA file, sorts by length (descending),
 * and clusters them greedily at a given identity threshold. Output: UC format.
 *
 * Part 2: Self-validating batch vs sequential comparison.
 *
 * Build:  g++ -std=c++11 -O2 -I../src -o example_cluster example_cluster.cc ../src/libvsearch.a -lpthread -ldl
 * Run:    ./example_cluster
 * Verify: diff <(grep "^[SH]" output.uc | sort) <(grep "^[SH]" data/expected_cluster.uc | sort)
 */

#include "vsearch_api.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>


static void read_fasta(const char * path,
                       std::vector<std::string> & labels,
                       std::vector<std::string> & sequences) {
    std::FILE * fp = std::fopen(path, "r");
    if (fp == nullptr) { return; }
    char line[65536];
    std::string label, seq;
    while (std::fgets(line, sizeof(line), fp) != nullptr) {
        char * nl = std::strchr(line, '\n');
        if (nl != nullptr) { *nl = '\0'; }
        nl = std::strchr(line, '\r');
        if (nl != nullptr) { *nl = '\0'; }
        if (line[0] == '>') {
            if (!label.empty()) { labels.push_back(label); sequences.push_back(seq); }
            label = line + 1; seq.clear();
        } else {
            seq += line;
        }
    }
    if (!label.empty()) { labels.push_back(label); sequences.push_back(seq); }
    std::fclose(fp);
}


/* --- Part 1: UC output for diff comparison --- */
static int run_cluster_uc() {
    vsearch_init_defaults();
    opt_wordlength = 8;
    opt_id = 0.70;
    opt_maxaccepts = 1;
    opt_maxrejects = 32;
    vsearch_apply_defaults_fixups();

    std::vector<std::string> labels, seqs;
    read_fasta("data/chimera_ref.fasta", labels, seqs);

    db_init();
    for (size_t i = 0; i < labels.size(); i++) {
        db_add(false, labels[i].c_str(), seqs[i].c_str(),
               nullptr, labels[i].size(), seqs[i].size(), 1);
    }
    dust_all();
    db_sortbylength();
    dbindex_prepare(1, opt_qmask);

    struct cluster_session_s * cs = cluster_session_alloc();
    cluster_session_init(cs);

    int seqcount = static_cast<int>(db_getsequencecount());
    std::vector<struct cluster_result_s> results(seqcount);
    std::vector<int> cluster_sizes;

    for (int i = 0; i < seqcount; i++) {
        cluster_assign_single(cs, i, &results[i]);
        if (results[i].is_centroid) {
            cluster_sizes.push_back(1);
        } else {
            cluster_sizes[results[i].cluster_id]++;
        }
    }

    for (int i = 0; i < seqcount; i++) {
        auto const & r = results[i];
        if (r.is_centroid) {
            std::printf("S\t%d\t%lu\t*\t*\t*\t*\t*\t%s\t*\n",
                        r.cluster_id,
                        (unsigned long) db_getsequencelen(i),
                        db_getheader(i));
        } else {
            std::printf("H\t%d\t%lu\t%.1f\t+\t0\t0\t%s\t%s\t%s\n",
                        r.cluster_id,
                        (unsigned long) db_getsequencelen(i),
                        r.identity,
                        r.cigar[0] ? r.cigar : "*",
                        db_getheader(i),
                        r.centroid_label);
        }
    }

    for (int c = 0; c < static_cast<int>(cluster_sizes.size()); c++) {
        for (int i = 0; i < seqcount; i++) {
            if (results[i].is_centroid && results[i].cluster_id == c) {
                std::printf("C\t%d\t%d\t*\t*\t*\t*\t*\t%s\t*\n",
                            c, cluster_sizes[c], db_getheader(i));
                break;
            }
        }
    }

    cluster_session_cleanup(cs);
    cluster_session_free(cs);
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
  opt_id = 0.70;
  opt_maxaccepts = 1;
  opt_maxrejects = 32;
  opt_threads = 2;
  vsearch_apply_defaults_fixups();

  std::vector<std::string> labels, seqs;
  read_fasta("data/chimera_ref.fasta", labels, seqs);

  db_init();
  for (size_t i = 0; i < labels.size(); i++)
    {
      db_add(false, labels[i].c_str(), seqs[i].c_str(),
             nullptr, labels[i].size(), seqs[i].size(), 1);
    }
  dust_all();
  db_sortbylength();

  int const sc = static_cast<int>(db_getsequencecount());

  /* Sequential: use cluster_assign_single one at a time */
  dbindex_prepare(1, opt_qmask);
  struct cluster_session_s * cs_seq = cluster_session_alloc();
  cluster_session_init(cs_seq);

  std::vector<struct cluster_result_s> seq_results(sc);
  for (int i = 0; i < sc; i++)
    {
      cluster_assign_single(cs_seq, i, &seq_results[i]);
    }

  cluster_session_cleanup(cs_seq);
  cluster_session_free(cs_seq);
  dbindex_free();

  /* Batch: use cluster_assign_batch for all at once */
  dbindex_prepare(1, opt_qmask);
  struct cluster_session_s * cs_batch = cluster_session_alloc();
  cluster_session_init(cs_batch);

  std::vector<struct cluster_result_s> batch_results(sc);
  cluster_assign_batch(cs_batch, 0, sc, batch_results.data());

  cluster_session_cleanup(cs_batch);
  cluster_session_free(cs_batch);
  dbindex_free();

  /* Compare results */
  for (int i = 0; i < sc; i++)
    {
      auto const & sr = seq_results[i];
      auto const & br = batch_results[i];

      bool mismatch = (sr.is_centroid != br.is_centroid) ||
                      (sr.cluster_id != br.cluster_id) ||
                      (sr.centroid_seqno != br.centroid_seqno);
      if (!sr.is_centroid)
        {
          mismatch = mismatch ||
            (std::fabs(sr.identity - br.identity) > 0.01) ||
            (std::strcmp(sr.cigar, br.cigar) != 0);
        }
      if (mismatch)
        {
          std::fprintf(stderr,
                       "FAIL: batch cluster seq %d: "
                       "centroid=%d/cid=%d/cseq=%d/id=%.1f/cigar=%s "
                       "!= centroid=%d/cid=%d/cseq=%d/id=%.1f/cigar=%s\n",
                       i,
                       br.is_centroid, br.cluster_id,
                       br.centroid_seqno, br.identity, br.cigar,
                       sr.is_centroid, sr.cluster_id,
                       sr.centroid_seqno, sr.identity, sr.cigar);
          ++failures;
        }
    }

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: batch cluster matches sequential "
                   "(%d sequences, %ld threads)\n", sc, (long) opt_threads);
    }

  db_free();
  vsearch_session_end();

  return failures;
}


int main() {
    int rc = run_cluster_uc();
    if (rc != 0)
      {
        return rc;
      }

    int failures = run_batch_tests();
    return failures == 0 ? 0 : 1;
}
