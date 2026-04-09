/*
 * example_search.cc — Global sequence search using the vsearch library API.
 *
 * Part 1: Loads a reference database, then searches each query sequence
 * against it, reporting the top hits above a minimum identity threshold.
 * Output: tab-separated query, target, identity (matching vsearch
 * --usearch_global --userfields query+target+id).
 *
 * Part 2: Self-validating strand test. Verifies that opt_strand = 2 finds
 * reverse-complement hits, and opt_strand = 1 does not.
 *
 * Build:  g++ -std=c++11 -O2 -I../src -o example_search example_search.cc ../src/libvsearch.a -lpthread -ldl
 * Run:    ./example_search
 * Verify: diff <(sort output.tsv) <(sort data/expected_search.tsv)
 */

#include "vsearch_api.h"

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


static std::string make_reverse_complement(const std::string & seq)
{
  std::string result;
  result.reserve(seq.size());
  for (int i = static_cast<int>(seq.size()) - 1; i >= 0; --i)
    {
      switch (seq[i])
        {
        case 'A': result += 'T'; break;
        case 'T': result += 'A'; break;
        case 'C': result += 'G'; break;
        case 'G': result += 'C'; break;
        default:  result += 'N'; break;
        }
    }
  return result;
}


/* --- Part 1: TSV output for diff comparison --- */
static int run_search_tsv() {
    vsearch_init_defaults();

    opt_wordlength = 8;
    opt_id = 0.5;
    opt_maxaccepts = 3;
    opt_maxrejects = 16;

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

    struct search_session_s * ss = search_session_alloc();
    search_session_init(ss);

    std::vector<std::string> query_labels, query_seqs;
    read_fasta("data/chimera_queries.fasta", query_labels, query_seqs);

    for (size_t i = 0; i < query_labels.size(); i++) {
        struct search_result_s results[16];
        int count = 0;

        search_single(ss,
                      query_seqs[i].c_str(),
                      query_labels[i].c_str(),
                      static_cast<int>(query_seqs[i].size()),
                      1,
                      results,
                      3,
                      &count);

        for (int j = 0; j < count; j++) {
            std::printf("%s\t%s\t%.1f\n",
                        query_labels[i].c_str(),
                        results[j].target_label,
                        results[j].id);
        }
    }

    search_session_cleanup(ss);
    search_session_free(ss);
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
  opt_id = 0.5;
  opt_maxaccepts = 3;
  opt_maxrejects = 16;
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

  /* Sequential: search each query one-by-one */
  std::vector<std::string> query_labels, query_seqs;
  read_fasta("data/chimera_queries.fasta", query_labels, query_seqs);
  int const nq = static_cast<int>(query_labels.size());

  int const max_per_query = 3;
  std::vector<struct search_result_s> seq_results(nq * max_per_query);
  std::vector<int> seq_counts(nq, 0);

  struct search_session_s * ss = search_session_alloc();
  search_session_init(ss);

  for (int i = 0; i < nq; i++)
    {
      search_single(ss,
                    query_seqs[i].c_str(),
                    query_labels[i].c_str(),
                    static_cast<int>(query_seqs[i].size()),
                    1,
                    &seq_results[i * max_per_query],
                    max_per_query,
                    &seq_counts[i]);
    }

  search_session_cleanup(ss);
  search_session_free(ss);

  /* Batch: search all queries at once */
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

  std::vector<struct search_result_s> batch_results(nq * max_per_query);
  std::vector<int> batch_counts(nq, 0);

  search_batch(q_seqs.data(), q_heads.data(), q_lens.data(), q_sizes.data(),
               nq, batch_results.data(), max_per_query, batch_counts.data());

  /* Compare */
  for (int i = 0; i < nq; i++)
    {
      if (seq_counts[i] != batch_counts[i])
        {
          std::fprintf(stderr,
                       "FAIL: batch query %d: count %d != sequential %d\n",
                       i, batch_counts[i], seq_counts[i]);
          ++failures;
          continue;
        }
      for (int j = 0; j < seq_counts[i]; j++)
        {
          auto const & sr = seq_results[i * max_per_query + j];
          auto const & br = batch_results[i * max_per_query + j];
          if (sr.target != br.target || sr.id != br.id ||
              sr.accepted != br.accepted || sr.strand != br.strand)
            {
              std::fprintf(stderr,
                           "FAIL: batch query %d hit %d: "
                           "target %d/%.1f/acc=%d/strand=%d "
                           "!= %d/%.1f/acc=%d/strand=%d\n",
                           i, j,
                           br.target, br.id, br.accepted, br.strand,
                           sr.target, sr.id, sr.accepted, sr.strand);
              ++failures;
            }
        }
    }

  if (failures == 0)
    {
      std::fprintf(stderr, "PASS: batch search matches sequential search "
                   "(%d queries, %ld threads)\n", nq, (long) opt_threads);
    }

  dbindex_free();
  db_free();
  vsearch_session_end();

  return failures;
}


/* --- Part 3: Self-validating strand tests --- */
static bool search_rc_finds_hit(const std::string & fwd,
                                const std::string & rev,
                                int strand_opt,
                                int * out_strand)
{
  vsearch_init_defaults();
  opt_wordlength = 8;
  opt_id = 0.97;
  opt_strand = strand_opt;
  opt_maxaccepts = 1;
  opt_maxrejects = 32;
  vsearch_apply_defaults_fixups();

  db_init();
  db_add(false, "fwd", fwd.c_str(), nullptr,
         3, static_cast<int>(fwd.size()), 1);
  dust_all();
  dbindex_prepare(1, opt_dbmask);
  dbindex_addallsequences(opt_dbmask);

  struct search_session_s * ss = search_session_alloc();
  search_session_init(ss);

  struct search_result_s results[4];
  int count = 0;

  search_single(ss,
                rev.c_str(),
                "rev_query",
                static_cast<int>(rev.size()),
                1,
                results,
                1,
                &count);

  bool found = (count > 0) && results[0].accepted;
  if (found && out_strand != nullptr)
    {
      *out_strand = results[0].strand;
    }

  search_session_cleanup(ss);
  search_session_free(ss);
  dbindex_free();
  db_free();
  vsearch_session_end();

  return found;
}


static int run_strand_tests()
{
  int failures = 0;

  std::string fwd =
    "AAGCCCAATCAACCACTCTCACTGGACGATTGCGGATATTGGCAACGAATTGGGAGGCGACCCGGACGACAG"
    "TCACGCCTTCTCGTTTGCGTACAGCTATTTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCCGCCG"
    "TGATACCGGACCAAACAAGACGTCCACTTCAATGTTTAAATGACCTACGCGTCAGAACACCTTTCTACTATG"
    "TGTTCTCCCAGAATCATCTAGTACAATGGCGCGTCGTCATTAAAGCACCGGATGCGACGAACGGAGCGTGAA"
    "TGAAGCTACTAC";
  std::string rev = make_reverse_complement(fwd);

  /* Test 1: plus-strand only — RC query should NOT find a hit */
  int strand_val = -1;
  bool strand1_found = search_rc_finds_hit(fwd, rev, 1, &strand_val);
  if (strand1_found)
    {
      std::fprintf(stderr, "FAIL: opt_strand=1: RC query unexpectedly matched fwd\n");
      ++failures;
    }
  else
    {
      std::fprintf(stderr, "PASS: opt_strand=1: RC query correctly found no hit\n");
    }

  /* Test 2: both strands — RC query MUST find a hit on minus strand */
  strand_val = -1;
  bool strand2_found = search_rc_finds_hit(fwd, rev, 2, &strand_val);
  if (!strand2_found)
    {
      std::fprintf(stderr, "FAIL: opt_strand=2: RC query did not match fwd\n");
      ++failures;
    }
  else if (strand_val != 1)
    {
      std::fprintf(stderr, "FAIL: opt_strand=2: hit strand=%d (expected 1)\n",
                  strand_val);
      ++failures;
    }
  else
    {
      std::fprintf(stderr, "PASS: opt_strand=2: RC query matched on minus strand\n");
    }

  return failures;
}


int main() {
    int rc = run_search_tsv();
    if (rc != 0)
      {
        return rc;
      }

    int failures = 0;
    failures += run_batch_tests();
    failures += run_strand_tests();
    return failures == 0 ? 0 : 1;
}
