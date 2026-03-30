/*
 * example_search.cc — Global sequence search using the vsearch library API.
 *
 * Loads a reference database, then searches each query sequence against it,
 * reporting the top hits above a minimum identity threshold.
 *
 * Output: tab-separated query, target, identity (matching vsearch
 * --usearch_global --userfields query+target+id).
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


int main() {
    /* 1. Initialize vsearch globals */
    vsearch_init_defaults();

    /* Search parameters */
    opt_wordlength = 8;
    opt_id = 0.5;              /* minimum 50% identity */
    opt_maxaccepts = 3;        /* return up to 3 hits per query */
    opt_maxrejects = 16;       /* give up after 16 rejected candidates */

    vsearch_apply_defaults_fixups();

    /* 2. Load reference database */
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

    /* 3. Initialize search state */
    struct searchinfo_s * si = search_info_alloc();
    search_init(si);

    /* 4. Search each query */
    std::vector<std::string> query_labels, query_seqs;
    read_fasta("data/chimera_queries.fasta", query_labels, query_seqs);

    for (size_t i = 0; i < query_labels.size(); i++) {
        struct search_result_s results[16];
        int count = 0;

        search_single(si,
                      query_seqs[i].c_str(),
                      query_labels[i].c_str(),
                      static_cast<int>(query_seqs[i].size()),
                      1,  /* abundance */
                      results,
                      3,  /* max results */
                      &count);

        for (int j = 0; j < count; j++) {
            std::printf("%s\t%s\t%.1f\n",
                        query_labels[i].c_str(),
                        results[j].target_label,
                        results[j].id);
        }
    }

    /* 5. Cleanup */
    search_cleanup(si);
    search_info_free(si);
    dbindex_free();
    db_free();
    vsearch_session_end();

    return 0;
}
