/*
 * example_reinit.cc — Regression test for library API re-initialization.
 *
 * Verifies that calling vsearch_init_defaults() + vsearch_apply_defaults_fixups()
 * multiple times in the same process produces correct gap penalties and
 * identical chimera detection results each time.
 *
 * This test catches the bug where a static bool in vsearch_apply_defaults_fixups()
 * permanently prevented gap-open penalty adjustment after the first call.
 *
 * Build:  g++ -std=c++11 -O2 -I../src -o example_reinit example_reinit.cc ../src/libvsearch.a -lpthread -ldl
 * Run:    ./example_reinit
 */

#include "vsearch_api.h"

#include <cassert>
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


struct session_results {
    int gap_open_query_interior;
    int gap_open_query_left;
    int gap_open_query_right;
    std::vector<char> flags;
    std::vector<double> scores;
};


static session_results run_session(
    const std::vector<std::string> & ref_labels,
    const std::vector<std::string> & ref_seqs,
    const std::vector<std::string> & query_labels,
    const std::vector<std::string> & query_seqs) {

    session_results results;

    /* Initialize */
    vsearch_init_defaults();
    opt_wordlength = 8;
    vsearch_apply_defaults_fixups();

    /* Record gap penalties after fixups */
    results.gap_open_query_interior = opt_gap_open_query_interior;
    results.gap_open_query_left = opt_gap_open_query_left;
    results.gap_open_query_right = opt_gap_open_query_right;

    /* Load reference database */
    db_init();
    for (size_t i = 0; i < ref_labels.size(); i++) {
        db_add(false, ref_labels[i].c_str(), ref_seqs[i].c_str(),
               nullptr, ref_labels[i].size(), ref_seqs[i].size(), 1);
    }
    dust_all();
    dbindex_prepare(1, opt_dbmask);
    dbindex_addallsequences(opt_dbmask);

    /* Detect chimeras */
    struct chimera_info_s * ci = chimera_info_alloc();
    chimera_detect_init(ci);

    for (size_t i = 0; i < query_labels.size(); i++) {
        struct chimera_result_s result;
        chimera_detect_single(ci,
                              query_seqs[i].c_str(),
                              query_labels[i].c_str(),
                              static_cast<int>(query_seqs[i].size()),
                              1,
                              &result);
        results.flags.push_back(result.flag);
        results.scores.push_back(result.score);
    }

    /* Cleanup */
    chimera_detect_cleanup(ci);
    chimera_info_free(ci);
    dbindex_free();
    db_free();

    return results;
}


int main() {
    /* Load test data once */
    std::vector<std::string> ref_labels, ref_seqs;
    read_fasta("data/chimera_ref.fasta", ref_labels, ref_seqs);

    std::vector<std::string> query_labels, query_seqs;
    read_fasta("data/chimera_queries.fasta", query_labels, query_seqs);

    if (ref_labels.empty()) {
        std::fprintf(stderr, "Error: could not load data/chimera_ref.fasta\n");
        return 1;
    }
    if (query_labels.empty()) {
        std::fprintf(stderr, "Error: could not load data/chimera_queries.fasta\n");
        return 1;
    }

    /* Run two identical sessions */
    std::fprintf(stderr, "Session 1...\n");
    session_results s1 = run_session(ref_labels, ref_seqs, query_labels, query_seqs);

    std::fprintf(stderr, "Session 2...\n");
    session_results s2 = run_session(ref_labels, ref_seqs, query_labels, query_seqs);

    /* Verify gap penalties are correct in both sessions */
    std::fprintf(stderr, "Checking gap penalties...\n");

    /* opt_gap_open_query_interior: raw=20, extension=2, expected=18 */
    assert(s1.gap_open_query_interior == 18);
    assert(s2.gap_open_query_interior == 18);

    /* opt_gap_open_query_left: raw=2, extension=1, expected=1 */
    assert(s1.gap_open_query_left == 1);
    assert(s2.gap_open_query_left == 1);

    /* opt_gap_open_query_right: raw=2, extension=1, expected=1 */
    assert(s1.gap_open_query_right == 1);
    assert(s2.gap_open_query_right == 1);

    /* Verify detection results match between sessions.
       Exact float equality is safe here: chimera_detect_single() uses only
       per-thread state (chimera_info_s) and the read-only global DB/index.
       The file-scope counters in chimera.cc (chimera_count, total_count, etc.)
       are CLI-only — never touched by the library API path. Both sessions
       execute identical computations on identical input with no shared
       mutable state, so scores are bit-for-bit reproducible. */
    std::fprintf(stderr, "Checking detection results...\n");
    assert(s1.flags.size() == s2.flags.size());
    for (size_t i = 0; i < s1.flags.size(); i++) {
        if (s1.flags[i] != s2.flags[i]) {
            std::fprintf(stderr,
                "FAIL: query %zu flag mismatch: session1=%c session2=%c\n",
                i, s1.flags[i], s2.flags[i]);
            return 1;
        }
        if (s1.scores[i] != s2.scores[i]) {
            std::fprintf(stderr,
                "FAIL: query %zu score mismatch: session1=%.4f session2=%.4f\n",
                i, s1.scores[i], s2.scores[i]);
            return 1;
        }
    }

    std::fprintf(stderr,
        "PASS: %zu queries, gap penalties correct, results identical across sessions\n",
        s1.flags.size());

    return 0;
}
