/*
 * example_reinit.cc — Regression test for library API re-initialization.
 *
 * Verifies that beginning a session with vsearch_session_begin(Parameters &)
 * multiple times in the same process produces correct gap penalties and
 * identical chimera detection results each time.
 *
 * This test catches the bug where a static bool in the fixups step
 * permanently prevented gap-open penalty adjustment after the first call.
 *
 * Build:  g++ -std=c++11 -O3 -I../src -o example_reinit example_reinit.cc ../src/libvsearch.a -lpthread -ldl
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
    struct Parameters parameters;
    parameters.opt_wordlength = 8;
    vsearch_session_begin(parameters);

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
    vsearch_session_end();

    return results;
}


/* Adversarial global-state reset check (C1a regression net).

   Snapshot the defaults from a fresh init, garbage-fill a representative set of
   opt_* globals spanning every type (including the four that C1a fixed:
   opt_notmatchedfq, opt_bzip2_decompress, opt_clusterout_id/sort, and
   opt_strand, now a bool — the former C1b drift hazard), then begin a fresh
   vsearch_session_begin(Parameters{}) session and assert every one is reset. Generalizes the
   "init forgot to reset X" bug: a global added without a matching reset, or a
   dropped reset, is caught here as long as it is in one of the arrays below.
   The garbage values are chosen to differ from every default so a missing
   reset cannot coincidentally pass. Returns the number of failures. */
static int test_global_state_reset() {
    char ** const char_globals[] = {
        &opt_notmatchedfq, &opt_notmatched, &opt_matched,
        &opt_fastaout, &opt_alnout, &opt_userout, &opt_output };
    bool * const bool_globals[] = {
        &opt_clusterout_id, &opt_clusterout_sort, &opt_bzip2_decompress,
        &opt_gzip_decompress, &opt_sizein, &opt_sizeout, &opt_strand };
    int64_t * const i64_globals[] = {
        &opt_maxaccepts, &opt_wordlength, &opt_iddef, &opt_topn };
    double * const dbl_globals[] = {
        &opt_id, &opt_weak_id, &opt_xn, &opt_abskew };

    char * char_snap[sizeof(char_globals) / sizeof(char_globals[0])];
    bool bool_snap[sizeof(bool_globals) / sizeof(bool_globals[0])];
    int64_t i64_snap[sizeof(i64_globals) / sizeof(i64_globals[0])];
    double dbl_snap[sizeof(dbl_globals) / sizeof(dbl_globals[0])];
    const size_t nc = sizeof(char_globals) / sizeof(char_globals[0]);
    const size_t nb = sizeof(bool_globals) / sizeof(bool_globals[0]);
    const size_t ni = sizeof(i64_globals) / sizeof(i64_globals[0]);
    const size_t nd = sizeof(dbl_globals) / sizeof(dbl_globals[0]);

    /* 1. snapshot the defaults */
    struct Parameters snapshot_params;
    vsearch_session_begin(snapshot_params);
    for (size_t i = 0; i < nc; i++) { char_snap[i] = *char_globals[i]; }
    for (size_t i = 0; i < nb; i++) { bool_snap[i] = *bool_globals[i]; }
    for (size_t i = 0; i < ni; i++) { i64_snap[i] = *i64_globals[i]; }
    for (size_t i = 0; i < nd; i++) { dbl_snap[i] = *dbl_globals[i]; }
    vsearch_session_end();

    /* 2. garbage-fill (each value differs from the corresponding default) */
    static char junk[] = "garbage-not-a-default";
    for (size_t i = 0; i < nc; i++) { *char_globals[i] = junk; }
    for (size_t i = 0; i < nb; i++) { *bool_globals[i] = true; }
    for (size_t i = 0; i < ni; i++) { *i64_globals[i] = 0x5A5A5A5A; }
    for (size_t i = 0; i < nd; i++) { *dbl_globals[i] = -99999.0; }

    /* 3. re-init and assert every one is reset to its default */
    struct Parameters reinit_params;
    vsearch_session_begin(reinit_params);
    int failures = 0;
    for (size_t i = 0; i < nc; i++) {
        if (*char_globals[i] != char_snap[i]) {
            std::fprintf(stderr, "FAIL: char opt_* global #%zu not reset by init_defaults\n", i);
            ++failures;
        }
    }
    for (size_t i = 0; i < nb; i++) {
        if (*bool_globals[i] != bool_snap[i]) {
            std::fprintf(stderr, "FAIL: bool opt_* global #%zu not reset by init_defaults\n", i);
            ++failures;
        }
    }
    for (size_t i = 0; i < ni; i++) {
        if (*i64_globals[i] != i64_snap[i]) {
            std::fprintf(stderr, "FAIL: int64 opt_* global #%zu not reset by init_defaults\n", i);
            ++failures;
        }
    }
    for (size_t i = 0; i < nd; i++) {
        if (*dbl_globals[i] != dbl_snap[i]) {
            std::fprintf(stderr, "FAIL: double opt_* global #%zu not reset by init_defaults\n", i);
            ++failures;
        }
    }
    vsearch_session_end();
    return failures;
}


/* Version-consistency check: the accessors must agree with the compile-time
   macros. This is the one public entry point with no CLI equivalent to diff
   against, so it has no golden test. */
static int test_api_version() {
    int failures = 0;
    if (vsearch_api_version() != VSEARCH_API_VERSION) {
        std::fprintf(stderr, "FAIL: vsearch_api_version() != VSEARCH_API_VERSION\n");
        ++failures;
    }
    if (std::strcmp(vsearch_api_version_string(), VSEARCH_API_VERSION_STRING) != 0) {
        std::fprintf(stderr, "FAIL: vsearch_api_version_string() != VSEARCH_API_VERSION_STRING\n");
        ++failures;
    }
    return failures;
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

    /* === Test 2: Multi-handle via split API ===
       Create two chimera_info_s handles using the session/thread split.
       This was impossible before the fix — chimera_detect_init() would
       fatal on the second call due to re-initializing an active mutex. */
    std::fprintf(stderr, "Multi-handle test...\n");

    struct Parameters parameters;
    parameters.opt_wordlength = 8;
    vsearch_session_begin(parameters);

    db_init();
    for (size_t i = 0; i < ref_labels.size(); i++) {
        db_add(false, ref_labels[i].c_str(), ref_seqs[i].c_str(),
               nullptr, ref_labels[i].size(), ref_seqs[i].size(), 1);
    }
    dust_all();
    dbindex_prepare(1, opt_dbmask);
    dbindex_addallsequences(opt_dbmask);

    /* Session init once, then two per-thread handles */
    chimera_session_init();

    struct chimera_info_s * ci1 = chimera_info_alloc();
    chimera_detect_thread_init(ci1);

    struct chimera_info_s * ci2 = chimera_info_alloc();
    chimera_detect_thread_init(ci2);

    /* Run same queries through both handles and compare with session 1 */
    for (size_t i = 0; i < query_labels.size(); i++) {
        struct chimera_result_s r1, r2;
        chimera_detect_single(ci1,
                              query_seqs[i].c_str(),
                              query_labels[i].c_str(),
                              static_cast<int>(query_seqs[i].size()),
                              1, &r1);
        chimera_detect_single(ci2,
                              query_seqs[i].c_str(),
                              query_labels[i].c_str(),
                              static_cast<int>(query_seqs[i].size()),
                              1, &r2);

        /* Both handles must produce same result as session 1 */
        if (r1.flag != s1.flags[i] || r2.flag != s1.flags[i]) {
            std::fprintf(stderr,
                "FAIL: multi-handle query %zu flag mismatch: "
                "ci1=%c ci2=%c expected=%c\n",
                i, r1.flag, r2.flag, s1.flags[i]);
            return 1;
        }
        if (r1.score != s1.scores[i] || r2.score != s1.scores[i]) {
            std::fprintf(stderr,
                "FAIL: multi-handle query %zu score mismatch: "
                "ci1=%.4f ci2=%.4f expected=%.4f\n",
                i, r1.score, r2.score, s1.scores[i]);
            return 1;
        }
    }

    /* Per-thread cleanup, then session cleanup */
    chimera_detect_thread_cleanup(ci1);
    chimera_info_free(ci1);
    chimera_detect_thread_cleanup(ci2);
    chimera_info_free(ci2);
    chimera_session_cleanup();

    dbindex_free();
    db_free();
    vsearch_session_end();

    std::fprintf(stderr,
        "PASS: multi-handle detection, %zu queries identical across 2 handles\n",
        query_labels.size());

    /* === Test 3: apply_defaults_fixups() is idempotent (C1d) ===
       The gap-open penalties are adjusted in place by subtracting the
       gap-extension cost. Calling fixups twice without an intervening
       init must NOT double-subtract. */
    std::fprintf(stderr, "Idempotent-fixups test...\n");
    struct Parameters idem;
    idem.opt_wordlength = 8;
    vsearch_apply_defaults_fixups(idem);
    int const gap_after_first = idem.opt_gap_open_query_interior;   /* raw 20 - ext 2 = 18 */
    vsearch_apply_defaults_fixups(idem);                            /* again, no re-init */
    int const gap_after_second = idem.opt_gap_open_query_interior;
    if (gap_after_first != 18 || gap_after_second != 18) {
        std::fprintf(stderr,
            "FAIL: gap-open interior after 1st fixups=%d, after 2nd=%d (expected 18, 18)\n",
            gap_after_first, gap_after_second);
        return 1;
    }
    std::fprintf(stderr, "PASS: apply_defaults_fixups is idempotent (gap-open stays 18)\n");

    /* === Test 4: adversarial global-state reset (C1a) ===
       Garbage-fill a representative set of opt_* globals across every type,
       re-init, and assert all are reset. Supersedes the earlier hand-picked
       four-global check (which those globals are a subset of). */
    std::fprintf(stderr, "Global-reset test...\n");
    if (test_global_state_reset() != 0) {
        std::fprintf(stderr, "FAIL: init_defaults did not reset all sampled globals\n");
        return 1;
    }
    std::fprintf(stderr, "PASS: init_defaults resets sampled opt_* globals (all types)\n");

    /* === Test 5: API version consistency === */
    std::fprintf(stderr, "API-version test...\n");
    if (test_api_version() != 0) {
        return 1;
    }
    std::fprintf(stderr, "PASS: vsearch_api_version()/string match the macros\n");

    return 0;
}
