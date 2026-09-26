/*
 * example_batch_session.cc — the batch functions reuse the session's threads
 * and per-thread state from one call to the next.
 *
 * search_batch, chimera_detect_batch and cluster_assign_batch share one pool
 * of worker threads, owned by the VsearchSession, and keep their per-thread
 * state between calls. Their results must not depend on it: this program
 * compares, against the single-query API,
 *   1. repeated and interleaved batch calls of several sizes in one session,
 *      with opt_threads changed between calls;
 *   2. batch calls after a fatal error was caught in the same session;
 *   3. a batch call from a thread with no session of its own (it creates its
 *      threads and state for the call);
 *   4. batch calls in a nested session, then again in the outer one;
 *   5. cluster_assign_batch in chunks, with opt_threads changed mid-run.
 *
 * Build:  g++ -std=c++11 -O3 -I../src -o example_batch_session example_batch_session.cc ../src/libvsearch.a -lpthread -ldl
 * Run:    ./example_batch_session   (exit status 0 when every check passes)
 */

#include "vsearch_api.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <thread>
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


static int const max_hits = 4;

struct Fixture {
    std::vector<std::string> ref_labels, ref_seqs, query_labels, query_seqs;
    std::vector<struct query_record_s> queries;
};


static void load_db(Database & db, Fixture const & fx, struct Parameters const & parameters) {
    db.init();
    for (size_t i = 0; i < fx.ref_labels.size(); i++) {
        db.add(false, SeqRecord{View<char>{fx.ref_labels[i].c_str(), fx.ref_labels[i].size()},
                                View<char>{fx.ref_seqs[i].c_str(), fx.ref_seqs[i].size()},
                                View<char>{}}, 1);
    }
    dust_all(db, parameters);
}


static bool same_chimera(struct chimera_result_s const & a, struct chimera_result_s const & b) {
    return a.flag == b.flag && a.score == b.score &&
        std::strcmp(a.query_label.data(), b.query_label.data()) == 0 &&
        std::strcmp(a.parent_a_label.data(), b.parent_a_label.data()) == 0 &&
        std::strcmp(a.parent_b_label.data(), b.parent_b_label.data()) == 0 &&
        a.id_query_model == b.id_query_model && a.divergence == b.divergence &&
        a.left_yes == b.left_yes && a.right_yes == b.right_yes;
}


static bool same_hit(struct search_result_s const & a, struct search_result_s const & b) {
    return a.target == b.target && a.id == b.id && a.matches == b.matches &&
        a.mismatches == b.mismatches && a.gaps == b.gaps &&
        a.alignment_length == b.alignment_length && a.strand == b.strand;
}


/* The references, from the single-query API */
struct Reference {
    std::vector<struct chimera_result_s> chimera;
    std::vector<struct search_result_s> hits;
    std::vector<int> counts;
};


static Reference single_query_reference(Fixture const & fx, struct Parameters const & parameters,
                                        Dbindex const & dbindex, Database const & db) {
    Reference ref;
    size_t const nq = fx.queries.size();
    ref.chimera.resize(nq);
    struct chimera_info_s * ci = chimera_info_alloc();
    chimera_detect_init(ci, parameters, dbindex, db);
    for (size_t i = 0; i < nq; i++) { chimera_detect_single(ci, fx.queries[i], &ref.chimera[i]); }
    chimera_detect_cleanup(ci);
    chimera_info_free(ci);

    ref.hits.resize(nq * max_hits);
    ref.counts.resize(nq);
    struct search_session_s * ss = search_session_alloc();
    search_session_init(ss, parameters, dbindex, db);
    for (size_t i = 0; i < nq; i++) {
        ref.counts[i] = search_session_single(ss, fx.queries[i],
                                              make_span(ref.hits).subspan(i * max_hits, max_hits));
    }
    search_session_cleanup(ss);
    search_session_free(ss);
    return ref;
}


/* Run both batch functions over all queries in batches of `batch`, alternating
   them call by call, and compare with the reference. */
static int batches_match(char const * what, Fixture const & fx, struct Parameters const & parameters,
                         Dbindex const & dbindex, Database const & db, Reference const & ref,
                         size_t const batch) {
    size_t const nq = fx.queries.size();
    std::vector<struct chimera_result_s> chimera(nq);
    std::vector<struct search_result_s> hits(nq * max_hits);
    std::vector<int> counts(nq, -1);
    for (size_t start = 0; start < nq; start += batch) {
        size_t const n = std::min(batch, nq - start);
        chimera_detect_batch(parameters, dbindex, db, make_view(fx.queries).subspan(start, n),
                             make_span(chimera).subspan(start, n));
        search_batch(parameters, dbindex, db, make_view(fx.queries).subspan(start, n),
                     make_span(hits).subspan(start * max_hits, n * max_hits), max_hits,
                     make_span(counts).subspan(start, n));
    }
    int failures = 0;
    for (size_t i = 0; i < nq; i++) {
        if (!same_chimera(chimera[i], ref.chimera[i])) {
            std::fprintf(stderr, "FAIL: %s: chimera query %zu differs (batch %zu, %lld threads)\n",
                         what, i, batch, static_cast<long long>(parameters.opt_threads));
            ++failures;
        }
        if (counts[i] != ref.counts[i]) {
            std::fprintf(stderr, "FAIL: %s: search query %zu: %d hits, expected %d\n",
                         what, i, counts[i], ref.counts[i]);
            ++failures;
            continue;
        }
        for (int j = 0; j < counts[i]; j++) {
            if (!same_hit(hits[i * max_hits + j], ref.hits[i * max_hits + j])) {
                std::fprintf(stderr, "FAIL: %s: search query %zu hit %d differs\n", what, i, j);
                ++failures;
            }
        }
    }
    return failures;
}


static int test_search_and_chimera(Fixture const & fx) {
    int failures = 0;
    struct Parameters parameters;
    parameters.opt_wordlength = 8;
    parameters.opt_id = 0.70;
    parameters.opt_maxaccepts = 1;
    parameters.opt_maxrejects = 32;
    parameters.opt_threads = 2;
    VsearchSession const session(parameters);

    Database db;
    load_db(db, fx, parameters);
    Dbindex dbindex;
    dbindex.prepare(parameters.opt_dbmask, db, parameters);
    dbindex.add_all_sequences(parameters.opt_dbmask, db, parameters);
    Reference const ref = single_query_reference(fx, parameters, dbindex, db);
    bool const has_chimera = std::any_of(ref.chimera.begin(), ref.chimera.end(),
                                         [](struct chimera_result_s const & r) { return r.flag == 'Y'; });
    bool const has_hits = std::any_of(ref.counts.begin(), ref.counts.end(),
                                      [](int const count) { return count > 0; });
    if (!has_chimera || !has_hits) {
        std::fprintf(stderr, "FAIL: the reference has no chimera or no hit (the test is vacuous)\n");
        ++failures;
    }

    /* 1. repeated calls, several batch sizes, thread count changed between calls */
    int const threads[] = {2, 3, 1, 4, 2};
    size_t const batches[] = {1, 3, fx.queries.size()};
    int before = failures;
    for (int const t : threads) {
        parameters.opt_threads = t;
        for (size_t const b : batches) {
            failures += batches_match("repeated calls", fx, parameters, dbindex, db, ref, b);
        }
    }
    if (failures == before) {
        std::fprintf(stderr, "PASS: repeated and interleaved batch calls match the single-query API\n");
    }

    /* 2. a fatal caught in the session leaves the pool and the states usable */
    before = failures;
    bool caught = false;
    try {
        Database bad;
        bad.read("data/this_file_does_not_exist.fasta", 0, parameters);
    } catch (VsearchError const &) {
        caught = true;
    }
    if (!caught) {
        std::fprintf(stderr, "FAIL: reading a missing file did not raise VsearchError\n");
        ++failures;
    }
    failures += batches_match("after a caught fatal", fx, parameters, dbindex, db, ref, 2);
    if (failures == before) {
        std::fprintf(stderr, "PASS: batch calls after a caught fatal match the single-query API\n");
    }

    /* 3. a thread with no session of its own */
    before = failures;
    int thread_failures = 0;
    std::thread other([&]() {
        thread_failures = batches_match("thread without a session", fx, parameters, dbindex, db, ref, 3);
    });
    other.join();
    failures += thread_failures;
    if (failures == before) {
        std::fprintf(stderr, "PASS: batch calls from a thread without a session match\n");
    }

    /* 4. a nested session, then the outer one again */
    before = failures;
    {
        /* a fresh struct: the session constructor resolves the gap penalties
           again, so a copy of an opened session's struct would not do */
        struct Parameters inner;
        inner.opt_wordlength = 8;
        inner.opt_id = 0.70;
        inner.opt_maxaccepts = 1;
        inner.opt_maxrejects = 32;
        VsearchSession const inner_session(inner);
        inner.opt_threads = 3;
        failures += batches_match("nested session", fx, inner, dbindex, db, ref, 2);
    }
    failures += batches_match("outer session after a nested one", fx, parameters, dbindex, db, ref, 2);
    if (failures == before) {
        std::fprintf(stderr, "PASS: batch calls in a nested session, then in the outer one, match\n");
    }

    dbindex.clear();
    db.clear();
    return failures;
}


static int test_cluster(Fixture const & fx) {
    int failures = 0;
    struct Parameters parameters;
    parameters.opt_wordlength = 8;
    parameters.opt_id = 0.70;
    parameters.opt_maxaccepts = 1;
    parameters.opt_maxrejects = 32;
    parameters.opt_threads = 2;
    VsearchSession const session(parameters);

    Database db;
    load_db(db, fx, parameters);
    db.sortbylength(parameters);
    int const sc = static_cast<int>(db.getsequencecount());

    auto const run = [&](int const chunk) -> std::vector<struct cluster_result_s> {
        std::vector<struct cluster_result_s> results(sc);
        Dbindex dbindex;
        dbindex.prepare(parameters.opt_qmask, db, parameters);
        struct cluster_session_s * cs = cluster_session_alloc();
        cluster_session_init(cs, parameters, dbindex, db);
        if (chunk == 0) {
            for (int i = 0; i < sc; i++) { cluster_assign_single(cs, i, &results[i]); }
        } else {
            for (int start = 0, call = 0; start < sc; start += chunk, ++call) {
                parameters.opt_threads = 1 + call % 3;  // 1, 2, 3, 1, ...: slots rebuilt
                int const n = std::min(chunk, sc - start);
                cluster_assign_batch(cs, start, make_span(results).subspan(start, n));
            }
        }
        cluster_session_cleanup(cs);
        cluster_session_free(cs);
        dbindex.clear();
        return results;
    };

    auto const single = run(0);
    for (int const chunk : {1, 2, sc}) {
        auto const batch = run(chunk);
        for (int i = 0; i < sc; i++) {
            auto const & a = single[i];
            auto const & b = batch[i];
            if (a.is_centroid != b.is_centroid || a.cluster_id != b.cluster_id ||
                a.centroid_seqno != b.centroid_seqno ||
                (!a.is_centroid && (std::fabs(a.identity - b.identity) > 0.0 ||
                                    std::strcmp(a.cigar.data(), b.cigar.data()) != 0))) {
                std::fprintf(stderr, "FAIL: cluster_assign_batch in chunks of %d: sequence %d differs\n",
                             chunk, i);
                ++failures;
            }
        }
    }
    if (failures == 0) {
        std::fprintf(stderr, "PASS: cluster_assign_batch in chunks, threads changed mid-run, matches\n");
    }
    db.clear();
    return failures;
}


int main() {
    Fixture fx;
    read_fasta("data/chimera_ref.fasta", fx.ref_labels, fx.ref_seqs);
    read_fasta("data/chimera_queries.fasta", fx.query_labels, fx.query_seqs);
    if (fx.ref_labels.empty() || fx.query_labels.empty()) {
        std::fprintf(stderr, "FAIL: cannot read the test data (run from api_examples/)\n");
        return 1;
    }
    for (size_t i = 0; i < fx.query_labels.size(); i++) {
        fx.queries.push_back(query_record_s{make_view(fx.query_labels[i]),
                                            make_view(fx.query_seqs[i]), 1});
    }

    int failures = 0;
    failures += test_search_and_chimera(fx);
    failures += test_cluster(fx);
    return failures == 0 ? 0 : 1;
}
