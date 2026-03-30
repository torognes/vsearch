/*
 * example_cluster.cc — Sequence clustering using the vsearch library API.
 *
 * Loads sequences from a FASTA file, sorts by length (descending), and
 * clusters them greedily at a given identity threshold. Each unmatched
 * sequence becomes a new centroid; matched sequences join the best
 * centroid's cluster.
 *
 * Output: UC format (S/H/C records) matching vsearch --cluster_fast --uc.
 *
 * Build:  g++ -std=c++11 -O2 -I../src -o example_cluster example_cluster.cc ../src/libvsearch.a -lpthread -ldl
 * Run:    ./example_cluster
 * Verify: diff <(grep "^[SH]" output.uc | sort) <(grep "^[SH]" data/expected_cluster.uc | sort)
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
    opt_wordlength = 8;
    opt_id = 0.70;             /* 70% identity threshold */
    opt_maxaccepts = 1;        /* accept first hit above threshold */
    opt_maxrejects = 32;       /* search depth */
    vsearch_apply_defaults_fixups();

    /* 2. Load sequences — vsearch --cluster_fast sorts by length */
    std::vector<std::string> labels, seqs;
    read_fasta("data/chimera_ref.fasta", labels, seqs);

    db_init();
    for (size_t i = 0; i < labels.size(); i++) {
        db_add(false, labels[i].c_str(), seqs[i].c_str(),
               nullptr, labels[i].size(), seqs[i].size(), 1);
    }
    dust_all();

    /* Sort by length (descending) — matches --cluster_fast behavior */
    db_sortbylength();

    /* Prepare k-mer index but do NOT index all sequences.
       Centroids are indexed incrementally during clustering. */
    dbindex_prepare(1, opt_qmask);

    /* 3. Initialize clustering session */
    struct cluster_session_s * cs = cluster_session_alloc();
    cluster_session_init(cs);

    /* 4. Cluster all sequences sequentially */
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

    /* 5. Output in UC format (S/H records, then C records) */
    for (int i = 0; i < seqcount; i++) {
        auto const & r = results[i];
        if (r.is_centroid) {
            /* S record: cluster seed */
            std::printf("S\t%d\t%lu\t*\t*\t*\t*\t*\t%s\t*\n",
                        r.cluster_id,
                        (unsigned long) db_getsequencelen(i),
                        db_getheader(i));
        } else {
            /* H record: hit assigned to cluster */
            std::printf("H\t%d\t%lu\t%.1f\t+\t0\t0\t%s\t%s\t%s\n",
                        r.cluster_id,
                        (unsigned long) db_getsequencelen(i),
                        r.identity,
                        r.cigar[0] ? r.cigar : "*",
                        db_getheader(i),
                        r.centroid_label);
        }
    }

    /* C records: cluster summaries */
    for (int c = 0; c < static_cast<int>(cluster_sizes.size()); c++) {
        /* Find centroid label for this cluster */
        for (int i = 0; i < seqcount; i++) {
            if (results[i].is_centroid && results[i].cluster_id == c) {
                std::printf("C\t%d\t%d\t*\t*\t*\t*\t*\t%s\t*\n",
                            c, cluster_sizes[c], db_getheader(i));
                break;
            }
        }
    }

    /* 6. Cleanup */
    cluster_session_cleanup(cs);
    cluster_session_free(cs);
    dbindex_free();
    db_free();
    vsearch_session_end();

    return 0;
}
