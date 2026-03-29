/*
 * example_derep.cc -- Full-length dereplication using the vsearch library API.
 *
 * Reads a FASTA file, deduplicates identical sequences, and outputs
 * unique sequences sorted by abundance (descending) with ;size=N annotation.
 *
 * Output matches vsearch --derep_fulllength --output --sizeout format.
 *
 * Build:  g++ -std=c++11 -O2 -I../src -o example_derep example_derep.cc ../src/libvsearch.a -lpthread -ldl
 * Run:    ./example_derep
 * Verify: diff output.fasta data/expected_derep.fasta
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
    vsearch_apply_defaults_fixups();

    /* 2. Read input sequences */
    std::vector<std::string> labels, seqs;
    read_fasta("data/derep_test.fasta", labels, seqs);

    /* 3. Create dereplication session and add all sequences */
    struct derep_session_s * ds = derep_session_alloc();
    derep_session_init(ds);

    for (size_t i = 0; i < labels.size(); i++) {
        derep_add_sequence(ds,
                           labels[i].c_str(),
                           seqs[i].c_str(),
                           static_cast<int>(seqs[i].size()),
                           1);  /* abundance = 1 per input sequence */
    }

    /* 4. Get results (sorted by abundance descending) */
    std::vector<struct derep_result_s> results(labels.size());
    int result_count = 0;
    derep_get_results(ds, results.data(),
                      static_cast<int>(results.size()), &result_count);

    /* 5. Output in FASTA format with ;size=N */
    for (int i = 0; i < result_count; i++) {
        auto const & r = results[i];
        std::printf(">%s;size=%lu\n%s\n",
                    r.header,
                    (unsigned long) r.abundance,
                    r.sequence);
    }

    /* 6. Cleanup */
    derep_session_cleanup(ds);
    derep_session_free(ds);

    return 0;
}
