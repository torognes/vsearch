/*
 * example_dust.cc -- DUST low-complexity masking using the vsearch library API.
 *
 * Reads a FASTA file and applies DUST soft masking (lowercase) to each
 * sequence independently. No database loading required.
 *
 * Output matches vsearch --fastx_mask --qmask dust --fastaout format.
 *
 * Build:  g++ -std=c++11 -O2 -I../src -o example_dust example_dust.cc ../src/libvsearch.a -lpthread -ldl
 * Run:    ./example_dust
 * Verify: diff output.fasta data/expected_dust.fasta
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
    /* 1. Initialize vsearch globals (needed for opt_hardmask) */
    vsearch_init_defaults();
    vsearch_apply_defaults_fixups();

    /* 2. Read sequences */
    std::vector<std::string> labels, seqs;
    read_fasta("data/dust_test.fasta", labels, seqs);

    /* 3. Apply DUST masking to each sequence and print */
    for (size_t i = 0; i < labels.size(); i++) {
        /* dust_single modifies in-place; work on a copy */
        std::vector<char> buf(seqs[i].begin(), seqs[i].end());
        buf.push_back('\0');

        dust_single(buf.data(), static_cast<int>(seqs[i].size()), false);

        std::printf(">%s\n%s\n", labels[i].c_str(), buf.data());
    }

    return 0;
}
