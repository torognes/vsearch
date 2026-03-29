/*
 * example_merge.cc — Paired-end read merging using the vsearch library API.
 *
 * Reads one forward and one reverse FASTQ read, merges them based on
 * overlap detection, and outputs the merged sequence in FASTA format.
 *
 * Build:  g++ -std=c++11 -O2 -I../src -o example_merge example_merge.cc ../src/libvsearch.a -lpthread -ldl
 * Run:    ./example_merge
 * Verify: diff output.fasta data/expected_merge.fasta
 */

#include "vsearch_api.h"

#include <cstdio>
#include <cstring>
#include <string>


/* Read one FASTQ record. Returns false on EOF/error. */
static bool read_fastq_record(const char * path,
                              std::string & header,
                              std::string & sequence,
                              std::string & quality) {
    std::FILE * fp = std::fopen(path, "r");
    if (fp == nullptr) { return false; }
    char line[65536];

    /* Line 1: @header */
    if (std::fgets(line, sizeof(line), fp) == nullptr) { std::fclose(fp); return false; }
    char * nl = std::strchr(line, '\n');
    if (nl != nullptr) { *nl = '\0'; }
    header = line + 1;  /* skip @ */

    /* Line 2: sequence */
    if (std::fgets(line, sizeof(line), fp) == nullptr) { std::fclose(fp); return false; }
    nl = std::strchr(line, '\n');
    if (nl != nullptr) { *nl = '\0'; }
    sequence = line;

    /* Line 3: + separator */
    if (std::fgets(line, sizeof(line), fp) == nullptr) { std::fclose(fp); return false; }

    /* Line 4: quality */
    if (std::fgets(line, sizeof(line), fp) == nullptr) { std::fclose(fp); return false; }
    nl = std::strchr(line, '\n');
    if (nl != nullptr) { *nl = '\0'; }
    quality = line;

    std::fclose(fp);
    return true;
}


int main() {
    /* 1. Initialize vsearch globals */
    vsearch_init_defaults();
    vsearch_apply_defaults_fixups();

    /* 2. Initialize merge quality lookup table */
    mergepairs_init();

    /* 3. Read paired-end reads */
    std::string fwd_header, fwd_seq, fwd_qual;
    std::string rev_header, rev_seq, rev_qual;

    if (!read_fastq_record("data/merge_fwd.fastq", fwd_header, fwd_seq, fwd_qual)) {
        std::fprintf(stderr, "Error reading forward read\n");
        return 1;
    }
    if (!read_fastq_record("data/merge_rev.fastq", rev_header, rev_seq, rev_qual)) {
        std::fprintf(stderr, "Error reading reverse read\n");
        return 1;
    }

    /* 4. Merge the pair */
    struct merge_result_s result;
    int rc = mergepairs_single(fwd_seq.c_str(), fwd_qual.c_str(),
                               static_cast<int>(fwd_seq.size()),
                               rev_seq.c_str(), rev_qual.c_str(),
                               static_cast<int>(rev_seq.size()),
                               fwd_header.c_str(), rev_header.c_str(),
                               &result);

    if (rc == 0 && result.merged) {
        /* Output in FASTA format matching vsearch --fastaout (80-char lines) */
        std::printf(">%s\n", fwd_header.c_str());
        int len = result.merged_length;
        const char * seq = result.merged_sequence;
        for (int i = 0; i < len; i += 80) {
            int chunk = (len - i < 80) ? len - i : 80;
            std::printf("%.*s\n", chunk, seq + i);
        }
    } else {
        std::fprintf(stderr, "Merge failed\n");
        return 1;
    }

    return 0;
}
