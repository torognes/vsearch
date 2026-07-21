/*
 * example_merge_error.cc — Distinguishing a hard input error from an ordinary
 * non-merge with the vsearch library API (MergeResult::error, since API 0.16.0).
 *
 * MergePairs::merge() returns merged == false both for an ordinary non-merge
 * (poor overlap, too many differences, ...) and for a hard input error: a FASTQ
 * quality symbol outside [fastq_qmin, fastq_qmax]. The CLI treats the latter as
 * fatal; the library reports it on the result via MergeResult::error /
 * error_value instead of throwing, so a caller can tell the two apart and decide
 * whether to skip the pair or stop its batch. This example asserts that contract.
 *
 * Build:  make example_merge_error   (or: g++ -std=c++11 -O3 -I../src \
 *                 -o example_merge_error example_merge_error.cc ../src/libvsearch.a -lpthread -ldl)
 * Run:    ./example_merge_error       (self-checking; exit 0 = all checks pass)
 */

#include "vsearch_api.h"

#include <cstdio>
#include <string>

namespace {

int checks_failed = 0;

void check(char const * label, bool ok) {
    std::printf("  %-38s %s\n", label, ok ? "PASS" : "FAIL");
    if (!ok) { ++checks_failed; }
}

/* Merge one forward/reverse pair given raw sequence + quality strings. */
MergeResult merge_pair(MergePairs const & merger,
                       struct Parameters const & parameters,
                       std::string const & fwd_seq, std::string const & fwd_qual,
                       std::string const & rev_seq, std::string const & rev_qual) {
    return merger.merge(
        parameters,
        MergeInput{View<char>{fwd_seq.data(), fwd_seq.size()},
                   View<char>{fwd_qual.data(), fwd_qual.size()}},
        MergeInput{View<char>{rev_seq.data(), rev_seq.size()},
                   View<char>{rev_qual.data(), rev_qual.size()}});
}

}  // namespace

int main() {
    /* One 40 bp read pair. The sequences do not overlap meaningfully, so a
       clean pair yields an ordinary non-merge (merged == false, error == none);
       what this example checks is the error channel, not the merge itself. */
    std::string const seq(40, 'A');
    std::string const qual_ok(40, 'I');   /* 'I' = ASCII 73 -> quality value 40, within [0, 41] */

    /* 1. Clean input: an ordinary non-merge must NOT be flagged as an error. */
    {
        struct Parameters parameters;
        VsearchSession const session(parameters);
        MergePairs const merger(parameters);

        MergeResult const result = merge_pair(merger, parameters, seq, qual_ok, seq, qual_ok);
        check("in-range quality -> error none",
              result.error == MergeError::none);
    }

    /* 2. A quality symbol above qmax: flagged as quality_above_qmax, and
          error_value carries the offending value. 'K' = ASCII 75 -> value 42,
          just above the default qmax of 41. */
    {
        struct Parameters parameters;
        VsearchSession const session(parameters);
        MergePairs const merger(parameters);

        std::string qual_bad = qual_ok;
        qual_bad[5] = 'K';
        MergeResult const result = merge_pair(merger, parameters, seq, qual_bad, seq, qual_ok);
        check("above qmax -> not merged",        !result.merged);
        check("above qmax -> error above_qmax",  result.error == MergeError::quality_above_qmax);
        check("above qmax -> error_value 42",    result.error_value == 42);
    }

    /* 3. A quality symbol below qmin: raise fastq_qmin so the in-range value 40
          now falls below it, and check it is flagged as quality_below_qmin. */
    {
        struct Parameters parameters;
        parameters.opt_fastq_qmin = 41;
        VsearchSession const session(parameters);
        MergePairs const merger(parameters);

        MergeResult const result = merge_pair(merger, parameters, seq, qual_ok, seq, qual_ok);
        check("below qmin -> error below_qmin",  result.error == MergeError::quality_below_qmin);
        check("below qmin -> error_value 40",    result.error_value == 40);
    }

    if (checks_failed == 0) {
        std::printf("All merge-error checks passed.\n");
        return 0;
    }
    std::fprintf(stderr, "%d merge-error check(s) failed.\n", checks_failed);
    return 1;
}
