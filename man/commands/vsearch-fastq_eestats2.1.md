% vsearch-fastq_eestats2(1) version 2.30.4 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastq_eestats2 --- report read retention at combinations of length and expected error cutoffs


# SYNOPSIS

| **vsearch** **\-\-fastq_eestats2** _fastqfile_ **\-\-output** _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--fastq_eestats2` analyzes a fastq file and reports the
number (and percentage) of sequences that would be retained at each
combination of truncation length and maximum expected error threshold. The
output is useful for choosing values for `--fastq_trunclen` and
`--fastq_maxee` when running `--fastx_filter` (see
[`vsearch-fastx_filter(1)`](./vsearch-fastx_filter.1.md)).

The output is a tab-separated table written to `--output`. Each row
corresponds to a truncation length and each column (after the first) to an EE
threshold. Each cell shows the number of reads that would be retained and, in
parentheses, the corresponding percentage.

To illustrate with the default cutoffs (`--length_cutoffs 50,*,50` and
`--ee_cutoffs 0.5,1.0,2.0`):

```text
Length  0.500              1.000              2.000
    50  28595(100.0%)  28595(100.0%)  28595(100.0%)
   100   9854( 34.5%)  15240( 53.3%)  21567( 75.4%)
   150   3022( 10.6%)   6250( 21.9%)  11230( 39.3%)
```

The length cutoffs are set with `--length_cutoffs` (default `50,*,50`,
meaning 50, 100, 150, ..., up to the longest sequence in the file). The EE
cutoffs are set with `--ee_cutoffs` (default `0.5,1.0,2.0`).

The quality encoding and accepted score range can be set with `--fastq_ascii`,
`--fastq_qmin`, and `--fastq_qmax`.

See [`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md) for a
description of the expected error metric.


# OPTIONS

## mandatory options

#(./fragments/option_output_eestats.md)


## core options

#(./fragments/option_ee_cutoffs.md)

#(./fragments/option_fastq_ascii.md)

#(./fragments/option_fastq_qmax.md)

#(./fragments/option_fastq_qmin.md)

#(./fragments/option_length_cutoffs.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_quiet.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Analyze a fastq file with the default cutoffs and write the retention table:

```sh
vsearch \
    --fastq_eestats2 input.fastq \
    --output eestats2.tsv
```

Use custom length and EE cutoffs to focus on shorter reads:

```sh
vsearch \
    --fastq_eestats2 input.fastq \
    --length_cutoffs 100,300,50 \
    --ee_cutoffs 0.5,1.0,2.0,3.0 \
    --output eestats2.tsv
```


# SEE ALSO

[`vsearch-fastq_eestats(1)`](./vsearch-fastq_eestats.1.md),
[`vsearch-fastq_stats(1)`](./vsearch-fastq_stats.1.md),
[`vsearch-fastx_filter(1)`](./vsearch-fastx_filter.1.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md),
[`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md)


#(./fragments/footer.md)
