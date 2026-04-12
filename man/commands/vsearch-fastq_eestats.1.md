% vsearch-fastq_eestats(1) version 2.30.6 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastq_eestats --- report per-position quality and expected error statistics


# SYNOPSIS

| **vsearch** **\-\-fastq_eestats** _fastqfile_ **\-\-output** _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--fastq_eestats` analyzes a fastq file and reports
per-position statistics on quality scores, per-base error probabilities, and
cumulative expected errors. The output is a 21-column tab-separated table
written to `--output`, with one row per read position. The columns are:

- **Pos**: position in the reads (1-based).
- **Reads**: number of reads extending to or past this position.
- **PctRecs**: percentage of reads reaching this position.
- **Q\_Min**, **Q\_Low**, **Q\_Med**, **Q\_Mean**, **Q\_Hi**, **Q\_Max**:
  distribution of the Phred quality score at this position (minimum, lower
  quartile, median, mean, upper quartile, maximum).
- **Pe\_Min**, **Pe\_Low**, **Pe\_Med**, **Pe\_Mean**, **Pe\_Hi**, **Pe\_Max**:
  distribution of the per-base error probability at this position.
- **EE\_Min**, **EE\_Low**, **EE\_Med**, **EE\_Mean**, **EE\_Hi**, **EE\_Max**:
  distribution of the cumulative expected errors accumulated from position 1
  to this position.

The quality encoding and accepted score range can be set with `--fastq_ascii`,
`--fastq_qmin`, and `--fastq_qmax`.

See [`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md) for a
description of the expected error metric.


# OPTIONS

## mandatory options

#(./fragments/option_output_eestats.md)


## core options

#(./fragments/option_fastq_ascii.md)

#(./fragments/option_fastq_qmax.md)

#(./fragments/option_fastq_qmin.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_quiet.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Analyze a fastq file and write per-position statistics to a file:

```sh
vsearch \
    --fastq_eestats input.fastq \
    --output eestats.tsv
```

Specify a Phred+64 quality encoding and write results quietly:

```sh
vsearch \
    --fastq_eestats input.fastq \
    --fastq_ascii 64 \
    --quiet \
    --output eestats.tsv
```


# SEE ALSO

[`vsearch-fastq_eestats2(1)`](./vsearch-fastq_eestats2.1.md),
[`vsearch-fastq_stats(1)`](./vsearch-fastq_stats.1.md),
[`vsearch-fastx_filter(1)`](./vsearch-fastx_filter.1.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md),
[`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md)


#(./fragments/footer.md)
