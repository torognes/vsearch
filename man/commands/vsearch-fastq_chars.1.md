% vsearch-fastq_chars(1) version 2.30.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastq_chars --- analyze fastq files to identify the
quality encoding and the range of quality score values used.


# SYNOPSIS

| **vsearch** **\-\-fastq_chars** _fastqfile_ \[_options_]


# DESCRIPTION

The `vsearch --fastq_chars` command summarizes the number and
composition of sequences and quality strings contained in the input
fastq file. Results are written to the standard error, and to
*filename* if option `--log` *filename* is used.

For each sequence symbol, `--fastq_chars` gives the number of
occurrences of the symbol, its relative frequency and the length of
the longest run of that symbol (lowercase symbols are converted to
uppercase). For each quality symbol, `--fastq_chars` gives the ASCII
value of the symbol, its relative frequency, and the number of times a
*k*-mer of that symbol appears at the end of quality strings. The
length of the *k*-mer can be set with the option `--fastq_tail` (4 by
default).

The command `--fastq_chars` tries to automatically detect the quality
encoding (Solexa, Illumina 1.3+, Illumina 1.5+ or Illumina
1.8+/Sanger) by analyzing the range of observed quality score
values. In case of success, `--fastq_chars` suggests values for the
`--fastq_ascii` (33 or 64), as well as `--fastq_qmin` and
`--fastq_qmax` values to be used with the other commands that require
a fastq input file. If the quality encoding is ambiguous,
`--fastq_chars` will favor an offset of 33.


# OPTIONS

## core options

`--fastq_tail` *positive non-null integer*
: Count the number of times a series of characters of length *k* =
  *positive non-null integer*, a *k*-mer, appears at the end of
  quality strings. By default, *k* = 4.


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_log.md)
: A copy of the statistics computed by `--fastq_chars` is also recorded.

#(./fragments/option_no_progress.md)

#(./fragments/option_quiet.md)

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Read from *input.fastq*, count *5*-mers at the end of quality strings
(`--fastq_tail 5`), do not write to the standard error (`--quiet`),
and write results to *output.log* (`--log`):

```sh
vsearch \
    --fastq_chars input.fastq \
    --fastq_tail 5 \
    --quiet \
    --log output.log
```

#(./fragments/footer.md)
