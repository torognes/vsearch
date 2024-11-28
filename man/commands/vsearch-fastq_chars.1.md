% vsearch-fastq_chars(1) version 2.30.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastq_chars --- analyze fastq files to identify the
quality encoding and the range of quality score values used


# SYNOPSIS

| **vsearch** **\-\-fastq_chars** _fastqfile_ \[_options_]


# DESCRIPTION

The vsearch command `--fastq_chars` summarizes the number and
composition of sequences and quality strings contained in the input
fastq file. Results are written to the *standard error* `stderr(3)`,
and to *filename* if option `--log` *filename* is used.

The command `--fastq_chars` tries to automatically detect the quality
offset (33 or 64) and the fastq format (Solexa, Illumina 1.3+,
Illumina 1.5+ or Illumina 1.8+/Sanger) by analyzing the range of
observed quality score values. In case of success, `--fastq_chars`
suggests values for the quality offset `--fastq_ascii` (33 or 64), as
well as `--fastq_qmin` and `--fastq_qmax` values that could be used
with other commands. If the quality encoding is ambiguous, an offset
of 33 is favored. For example:

```text
Qmin 45, QMax 73, Range 29
Guess: -fastq_qmin 12 -fastq_qmax 40 -fastq_ascii 33
Guess: Original Sanger format (phred+33)
```

For each sequence symbol, `--fastq_chars` gives the number of
occurrences of the symbol, its relative frequency, and the length of
the longest run of that symbol (lowercase symbols are converted to
uppercase). For example:

```text
Letter          N   Freq MaxRun
------ ---------- ------ ------
     A    9050221  25.8%     19
     C    5657659  16.1%      9
     G    8731373  24.9%     11
     T   11610889  33.1%     24
```

For each quality symbol, `--fastq_chars` gives the ASCII value of the
symbol (see `ascii(7)`), its relative frequency, and the number of
times a *k*-mer of that symbol appears at the end of quality
strings. The length of the *k*-mer can be set with the option
`--fastq_tail` (4 by default). For example:

```text
Char  ASCII    Freq       Tails
----  -----  ------  ----------
 '-'     45    0.7%           8
 '.'     46    1.7%           9
 '/'     47    2.8%        2997
 '0'     48    2.4%         221
 '1'     49    0.9%           0
 '2'     50    0.4%           0
 '3'     51    0.4%           0
 '4'     52    0.2%           0
 '5'     53    0.3%           0
 '6'     54    0.0%           0
 '9'     57    1.3%          12
 ':'     58    0.7%           0
 ';'     59    1.3%           1
 '<'     60    0.6%           0
 '='     61    0.2%           0
 '>'     62    0.3%           0
 '?'     63    0.6%           0
 '@'     64    0.3%           0
 'A'     65    1.9%           0
 'B'     66    3.9%          91
 'C'     67    2.2%           0
 'D'     68    1.4%           0
 'E'     69    2.9%           0
 'F'     70   14.7%       24657
 'G'     71   23.5%        5890
 'H'     72   34.5%           9
 'I'     73    0.0%           0
```


# OPTIONS

## core options

`--fastq_tail` *positive non-null integer*
: Count the number of times a series of identical symbols of length
  *k* = *positive non-null integer*, a *k*-mer, appears at the end of
  quality strings. By default, *k* = 4.


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_log.md)
: A copy of the statistics computed by `--fastq_chars` is also recorded.

#(./fragments/option_no_progress.md)

#(./fragments/option_quiet.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Read from *input.fastq.gz*, count series of 5 identical symbols at the
end of quality strings (`--fastq_tail 5`), do not write to the
standard error (`--quiet`), and write results to *output.log*
(`--log`):

```sh
vsearch \
    --fastq_chars input.fastq.gz \
    --fastq_tail 5 \
    --quiet \
    --log output.log
```

# SEE ALSO

[`vsearch-fastq_stats(1)`](./commands/vsearch-stats.1.md)


#(./fragments/footer.md)
