% vsearch-fastq_convert(1) version 2.30.4 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastq_convert --- convert between fastq format variants


# SYNOPSIS

| **vsearch** **\-\-fastq_convert** _fastqfile_ **\-\-fastqout** _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--fastq_convert` converts a fastq file from one
quality-score encoding to another. The input encoding is specified with
`--fastq_ascii` (default: 33) and the output encoding with
`--fastq_asciiout` (default: 33). Both accept the values 33 (phred+33,
Sanger/Illumina 1.8+) and 64 (phred+64, Solexa/Illumina 1.3+/1.5+).

Quality scores are remapped during conversion. Output scores can be
clamped to a valid range with `--fastq_qminout` and `--fastq_qmaxout`.
The input score range is validated against `--fastq_qmin` and
`--fastq_qmax`.

Use `--fastq_chars` (see
[`vsearch-fastq_chars(1)`](./vsearch-fastq_chars.1.md)) to detect the
encoding of an unknown fastq file before converting.

To illustrate a conversion from phred+64 to phred+33:

```text
Input (phred+64):   --fastq_ascii 64 --fastq_asciiout 33:

>s1                 >s1
ACGTACGT            ACGTACGT
hijk                IJKL       (same Phred scores, different ASCII offset)
```


# OPTIONS

## mandatory options

#(./fragments/option_fastqout_convert.md)


## core options

#(./fragments/option_fastq_ascii.md)

#(./fragments/option_fastq_asciiout.md)

#(./fragments/option_fastq_qmaxout.md)

#(./fragments/option_fastq_qminout.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fastq_qmax.md)

#(./fragments/option_fastq_qmin.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_quiet.md)

#(./fragments/option_relabel.md)

#(./fragments/option_relabel_keep.md)

#(./fragments/option_relabel_md5.md)

#(./fragments/option_relabel_self.md)

#(./fragments/option_relabel_sha1.md)

#(./fragments/option_sample.md)

#(./fragments/option_sizein.md)

#(./fragments/option_sizeout.md)

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Convert a phred+64 file to the standard phred+33 encoding:

```sh
vsearch \
    --fastq_convert input.fastq \
    --fastq_ascii 64 \
    --fastq_asciiout 33 \
    --fastqout converted.fastq
```

Convert and clamp output quality scores to the range 0--40 (an older
Illumina convention):

```sh
vsearch \
    --fastq_convert input.fastq \
    --fastq_ascii 64 \
    --fastq_asciiout 33 \
    --fastq_qmaxout 40 \
    --fastqout converted.fastq
```


# SEE ALSO

[`vsearch-fastq_chars(1)`](./vsearch-fastq_chars.1.md),
[`vsearch-fastq_stats(1)`](./vsearch-fastq_stats.1.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
