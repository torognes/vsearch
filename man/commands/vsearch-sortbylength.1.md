% vsearch-sortbylength(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-sortbylength --- sort fasta or fastq sequences by decreasing
length


# SYNOPSIS

| **vsearch** **\-\-sortbylength** _fastxfile_ \-\-output _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--sortbylength` sorts fasta or fastq sequences by
decreasing length and writes them to the file specified with `--output`,
in fasta format. To illustrate:

```text
>s1                >s3
AAAA               TTTTTTTT
>s2  --- sort -->  >s2
CCCCCC             CCCCCC
>s3                >s1
TTTTTTTT           AAAA
```

When two sequences have the same length, ties are broken first by
decreasing abundance (if present) and then by label in alphanumerical
order. Sequences can be filtered by length with `--minseqlength` and
`--maxseqlength`. The output can be limited to the first *n* sequences
with `--topn`.

See [`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md) and
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md) for more
information on input formats.


# OPTIONS

## mandatory options

#(./fragments/option_output.md)


## core options

#(./fragments/option_maxseqlength.md)

#(./fragments/option_minseqlength_1.md)

#(./fragments/option_topn.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_fastq_ascii.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_notrunclabels.md)

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

#(./fragments/option_fastq_qmax_ignored.md)

#(./fragments/option_fastq_qmin_ignored.md)

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Sort sequences in *input.fasta* by decreasing length and write to
*sorted.fasta*:

```sh
vsearch \
    --sortbylength input.fasta \
    --output sorted.fasta
```

Sort, discard sequences shorter than 200 nucleotides, and keep only
the 1000 longest:

```sh
vsearch \
    --sortbylength input.fasta \
    --minseqlength 200 \
    --topn 1000 \
    --output sorted.fasta
```


# SEE ALSO

[`vsearch-sortbysize(1)`](./vsearch-sortbysize.1.md),
[`vsearch-shuffle(1)`](./vsearch-shuffle.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
