% vsearch-sortbysize(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-sortbysize --- sort fasta or fastq sequences by decreasing
abundance


# SYNOPSIS

| **vsearch** **\-\-sortbysize** _fastxfile_ \-\-output _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--sortbysize` sorts fasta or fastq sequences by
decreasing abundance and writes them to the file specified with
`--output`, in fasta format. Sequences without abundance annotations
are assumed to have abundance 1. To illustrate:

```text
>s1;size=1                >s3;size=10
AAAA                      TTTT
>s2;size=5  --- sort -->  >s2;size=5
CCCC                      CCCC
>s3;size=10               >s1;size=1
TTTT                      AAAA
```

When two sequences have the same abundance, ties are broken by label
in alphanumerical order. Sequences can be filtered by abundance with
`--minsize` and `--maxsize`. The output can be limited to the first
*n* sequences with `--topn`.

See [`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md) and
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md) for more
information on input formats.


# OPTIONS

## mandatory options

#(./fragments/option_output.md)


## core options

#(./fragments/option_maxsize.md)

#(./fragments/option_minsize.md)

#(./fragments/option_topn.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_fastq_ascii.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_maxseqlength.md)

#(./fragments/option_minseqlength_1.md)

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

Sort sequences in *input.fasta* by decreasing abundance and write to
*sorted.fasta*:

```sh
vsearch \
    --sortbysize input.fasta \
    --output sorted.fasta
```

Sort, discard singletons (`--minsize 2`), and keep only the 1000 most
abundant sequences:

```sh
vsearch \
    --sortbysize input.fasta \
    --minsize 2 \
    --topn 1000 \
    --output sorted.fasta
```


# SEE ALSO

[`vsearch-sortbylength(1)`](./vsearch-sortbylength.1.md),
[`vsearch-shuffle(1)`](./vsearch-shuffle.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
