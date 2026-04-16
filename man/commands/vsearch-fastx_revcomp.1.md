% vsearch-fastx_revcomp(1) version 2.30.6 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastx_revcomp --- reverse-complement fasta or fastq sequences


# SYNOPSIS

| **vsearch** **\-\-fastx_revcomp** _fastxfile_ (**\-\-fastaout** | **\-\-fastqout**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--fastx_revcomp` reverse-complements the sequences in a
fasta or fastq file. Output is written with `--fastaout`, `--fastqout`, or
both. At least one output option is required.

If the input is in fasta format, `--fastqout` cannot be used, as there are no
base quality scores to carry over.

To illustrate with a short sequence:

```text
Input:             Reverse complement:

>s1                >s1
AACGT              ACGTT
```


# OPTIONS

## mandatory options

At least one of the following output options is required:

#(./fragments/option_fastaout_revcomp.md)

#(./fragments/option_fastqout_revcomp.md)


## core options

None


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

`--notrunclabels`
: Retain whole sequence headers in output files. With the vsearch
  command `--fastx_revcomp`, sequence headers are never truncated, so
  this option has no visible effect.

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

`--fastq_ascii` 33|64
: Option is ignored and has no effect.

#(./fragments/option_fastq_qmax_ignored.md)

#(./fragments/option_fastq_qmin_ignored.md)

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Reverse-complement a fasta file:

```sh
vsearch \
    --fastx_revcomp input.fasta \
    --fastaout revcomp.fasta
```

Reverse-complement a fastq file and write both fasta and fastq output:

```sh
vsearch \
    --fastx_revcomp input.fastq \
    --fastaout revcomp.fasta \
    --fastqout revcomp.fastq
```


# SEE ALSO

[`vsearch-fastx_filter(1)`](./vsearch-fastx_filter.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
