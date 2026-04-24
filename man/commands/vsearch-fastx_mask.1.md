% vsearch-fastx_mask(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastx_mask --- mask low-complexity regions in fasta or fastq sequences


# SYNOPSIS

| **vsearch** **\-\-fastx_mask** _fastxfile_ (**\-\-fastaout** | **\-\-fastqout**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--fastx_mask` masks low-complexity regions and
simple repeats in sequences from a fasta or fastq file. The input
format is detected automatically (see
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md) and
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)).

Masking is controlled with the `--qmask` option, which accepts three
values: `dust` (default), `soft`, or `none`. With `dust`,
low-complexity regions are identified with the DUST algorithm
(R. Tatusov and D.J. Lipman, unpublished) and represented as lowercase
letters in the output. With `soft`, lowercase letters in the input are
treated as masked and preserved as-is. With `none`, no masking is
applied.

By default, masked regions are written as lowercase letters (soft
masking). When `--hardmask` is specified, masked regions are replaced
with N's instead.

Sequences can be filtered by the proportion of unmasked residues after
masking. Use `--min_unmasked_pct` to discard sequences with too few
unmasked residues, and `--max_unmasked_pct` to discard those with too
many.

At least one output option must be specified: `--fastaout` and/or
`--fastqout`. The `--fastqout` option requires fastq input.

This command supersedes `--maskfasta` (see
[`vsearch-maskfasta(1)`](./vsearch-maskfasta.1.md)), which is
restricted to fasta input.

To illustrate the effect of DUST masking on a fasta file:

```text
Input:                    --qmask dust:             --qmask dust --hardmask:

>s1                       >s1                       >s1
AAAAAAAAACGTACGT   -->    aaaaaaaaACGTACGT   -->    NNNNNNNNACGTACGT
>s2                       >s2                       >s2
ACGTACGTACGTACGT          ACGTACGTACGTACGT          ACGTACGTACGTACGT
```


# OPTIONS

## mandatory options

At least one of the following output options must be specified:

#(./fragments/option_fastaout_mask.md)

#(./fragments/option_fastqout_mask.md)


## core options

### masking

#(./fragments/option_hardmask.md)

#(./fragments/option_qmask.md)


### filtering

#(./fragments/option_max_unmasked_pct.md)

#(./fragments/option_min_unmasked_pct.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_fastq_ascii.md)

#(./fragments/option_fastq_qmax.md)

#(./fragments/option_fastq_qmin.md)

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

#(./fragments/option_threads.md)

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


# EXAMPLES

Mask low-complexity regions in a fasta file using the DUST algorithm
(default), and write the results to *masked.fasta*:

```sh
vsearch \
    --fastx_mask input.fasta \
    --fastaout masked.fasta
```

Replace masked regions with N's instead of lowercase letters
(`--hardmask`), and write the results to *masked.fasta*:

```sh
vsearch \
    --fastx_mask input.fasta \
    --hardmask \
    --fastaout masked.fasta
```

Mask a fastq file and discard sequences where more than 50% of
residues are unmasked (i.e., retain only highly repetitive sequences):

```sh
vsearch \
    --fastx_mask input.fastq \
    --max_unmasked_pct 50.0 \
    --fastqout masked.fastq
```

Mask a fasta file and discard sequences where fewer than 80% of
residues are unmasked (i.e., retain only sequences with sufficient
complexity):

```sh
vsearch \
    --fastx_mask input.fasta \
    --min_unmasked_pct 80.0 \
    --fastaout masked.fasta
```


# SEE ALSO

[`vsearch-maskfasta(1)`](./vsearch-maskfasta.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
