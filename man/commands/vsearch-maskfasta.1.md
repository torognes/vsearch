% vsearch-maskfasta(1) version 2.30.6 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-maskfasta --- mask low-complexity regions in fasta sequences (deprecated)


# SYNOPSIS

| **vsearch** **\-\-maskfasta** _fastafile_ **\-\-output** _filename_ \[_options_]


# DESCRIPTION

**This command is deprecated. Use `--fastx_mask` instead** (see
[`vsearch-fastx_mask(1)`](./vsearch-fastx_mask.1.md)), which accepts
both fasta and fastq input.

The vsearch command `--maskfasta` masks low-complexity regions and
simple repeats in sequences from a fasta file (see
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md)). The output file
is specified with `--output`, in fasta format.

Masking is controlled with the `--qmask` option, which accepts three
values: `dust` (default), `soft`, or `none`. With `dust`, low-complexity
regions are identified with the DUST algorithm and represented as
lowercase letters in the output. With `soft`, lowercase letters in the
input are treated as masked and preserved as-is. With `none`, no masking
is applied.

By default, masked regions are written as lowercase letters (soft
masking). When `--hardmask` is specified, masked regions are replaced
with N's instead.


# OPTIONS

## mandatory options

#(./fragments/option_output.md)


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

#(./fragments/option_threads.md)

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


# EXAMPLES

Mask low-complexity regions in a fasta file using the DUST algorithm
(default), and write the results to *masked.fasta*:

```sh
vsearch \
    --maskfasta input.fasta \
    --qmask dust \
    --output masked.fasta
```

Replace masked regions with N's instead of lowercase letters
(`--hardmask`):

```sh
vsearch \
    --maskfasta input.fasta \
    --hardmask \
    --output masked.fasta
```


# SEE ALSO

[`vsearch-fastx_mask(1)`](./vsearch-fastx_mask.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md)


#(./fragments/footer.md)
