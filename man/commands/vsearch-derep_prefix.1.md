% vsearch-derep_prefix(1) version 2.30.4 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-derep_prefix --- merge fasta or fastq sequences with identical
prefixes


# SYNOPSIS

| **vsearch** **\-\-derep_prefix** _fastxfile_ \-\-output _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--derep_prefix` groups sequences based on prefix
identity. A shorter sequence that is identical to an initial segment
(prefix) of a longer sequence is considered a replicate of the longer
one. Comparison is case-insensitive; T and U are treated as
equivalent. To illustrate:

```text
Input:                Output (--sizeout):

>long                 >long;size=2
AAAACCCG              AAAACCCG
>short    -->         >other
AAAA                  TTTT
>other
TTTT
```

Here, 'short' (AAAA) is a prefix of 'long' (AAAACCCG), so it is
grouped with 'long'.

If a short sequence is a prefix of multiple longer sequences, it is
grouped with the shortest. Ties in length are resolved by decreasing
abundance, then by header in alphanumerical order, then by input
order. The unique sequences are written to `--output`, in fasta
format, sorted by decreasing abundance.

This command is not multithreaded. The inverse operation is
`--rereplicate` (see
[`vsearch-rereplicate(1)`](./vsearch-rereplicate.1.md)).

See [`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md) and
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md) for a description
of the input formats.


# OPTIONS

## mandatory options

#(./fragments/option_output.md)


## core options

#(./fragments/option_maxuniquesize.md)

#(./fragments/option_minuniquesize.md)

#(./fragments/option_sizein.md)

#(./fragments/option_sizeout.md)

#(./fragments/option_strand.md)

#(./fragments/option_topn.md)

#(./fragments/option_uc.md)


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

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Group sequences with identical prefixes, annotate with abundance, and
write to *derep.fasta*:

```sh
vsearch \
    --derep_prefix input.fasta \
    --sizeout \
    --output derep.fasta
```


# SEE ALSO

[`vsearch-derep_fulllength(1)`](./vsearch-derep_fulllength.1.md),
[`vsearch-derep_id(1)`](./vsearch-derep_id.1.md),
[`vsearch-derep_smallmem(1)`](./vsearch-derep_smallmem.1.md),
[`vsearch-fastx_uniques(1)`](./vsearch-fastx_uniques.1.md),
[`vsearch-rereplicate(1)`](./vsearch-rereplicate.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
