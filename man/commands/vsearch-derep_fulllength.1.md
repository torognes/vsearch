% vsearch-derep_fulllength(1) version 2.30.4 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-derep_fulllength --- merge strictly identical fasta or fastq
sequences


# SYNOPSIS

| **vsearch** **\-\-derep_fulllength** _fastxfile_ \-\-output _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--derep_fulllength` groups strictly identical
sequences from *fastxfile* and writes the unique sequences to
`--output`, in fasta format, sorted by decreasing abundance. Two
sequences are identical if they have the same length and the same
nucleotide string (comparison is case-insensitive; T and U are treated
as equivalent). Each unique sequence retains the header of the first
occurrence in the input. To illustrate:

```text
Input:           Output (--sizeout):

>s1              >s1;size=2
AAAA             AAAA
>s2    -->       >s3
AAAA             TTTT
>s3
TTTT
```

By default, only the plus strand is compared. Use `--strand both` to
also group sequences with their reverse complements.

Use `--sizein` to take existing abundance annotations into account when
counting. Use `--sizeout` to add abundance annotations to output
sequences. Use `--maxuniquesize` and `--minuniquesize` to filter on
post-dereplication abundance. Additional output in uclust-like format
is available with `--uc`.

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

Dereplicate sequences in *input.fasta*, annotating each unique sequence
with its abundance, and write results to *derep.fasta*:

```sh
vsearch \
    --derep_fulllength input.fasta \
    --sizeout \
    --output derep.fasta
```

Dereplicate, taking existing abundance annotations into account, and
discard singletons (`--minuniquesize 2`):

```sh
vsearch \
    --derep_fulllength input.fasta \
    --sizein \
    --sizeout \
    --minuniquesize 2 \
    --output derep.fasta
```


# SEE ALSO

[`vsearch-derep_id(1)`](./vsearch-derep_id.1.md),
[`vsearch-derep_prefix(1)`](./vsearch-derep_prefix.1.md),
[`vsearch-derep_smallmem(1)`](./vsearch-derep_smallmem.1.md),
[`vsearch-fastx_uniques(1)`](./vsearch-fastx_uniques.1.md),
[`vsearch-rereplicate(1)`](./vsearch-rereplicate.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
