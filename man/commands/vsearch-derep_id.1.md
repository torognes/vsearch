% vsearch-derep_id(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-derep_id --- merge identical fasta sequences sharing the
same label


# SYNOPSIS

| **vsearch** **\-\-derep_id** _fastafile_ \-\-output _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--derep_id` groups sequences that are both
sequence-identical and share the same label (identifier), and writes
the unique sequences to `--output`, in fasta format, sorted by
decreasing abundance. It behaves like `--derep_fulllength` (see
[`vsearch-derep_fulllength(1)`](./vsearch-derep_fulllength.1.md)),
with the additional requirement that the sequence label (the first
word of the header line, or the whole line under `--notrunclabels`)
must be identical across all sequences in a group. Abundance
annotations embedded in the label (such as `;size=`) count as part of
it: `>a;size=5` and `>a;size=9` are different labels. Sequences with
identical nucleotides but different labels are not grouped.

This command is not multithreaded. The inverse operation is
`--rereplicate` (see
[`vsearch-rereplicate(1)`](./vsearch-rereplicate.1.md)).

Fastq input is rejected: use `--fastx_uniques` (see
[`vsearch-fastx_uniques(1)`](./vsearch-fastx_uniques.1.md)) to
dereplicate fastq files. See
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md) for a description
of the input format.


# OPTIONS

## mandatory options

At least one of `--output` or `--uc` (see core options) must be
specified.

#(./fragments/option_output.md)


## core options

#(./fragments/option_maxuniquesize.md)

#(./fragments/option_minuniquesize.md)

#(./fragments/option_sizein.md)

#(./fragments/option_sizeout_derep.md)

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

#(./fragments/option_minseqlength_32.md)

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

Group sequences that are both sequence-identical and share the same
label:

```sh
vsearch \
    --derep_id input.fasta \
    --sizeout \
    --output derep.fasta
```


# SEE ALSO

[`vsearch-derep_fulllength(1)`](./vsearch-derep_fulllength.1.md),
[`vsearch-derep_prefix(1)`](./vsearch-derep_prefix.1.md),
[`vsearch-derep_smallmem(1)`](./vsearch-derep_smallmem.1.md),
[`vsearch-fastx_uniques(1)`](./vsearch-fastx_uniques.1.md),
[`vsearch-rereplicate(1)`](./vsearch-rereplicate.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
