% vsearch-derep_smallmem(1) version 2.30.4 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-derep_smallmem --- merge strictly identical sequences using
minimal memory


# SYNOPSIS

| **vsearch** **\-\-derep_smallmem** _fastxfile_ \-\-fastaout _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--derep_smallmem` groups strictly identical
sequences, like `--derep_fulllength` (see
[`vsearch-derep_fulllength(1)`](./vsearch-derep_fulllength.1.md)),
but uses significantly less memory (approximately 24 bytes per unique
sequence). The output is written to `--fastaout`, in fasta format.

Key differences from `--derep_fulllength`:

- Output is written in the order that unique sequences first appear
  in the input, not sorted by decreasing abundance.
- Cannot read from a pipe; the input must be a regular file (it is
  read twice).
- Can read fastq files, but output is always fasta.
- Dereplication uses a 128-bit hash function; grouped sequences are
  not explicitly verified to be identical. The false positive
  probability is approximately 1e-21 for a dataset of one billion
  unique sequences.
- Does not support `--topn` or `--uc`.

This command is not multithreaded. The inverse operation is
`--rereplicate` (see
[`vsearch-rereplicate(1)`](./vsearch-rereplicate.1.md)).

See [`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md) and
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md) for a description
of the input formats.


# OPTIONS

## mandatory options

`--fastaout` *filename*
: Write the dereplicated sequences to *filename*, in fasta format, in
  the order they first appear in the input.


## core options

#(./fragments/option_maxuniquesize.md)

#(./fragments/option_minuniquesize.md)

#(./fragments/option_sizein.md)

#(./fragments/option_sizeout.md)

#(./fragments/option_strand.md)


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

Dereplicate sequences in *input.fasta* with minimal memory usage:

```sh
vsearch \
    --derep_smallmem input.fasta \
    --sizeout \
    --fastaout derep.fasta
```


# SEE ALSO

[`vsearch-derep_fulllength(1)`](./vsearch-derep_fulllength.1.md),
[`vsearch-derep_id(1)`](./vsearch-derep_id.1.md),
[`vsearch-derep_prefix(1)`](./vsearch-derep_prefix.1.md),
[`vsearch-fastx_uniques(1)`](./vsearch-fastx_uniques.1.md),
[`vsearch-rereplicate(1)`](./vsearch-rereplicate.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
