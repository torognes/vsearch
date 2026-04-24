% vsearch-fastx_uniques(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastx_uniques --- merge strictly identical fasta or fastq
sequences


# SYNOPSIS

| **vsearch** **\-\-fastx_uniques** _fastxfile_ (**\-\-fastaout** | **\-\-fastqout**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--fastx_uniques` groups strictly identical
sequences from a fasta or fastq file, like `--derep_fulllength` (see
[`vsearch-derep_fulllength(1)`](./vsearch-derep_fulllength.1.md)).
The unique sequences are written to `--fastaout` and/or `--fastqout`,
sorted by decreasing abundance.

When the input is a fastq file, quality scores in the output
correspond by default to the average error probability at each
position across all grouped sequences (see
[`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md)).
Use `--fastq_qout_max` to use the best (highest) quality score
observed at each position instead.

This command is not multithreaded. The inverse operation is
`--rereplicate` (see
[`vsearch-rereplicate(1)`](./vsearch-rereplicate.1.md)).

See [`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md) and
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md) for a description
of the input formats.


# OPTIONS

## mandatory options

At least one of `--fastaout` or `--fastqout` must be specified.


## core options

`--fastaout` *filename*
: Write the dereplicated sequences to *filename*, in fasta format,
  sorted by decreasing abundance. Each unique sequence retains the
  header of the first occurrence in the input.

`--fastqout` *filename*
: Write the dereplicated sequences to *filename*, in fastq format (see
  [`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)), sorted by
  decreasing abundance. Quality scores correspond to the average error
  probability at each position across all grouped sequences, or to the
  best quality score if `--fastq_qout_max` is set.

`--fastq_qout_max`
: Use the best (highest) quality score observed at each position when
  computing quality scores for fastq output, instead of averaging
  error probabilities.

#(./fragments/option_maxuniquesize.md)

#(./fragments/option_minuniquesize.md)

#(./fragments/option_sizein.md)

#(./fragments/option_sizeout.md)

#(./fragments/option_strand.md)

#(./fragments/option_tabbedout.md)

#(./fragments/option_topn.md)

#(./fragments/option_uc.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_fastq_ascii.md)

#(./fragments/option_fastq_asciiout.md)

#(./fragments/option_fastq_qmax.md)

#(./fragments/option_fastq_qmaxout.md)

#(./fragments/option_fastq_qmin.md)

#(./fragments/option_fastq_qminout.md)

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

Dereplicate sequences in *input.fasta*, annotate with abundance, write
to *derep.fasta*, and record clustering details in *derep.uc*:

```sh
vsearch \
    --fastx_uniques input.fasta \
    --sizeout \
    --fastaout derep.fasta \
    --uc derep.uc
```

Dereplicate a fastq file, keeping only the best quality score at each
position:

```sh
vsearch \
    --fastx_uniques input.fastq \
    --fastq_qout_max \
    --sizeout \
    --fastqout derep.fastq
```


# SEE ALSO

[`vsearch-derep_fulllength(1)`](./vsearch-derep_fulllength.1.md),
[`vsearch-derep_id(1)`](./vsearch-derep_id.1.md),
[`vsearch-derep_prefix(1)`](./vsearch-derep_prefix.1.md),
[`vsearch-derep_smallmem(1)`](./vsearch-derep_smallmem.1.md),
[`vsearch-rereplicate(1)`](./vsearch-rereplicate.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md),
[`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md)


#(./fragments/footer.md)
