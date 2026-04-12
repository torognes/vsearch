% vsearch-fastx_subsample(1) version 2.30.6 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastx_subsample --- randomly subsample fasta or fastq sequences


# SYNOPSIS

| **vsearch** **\-\-fastx_subsample** _fastxfile_ (**\-\-fastaout** | **\-\-fastqout**) _filename_ (**\-\-sample_pct** _real_ | **\-\-sample_size** _integer_) \[_options_]


# DESCRIPTION

The vsearch command `--fastx_subsample` randomly extracts a subset of
sequences from a fasta or fastq file, using a uniform distribution
without replacement. The subset size is specified either as a
percentage with `--sample_pct`, or as an absolute count with
`--sample_size`. To illustrate:

```text
Input (4 sequences):   --sample_size 2:

>s1                    sampled (--fastaout):
AAAA                   >s1
>s2         -->        AAAA
CCCC                   >s4
>s3                    GGGG
TTTT
>s4                    not sampled (--fastaout_discarded):
GGGG                   >s2
                       >s3
```

When `--sizein` is used, abundance annotations are taken into account:
sampling is performed as if the input sequences were first rereplicated
according to their abundances, then subsampled, then dereplicated. For
example, a sequence with `;size=10` contributes 10 slots to the
sampling pool. The selected sequences are written to `--fastaout`
and/or `--fastqout`. Sequences not selected can be written to
`--fastaout_discarded` and/or `--fastqout_discarded`.

Use `--randseed` to set a fixed seed for reproducible results.

See [`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md) and
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md) for more
information on input formats.


# OPTIONS

## mandatory options

At least one of `--fastaout` or `--fastqout` must be specified.
Exactly one of `--sample_pct` or `--sample_size` must be specified.


## core options

`--fastaout` *filename*
: Write the selected sequences to *filename*, in fasta format.

`--fastaout_discarded` *filename*
: Write the sequences not selected to *filename*, in fasta format.

`--fastqout` *filename*
: Write the selected sequences to *filename*, in fastq format.
  Requires input in fastq format (see
  [`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)).

`--fastqout_discarded` *filename*
: Write the sequences not selected to *filename*, in fastq format.
  Requires input in fastq format.

#(./fragments/option_randseed.md)

#(./fragments/option_sample_pct.md)

#(./fragments/option_sample_size.md)

#(./fragments/option_sizein.md)

#(./fragments/option_sizeout.md)


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

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Extract 10% of the sequences from *input.fasta* and write them to
*sub.fasta*:

```sh
vsearch \
    --fastx_subsample input.fasta \
    --sample_pct 10.0 \
    --fastaout sub.fasta
```

Extract exactly 1000 sequences, using a fixed seed for reproducibility,
and write those not selected to a separate file:

```sh
vsearch \
    --fastx_subsample input.fasta \
    --sample_size 1000 \
    --randseed 42 \
    --fastaout sub.fasta \
    --fastaout_discarded rest.fasta
```

Subsample accounting for abundance (10% of the total read count):

```sh
vsearch \
    --fastx_subsample input.fasta \
    --sizein \
    --sizeout \
    --sample_pct 10.0 \
    --fastaout sub.fasta
```


# SEE ALSO

[`vsearch-rereplicate(1)`](./vsearch-rereplicate.1.md),
[`vsearch-sortbysize(1)`](./vsearch-sortbysize.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
