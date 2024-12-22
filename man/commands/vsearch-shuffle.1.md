% vsearch-shuffle(1) version 2.30.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-shuffle --- randomize the order of fasta or fastq entries


# SYNOPSIS

| **vsearch** **\-\-shuffle** _fastxfile_ \-\-output _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--shuffle` (pseudo-)randomizes the order of fasta
or fastq entries (*fastxfile*), and writes to `--output` *filename*,
in fasta format. To illustrate:

```text
>s1                     >s3
AAAA                    TTTT
>s2   --- shuffle -->   >s1
CCCC                    AAAA
>s3                     >s2
TTTT                    CCCC
```

The number of entries written to *filename* can be limited with the
option `--topn`. For reproducibility, the seed for the pseudo-random
generator can be set with the option `--randseed`.


# OPTIONS

## mandatory options

`--output` *filename*
: Write the fasta or fastq entries to *filename*, in fasta format and
  in a pseudo-randomized order.


## core options

#(./fragments/option_randseed.md)

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

Read from *input.fasta*, do not write to the standard error
(`--quiet`), use a pseudo-random seed (`--randseed` is not set), write
results to *output.fasta* (`--output`), and record shuffling
parameters in *output.log* (`--log`):

```sh
vsearch \
    --shuffle input.fasta \
    --quiet \
    --log output.log \
    --output output.fasta
```

# SEE ALSO

Inverse operation (sort):
[`vsearch-sortbylength(1)`](./commands/vsearch-sortbylength.1.md),
[`vsearch-sortbysize(1)`](../formats/vsearch-sortbysize.1.md)


#(./fragments/footer.md)
