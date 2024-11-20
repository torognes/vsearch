% vsearch-fasta2fastq(1) version 2.30.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch `--fasta2fastq` --- convert a fasta file to a fastq file with
fake quality scores


# SYNOPSIS

| **vsearch** **--fasta2fastq** _fastafile_ \-\-fastqout _fastqfile_ \[_options_]


# DESCRIPTION

The `vsearch --fasta2fastq` command can be used to add fake nucleotide
quality scores to the sequences in the given fasta file and write them
to the fastq file specified with `--fastqout`. The quality score can
be adjusted with the `--fastq_qmaxout` option (default is 41). The
fastq output quality offset, in ASCII base character, can be adjusted
with the `--fastq_asciiout` option (default is 33).


# OPTIONS

## mandatory parameters

`--fastqout` *filename*
: Write sequences to *filename* in fastq format, with fake quality
  strings. The default quality value (41) can be changed with option
  `--fastq_qmaxout`.


## core options

#(./fragments/option_fastq_asciiout.md)

`--fastq_qmaxout` *positive integer*
: Specify the fake quality score used when writing the fastq file. The
  default is 41, which is the usual maximal quality score for recent
  Sanger/Illumina 1.8+ files (maximal quality score was 40 in older
  formats). Accepted values range from 0 to 93 when the quality offset
  is set to 33 (see `--fastq_asciiout`), and from 0 to 62 when the
  quality offset is set to 64.


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_quiet.md)

#(./fragments/option_relabel.md)

#(./fragments/option_relabel_keep.md)

#(./fragments/option_relabel_md5.md)

#(./fragments/option_relabel_self.md)

#(./fragments/option_relabel_sha1.md)

#(./fragments/option_sample.md)

#(./fragments/option_sizein.md)

#(./fragments/option_sizeout.md)

#(./fragments/option_threads_not_multithreaded.md)

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


# EXAMPLES

Add strings of fake quality values Q = 40 to *input.fasta*
sequences. Write the results to *output.fastq*, in fastq format with a
+64 quality offset:

```sh
vsearch \
    --fasta2fastq input.fasta \
    --fastq_asciiout 64 \
    --fastq_qmaxout 40 \
    --lengthout \
    --fastqout output.fastq
```

For instance, the following fasta input:

```text
>s1
CACCGCGGTTATACGAGGGGCTCAAATTGATATT
>s2
CACCGCGGTTATACGAGGGGCTCAAATTGATATT
AATATCAATTTGAGCCCCTCGTATAACCGCGGTG
```

becomes:

```text
@s1;length=34
CACCGCGGTTATACGAGGGGCTCAAATTGATATT
+
hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
@s2;length=68
CACCGCGGTTATACGAGGGGCTCAAATTGATATTAATATCAATTTGAGCCCCTCGTATAACCGCGGTG
+
hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
```

(note that fastq output files are not wrapped)


#(./fragments/footer.md)
