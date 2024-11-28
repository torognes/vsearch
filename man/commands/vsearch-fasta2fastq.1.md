% vsearch-fasta2fastq(1) version 2.30.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fasta2fastq --- convert a fasta file to a fastq file with
fake quality scores


# SYNOPSIS

| **vsearch** **\-\-fasta2fastq** _fastafile_ \-\-fastqout _fastqfile_ \[_options_]


# DESCRIPTION

The vsearch command `--fasta2fastq` converts fasta sequences into
fastq sequences by adding fake quality scores. See
[`vsearch-fasta(5)`](./formats/vsearch-fasta.5.md) and
[`vsearch-fastq(5)`](./formats/vsearch-fastq.5.md) for more
information on these formats. Sequences are written to the file
specified with `--fastqout`. The quality score can be adjusted with
the `--fastq_qmaxout` option (default is 41). The quality offset can
be adjusted with the `--fastq_asciiout` option (either 33 or 64,
default is 33).


# OPTIONS

## mandatory options

`--fastqout` *filename*
: Write sequences to *filename* in fastq format, with fake quality
  scores. The default quality value (41) and offset (33) can be
  changed with options `--fastq_qmaxout` and `--fastq_asciiout`.


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

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Convert sequences in *input.fasta* into fastq sequences, with fake
quality values (Q = 40), and a quality offset of 64. Add sequence
length annotations (`--lengthout`). Write the results to
*output.fastq*:

```sh
vsearch \
    --fasta2fastq input.fasta \
    --fastq_qmaxout 40 \
    --fastq_asciiout 64 \
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
