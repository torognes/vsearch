% vsearch-sff_convert(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-sff_convert --- convert an SFF file to fastq


# SYNOPSIS

| **vsearch** **\-\-sff_convert** _sfffile_ **\-\-fastqout** _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--sff_convert` converts the reads stored in an SFF
(Standard Flowgram Format) file to fastq. SFF is a binary format used by
Roche 454 and early Ion Torrent PGM sequencing platforms (see
[`vsearch-sff(5)`](../formats/vsearch-sff.5.md)).

The output fastq file is written to `--fastqout`. The quality encoding
offset can be set with `--fastq_asciiout` (default 33, phred+33). Output
quality scores can be clamped with `--fastq_qminout` and `--fastq_qmaxout`.

Each SFF read stores clipping coordinates that indicate low-quality or
adapter regions at the ends of the sequence. By default, no clipping is
applied: the full sequence is written, with bases that would be clipped
converted to lower case and the rest in upper case. Use `--sff_clip` to
apply the clipping and write only the retained region in upper case.


# OPTIONS

## mandatory options

#(./fragments/option_fastqout_convert.md)


## core options

#(./fragments/option_fastq_asciiout.md)

#(./fragments/option_fastq_qmaxout.md)

#(./fragments/option_fastq_qminout.md)

#(./fragments/option_sff_clip.md)


## secondary options

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

#(./fragments/option_sizeout.md)


## ignored options

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Convert an SFF file to fastq using the default phred+33 encoding:

```sh
vsearch \
    --sff_convert reads.sff \
    --fastqout reads.fastq
```

Convert and apply clipping coordinates from the SFF file:

```sh
vsearch \
    --sff_convert reads.sff \
    --sff_clip \
    --fastqout reads_clipped.fastq
```

Convert with phred+64 encoding and restrict quality scores to 0--40:

```sh
vsearch \
    --sff_convert reads.sff \
    --fastq_asciiout 64 \
    --fastq_qmaxout 40 \
    --fastqout reads.fastq
```


# SEE ALSO

[`vsearch-sff(5)`](../formats/vsearch-sff.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md),
[`vsearch-fastq_convert(1)`](./vsearch-fastq_convert.1.md)


#(./fragments/footer.md)
