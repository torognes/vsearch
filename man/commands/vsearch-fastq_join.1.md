% vsearch-fastq_join(1) version 2.30.4 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastq_join --- join paired-end reads into one sequence with
a gap


# SYNOPSIS

| **vsearch** **\-\-fastq_join** _fastqfile_ **\-\-reverse** _fastqfile_ (**\-\-fastaout** | **\-\-fastqout**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--fastq_join` joins paired-end reads into a
single sequence. Forward reads are given as the argument to
`--fastq_join`, and the reverse reads are specified with `--reverse`.
The resulting sequence consists of the forward read, a padding
sequence, and the reverse complement of the reverse read. To
illustrate:

```text
>R1            >R2                  >R1
AAACCC    +    GGGAAA    -->   AAACCCNNNNNNNNTTTCCC
                               (R1) (padding) (revcomp R2)
```

The padding sequence defaults to `NNNNNNNN` (eight N's) and can be
changed with `--join_padgap`. The corresponding quality string
defaults to `IIIIIIII` (quality score 40, error probability 0.0001)
and can be changed with `--join_padgapq`. Joined sequences are written
to `--fastaout` and/or `--fastqout`. See
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md) and
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md) for more
information on these formats.

Note that reads are joined, not merged. To merge overlapping
paired-end reads, see
[`vsearch-fastq_mergepairs(1)`](./vsearch-fastq_mergepairs.1.md).


# OPTIONS

## mandatory options

#(./fragments/option_reverse.md)

At least one of `--fastaout` or `--fastqout` must also be specified.


## core options

`--fastaout` *filename*
: Write the joined sequences to *filename*, in fasta format (see
  [`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md)). Quality scores
  are not written.

`--fastqout` *filename*
: Write the joined sequences to *filename*, in fastq format (see
  [`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)).

#(./fragments/option_join_padgap.md)

#(./fragments/option_join_padgapq.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_fastq_ascii.md)

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

Join forward reads in *R1.fastq* with reverse reads in *R2.fastq*,
and write the joined sequences to *joined.fastq*:

```sh
vsearch \
    --fastq_join R1.fastq \
    --reverse R2.fastq \
    --fastqout joined.fastq
```

Use a custom padding sequence of twelve N's and a matching quality
string of twelve I's:

```sh
vsearch \
    --fastq_join R1.fastq \
    --reverse R2.fastq \
    --join_padgap NNNNNNNNNNNN \
    --join_padgapq IIIIIIIIIIII \
    --fastqout joined.fastq
```


# SEE ALSO

[`vsearch-fastq_mergepairs(1)`](./vsearch-fastq_mergepairs.1.md) for
merging overlapping paired-end reads.
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md),
[`vsearch-expected_error(7)`](../misc/vsearch-expected_error.7.md)


#(./fragments/footer.md)
