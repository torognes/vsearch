% vsearch-fastx_syncpairs(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-fastx_syncpairs --- reorder paired reads so that mates
appear in the same order in both files


# SYNOPSIS

| **vsearch** **\-\-fastx_syncpairs** _fastxfile_ **\-\-reverse** _fastxfile_ (**\-\-fastaout** | **\-\-fastqout**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--fastx_syncpairs` reorders paired-end reads so
that the two mates of each pair appear at the same position in the
forward and reverse files. Forward reads are given as the argument to
`--fastx_syncpairs`, and the reverse reads are specified with
`--reverse`. Both files must be in the same format, either fasta or
fastq (not a mix). See [`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md)
and [`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md) for more
information on these formats.

Two reads form a pair when they share the same *matching key*. The key
is the read header, truncated at the first whitespace, with a trailing
`/1` or `/2` mate marker removed (see `--read_separators`). For
example, the Casava 1.8+ headers `instrument:42 1:N:0:7` and
`instrument:42 2:N:0:7` share the key `instrument:42`, and the older
headers `instrument:42/1` and `instrument:42/2` share the key
`instrument:42`.

The reverse file is read once and kept in memory; the forward file is
then read entry by entry, and the matching reverse read is retrieved.
Synchronized pairs are written in the order of the forward file: the
forward reads to `--fastaout` and/or `--fastqout`, and their mates to
`--fastaout_rev` and/or `--fastqout_rev`. Neither input file needs to
be seekable, so both may be compressed or read from a pipe.

A read whose key appears in only one of the two files is called an
*orphan*. By default orphans are discarded. They can instead be
written out with `--fastaout_orphans` and `--fastqout_orphans`
(forward orphans, in forward-file order) and `--fastaout_orphans_rev`
and `--fastqout_orphans_rev` (reverse orphans, in reverse-file order).

Read labels are expected to be unique within each file. A duplicated
label is reported as an error when it makes the pairing ambiguous:
always for the reverse file, and for the forward file when two reads
claim the same reverse mate.

Note that reads are only reordered: they are not joined, merged, or
modified. To join non-overlapping paired-end reads, see
[`vsearch-fastq_join(1)`](./vsearch-fastq_join.1.md). To merge
overlapping paired-end reads, see
[`vsearch-fastq_mergepairs(1)`](./vsearch-fastq_mergepairs.1.md).


# OPTIONS

## mandatory options

`--reverse` *filename*
: Specify the fasta or fastq file containing the reverse reads. The
  format must match that of the forward reads.

At least one output file must also be specified.


## core options

`--fastaout` *filename*
: Write the synchronized forward reads to *filename*, in fasta format
  (see [`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md)). Quality
  scores are not written.

`--fastaout_rev` *filename*
: Write the synchronized reverse reads to *filename*, in fasta format.

`--fastqout` *filename*
: Write the synchronized forward reads to *filename*, in fastq format
  (see [`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)). The input
  must be in fastq format.

`--fastqout_rev` *filename*
: Write the synchronized reverse reads to *filename*, in fastq format.
  The input must be in fastq format.

`--fastaout_orphans` *filename*
: Write the forward reads that have no reverse mate to *filename*, in
  fasta format.

`--fastaout_orphans_rev` *filename*
: Write the reverse reads that have no forward mate to *filename*, in
  fasta format.

`--fastqout_orphans` *filename*
: Write the forward reads that have no reverse mate to *filename*, in
  fastq format. The input must be in fastq format.

`--fastqout_orphans_rev` *filename*
: Write the reverse reads that have no forward mate to *filename*, in
  fastq format. The input must be in fastq format.

`--read_separators` *string*
: Specify the set of characters that, when immediately followed by the
  mate number `1` or `2` at the end of the matching key, mark a
  trailing mate marker to be stripped. The default is `/`, which
  removes the `/1` and `/2` markers of older Illumina headers. A
  whitespace always ends the matching key, regardless of this option,
  so Casava 1.8+ headers are handled without any configuration.


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_quiet.md)


## ignored options

Reads are reordered without inspecting their quality scores, so the
following quality-related options are accepted (for instance in a
pipeline) but have no effect.

`--fastq_ascii` 33|64
: Option is ignored and has no effect.

#(./fragments/option_fastq_qmax_ignored.md)

#(./fragments/option_fastq_qmin_ignored.md)

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Reorder the reads of *R1.fastq* and *R2.fastq* so that mates line up,
writing the synchronized reads to *R1_sorted.fastq* and
*R2_sorted.fastq*:

```sh
vsearch \
    --fastx_syncpairs R1.fastq \
    --reverse R2.fastq \
    --fastqout R1_sorted.fastq \
    --fastqout_rev R2_sorted.fastq
```

Same as above, but also keep the unpaired reads in separate files:

```sh
vsearch \
    --fastx_syncpairs R1.fastq \
    --reverse R2.fastq \
    --fastqout R1_sorted.fastq \
    --fastqout_rev R2_sorted.fastq \
    --fastqout_orphans R1_orphans.fastq \
    --fastqout_orphans_rev R2_orphans.fastq
```

Synchronize reads whose headers carry the mate number after an
underscore (for example `readid_1` and `readid_2`):

```sh
vsearch \
    --fastx_syncpairs R1.fastq \
    --reverse R2.fastq \
    --read_separators _ \
    --fastqout R1_sorted.fastq \
    --fastqout_rev R2_sorted.fastq
```


# SEE ALSO

[`vsearch-fastq_join(1)`](./vsearch-fastq_join.1.md) for joining
non-overlapping paired-end reads,
[`vsearch-fastq_mergepairs(1)`](./vsearch-fastq_mergepairs.1.md) for
merging overlapping paired-end reads.
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
