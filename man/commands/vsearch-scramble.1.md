% vsearch-scramble(1) version 2.31.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-scramble --- randomize the nucleotide order within each
fasta or fastq entry


# SYNOPSIS

| **vsearch** **\-\-scramble** _fastxfile_ (**\-\-fastaout** | **\-\-fastqout**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--scramble` reads the entries of a fasta or
fastq file (*fastxfile*) one by one and (pseudo-)randomizes the order
of the nucleotides within each sequence. Output is written with
`--fastaout`, `--fastqout`, or both. At least one output option is
required. If the input is in fasta format, `--fastqout` cannot be
used, as there are no base quality scores to carry over.

Scrambling is the within-entry counterpart of shuffling: the vsearch
command `--shuffle` randomizes the order of the *entries* and leaves
each sequence untouched, whereas `--scramble` keeps the entries in
their input order and randomizes the nucleotides *inside* each
sequence (see
[`vsearch-shuffle(1)`](./vsearch-shuffle.1.md)). To illustrate:

```text
>s1                      >s2                      >s1
AAACGT                   TGGGCA                   ACATAG
>s2   --- shuffle -->    >s1   --- scramble -->   >s2
TGGGCA                   AAACGT                   GAGCGT
```

Scrambled datasets are useful as null models: they preserve the
number of entries, the entry order, the headers and abundance
annotations, the length of each sequence, and the exact nucleotide
composition of each sequence, while destroying all positional
structure. Each output sequence is a permutation of the bytes of the
input sequence, so uppercase and lowercase symbols, ambiguous
nucleotide symbols (see
[`vsearch-nucleotides(7)`](../misc/vsearch-nucleotides.7.md)), and
their per-sequence counts all pass through exactly. Note that
positional features are *not* preserved: runs of lowercase-masked
nucleotides are dispersed (only the masked fraction is preserved),
and the composition in words of two or more nucleotides changes (see
`--scramble_kmer`).

For fastq input, the quality string of each entry is always copied
through unchanged: the positional quality profile of each entry
(and therefore its expected error) is preserved exactly, and the
pairing between each nucleotide and its quality value is deliberately
broken.

Sequences of length 0 or 1, and sequences consisting of a single
repeated symbol, pass through unchanged, as they admit only one
arrangement.

For reproducibility, the seed for the pseudo-random generator can be
set with the option `--randseed`; a given seed yields the same result
on any platform. Each entry is scrambled with its own random
sub-stream derived from the seed and the entry ordinal, so for a
given seed an entry's scramble does not depend on the other entries,
nor on the presence or absence of quality values.


# OPTIONS

## mandatory options

At least one of the following output options is required:

#(./fragments/option_fastaout_scramble.md)

#(./fragments/option_fastqout_scramble.md)


## core options

#(./fragments/option_randseed.md)

#(./fragments/option_scramble_kmer.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_fasta_width.md)

#(./fragments/option_fastq_ascii.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_label_suffix.md)

#(./fragments/option_lengthout.md)

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

`--notrunclabels`
: Retain whole sequence headers in output files. With the vsearch
  command `--scramble`, sequence headers are never truncated, so
  this option has no visible effect.

#(./fragments/option_quiet.md)

#(./fragments/option_relabel.md)

#(./fragments/option_relabel_keep.md)

#(./fragments/option_relabel_md5.md)

#(./fragments/option_relabel_self.md)

#(./fragments/option_relabel_sha1.md)

#(./fragments/option_sample.md)

#(./fragments/option_sizein.md)
: Has no effect with this command: abundance annotations are always
  parsed (`--sizeout` therefore writes the true abundances whether or
  not `--sizein` is given).

#(./fragments/option_sizeout.md)

#(./fragments/option_xee.md)

#(./fragments/option_xlength.md)

#(./fragments/option_xsize.md)


## ignored options

#(./fragments/option_fastq_qmax_ignored.md)

#(./fragments/option_fastq_qmin_ignored.md)

#(./fragments/option_threads_not_multithreaded.md)


# EXAMPLES

Scramble a fasta file, with a fixed seed for reproducibility:

```sh
vsearch \
    --scramble input.fasta \
    --randseed 1 \
    --fastaout scrambled.fasta
```

Scramble a fastq file and write both fasta and fastq output (quality
strings are copied through unchanged):

```sh
vsearch \
    --scramble input.fastq \
    --fastaout scrambled.fasta \
    --fastqout scrambled.fastq
```

# SEE ALSO

[`vsearch-shuffle(1)`](./vsearch-shuffle.1.md),
[`vsearch-fastx_revcomp(1)`](./vsearch-fastx_revcomp.1.md),
[`vsearch-fastx_subsample(1)`](./vsearch-fastx_subsample.1.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(./fragments/footer.md)
