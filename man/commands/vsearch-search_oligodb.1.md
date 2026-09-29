% vsearch-search_oligodb(1) version 2.33.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./fragments/date.md)

# NAME

vsearch \-\-search_oligodb --- find every occurrence of oligos (primers, tags, barcodes) in sequences


# SYNOPSIS

| **vsearch** **\-\-search_oligodb** _fastxfile_ **\-\-db** _fastafile_ (**\-\-alnout** | **\-\-blast6out** | **\-\-userout**) _filename_ \[_options_]


# DESCRIPTION

The vsearch command `--search_oligodb` looks for short oligonucleotides
(primers, tags, barcodes, adapters, splints), given in a fasta file
with `--db`, in the sequences of a fasta or fastq file, and reports
**every occurrence** of every oligo, with its position, on both strands.
It was written for long reads that carry several copies of a target
(rolling-circle amplification, concatemers, a barcode at each end), but
works on sequences of any length, from a few nucleotides to hundreds of
kilobases. The sequences are read one after the other and never held in
memory all at once.

The oligo is aligned end to end, and the sequence is free at both ends
(a *semi-global* alignment). An occurrence is reported when its
alignment has at most `--maxdiffs` differences (2 by default), a
difference being a substitution, an inserted nucleotide or a deleted
nucleotide. Among the alignments with the fewest differences, the one
with the fewest gap openings is reported. `--maxgaps` bounds the gap
openings, and `--maxgaps 0` restricts the search to substitutions.

Nucleotides are compared by overlap of their IUPAC base sets, on both
sides: an oligo `Y` matches a `C` or a `T` in the sequence, a sequence
`R` matches an oligo `A`, `G`, `R`, `N`... In particular, an `N` in the
sequence matches any nucleotide. As a consequence, a run of Ns at least
as long as an oligo is an occurrence without difference of that oligo,
and a site with a few Ns needs fewer differences to be reported. Use
`--n_mismatch` to count any position involving an N as a mismatch.
Case is ignored, and `U` is read as `T`.

An occurrence of an oligo is one site: the alignments ending near one
another (closer than half the oligo length) are collapsed into the one
with the fewest differences, the leftmost in case of a tie. Occurrences
of different oligos are all reported, even when they overlap.

Occurrences truncated by a sequence end (a barcode on a read that
starts or ends in the middle of it, for instance) are reported, as long
as at least `--target_cov` of the oligo (0.75 by default) lies inside
the sequence. The part of the oligo off the end is not counted as
differences; `tlo` > 1 or `thi` < `tl` shows which end is missing.

Occurrences are written in the order of the sequences, and within a
sequence by position. This command is multi-threaded, and its output is
the same whatever the number of threads.


## Positions

For this command, the userfields `qlo` and `qhi` report the first and
last positions of the occurrence in the sequence (the query), and `tlo`
and `thi` the first and last positions of the oligo (the target) that
are aligned. The columns 7 to 10 of `--blast6out` report the same four
values. For an occurrence on the minus strand, `qlo` and `qhi` are
positions on the plus strand, swapped (`qlo` > `qhi`), as vsearch does
for all minus-strand hits. The alignment fields `caln`, `aln`, `qrow`
and `trow` cover the occurrence only. (In the other commands, where the
alignments are global, `qlo`, `qhi`, `tlo` and `thi` span the whole
sequences.) See
[`vsearch-userfields(7)`](../misc/vsearch-userfields.7.md).

The field `raw` reports the number of differences of the alignment (as
`diffs` does), `bits` is 0 and `evalue` is -1.


## Random occurrences

The number of occurrences found by chance grows very fast with
`--maxdiffs`. Occurrences of one primer in uniform random sequence, per
100 kb of sequence searched on both strands:

| primer                              | `--maxdiffs 2` | 3    | 4   |
|-------------------------------------|----------------|------|-----|
| 515F, 19 nt, 2 ambiguous positions  | 0.02           | 0.44 | 6.9 |
| 515F, 19 nt, no ambiguous position  | 0.00           | 0.13 | 2.6 |
| 806R, 20 nt, 3 ambiguous positions  | 0.02           | 0.55 | 8.4 |

Real sequences are not random (low-complexity regions, repeats), so
these rates are a floor, not an estimate. For a 20-nucleotide primer, 3
differences is a practical ceiling on long reads, and each ambiguous
position multiplies the random occurrences by 2 to 3. When possible,
search for a longer oligo (a primer and its tag together, for
instance), and use the `diffs` field to filter the output rather than a
large `--maxdiffs`.

Truncated occurrences are the other source: the part of an oligo left
inside a read can be short enough to match by chance. On a nanopore run
of 1.7 million reads searched for 96 barcodes of 24 nucleotides (5 of
them present) with `--maxdiffs 3`, truncated occurrences of the 91
absent barcodes, which can only be random, numbered 2,456,380 with
`--target_cov 0.5`, 1,181 with 0.75, and none with 1; the present
barcodes had 340,343, 38,879 and none. The whole-length occurrences of
the absent barcodes numbered 26 in all three cases.


## Differences with usearch

usearch has a `search_oligodb` command. vsearch's version differs on
purpose:

- gaps are allowed (usearch aligns without gaps, at any `-maxdiffs`);
  `--maxgaps 0` gives usearch's substitution-only matching;
- an N in the sequence always matches, where usearch drops a site that
  holds two Ns or more, whatever `-maxdiffs`;
- occurrences truncated by a sequence end are reported;
- occurrences are ordered by position within a sequence, not by score;
- on the minus strand, `qlo` and `qhi` are swapped (usearch reports
  them in increasing order, and swaps the oligo positions instead);
- `--blast6out` reports an e-value of -1 and a bit score of 0;
- the default masking is *none*, and the default strand is *both*.


# OPTIONS

## mandatory options

`--db` must be specified, along with at least one output option.

#(./fragments/option_db_search_oligodb.md)


## core options

#(./fragments/option_maxdiffs_search_oligodb.md)

#(./fragments/option_maxgaps_search_oligodb.md)

#(./fragments/option_n_mismatch.md)

#(./fragments/option_strand_search_oligodb.md)

#(./fragments/option_target_cov_search_oligodb.md)

#(./fragments/option_threads.md)


## output options

#(./fragments/option_alnout_search_oligodb.md)

#(./fragments/option_blast6out_search_oligodb.md)

#(./fragments/option_rowlen.md)

#(./fragments/option_userfields.md)

#(./fragments/option_userout.md)


## secondary options

#(./fragments/option_bzip2_decompress.md)

#(./fragments/option_gzip_decompress.md)

#(./fragments/option_hardmask.md)

#(./fragments/option_log.md)

#(./fragments/option_no_progress.md)

#(./fragments/option_notrunclabels.md)

#(./fragments/option_qmask_search_oligodb.md)

#(./fragments/option_quiet.md)

#(./fragments/option_sizein.md)


# EXAMPLES

Find every occurrence of a set of primers and barcodes in nanopore
reads, with up to three differences, and write one line per
occurrence:

```sh
vsearch \
    --search_oligodb reads.fastq.gz \
    --db oligos.fasta \
    --maxdiffs 3 \
    --userout hits.tsv \
    --userfields query+target+qstrand+qlo+qhi+tlo+thi+diffs
```

The same search, without gaps, as usearch would do it, and with the
alignments written out for checking:

```sh
vsearch \
    --search_oligodb reads.fastq.gz \
    --db oligos.fasta \
    --maxdiffs 3 \
    --maxgaps 0 \
    --alnout hits.aln \
    --blast6out hits.tsv
```

Convert the occurrences to BED (0-based start, end excluded, both
positions on the plus strand):

```sh
vsearch \
    --search_oligodb reads.fastq \
    --db oligos.fasta \
    --userout hits.tsv \
    --userfields query+qlo+qhi+target+diffs+qstrand

awk 'BEGIN {OFS = "\t"}
     {lo = ($2 < $3) ? $2 : $3 ; hi = ($2 < $3) ? $3 : $2 ;
      split($1, label, " ") ;
      print label[1], lo - 1, hi, $4, $5, $6}' hits.tsv > hits.bed
```


# SEE ALSO

[`vsearch-cut(1)`](./vsearch-cut.1.md),
[`vsearch-search_exact(1)`](./vsearch-search_exact.1.md),
[`vsearch-usearch_global(1)`](./vsearch-usearch_global.1.md),
[`vsearch-cigar(5)`](../formats/vsearch-cigar.5.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md),
[`vsearch-nucleotides(7)`](../misc/vsearch-nucleotides.7.md),
[`vsearch-usearch(7)`](../misc/vsearch-usearch.7.md),
[`vsearch-userfields(7)`](../misc/vsearch-userfields.7.md)


#(./fragments/footer.md)
