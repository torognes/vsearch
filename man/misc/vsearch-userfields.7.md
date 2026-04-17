% vsearch-userfields(7) version 2.30.6 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(../commands/fragments/date.md)

# NAME

userfields --- output fields available with the --userout option


# DESCRIPTION

The option `--userfields` selects and orders the columns written to a
`--userout` file. Fields are separated by `+`:

```sh
vsearch \
    --usearch_global queries.fasta \
    --db db.fasta \
    --id 0.97 \
    --userout results.tsv \
    --userfields query+target+id+alnlen+mism+opens
```

When a query has no match (reported only with `--output_no_hits`),
numeric fields are set to 0, string fields are empty, and `target` is
set to `*`, unless noted otherwise.

Fields are grouped thematically below. See also
[`vsearch-pairwise_alignment_parameters(7)`](./vsearch-pairwise_alignment_parameters.7.md)
for a description of the alignment model and identity definitions.


## Alignment representation

`aln`
: Pairwise alignment encoded as a string of operation characters: `M`
  (match or mismatch, i.e. not a gap), `D` (deletion, i.e. a gap in
  the query), `I` (insertion, i.e. a gap in the target). Empty if
  there is no alignment. See
  [`vsearch-cigar(5)`](../formats/vsearch-cigar.5.md) for a description
  of the operation alphabet.

`caln`
: Compact pairwise alignment in CIGAR format (Compact Idiosyncratic
  Gapped Alignment Report): `M` (match or mismatch), `D` (deletion),
  `I` (insertion). The equal sign `=` indicates that the query is
  identical to the centroid sequence (ignoring terminal gaps). Empty if
  there is no alignment. See
  [`vsearch-cigar(5)`](../formats/vsearch-cigar.5.md) for a complete
  description.

`qrow`
: Query segment as seen in the pairwise alignment, with gap characters
  inserted where the target has an insertion. Empty if there is no
  alignment.

`trow`
: Target segment as seen in the pairwise alignment, with gap characters
  inserted where the query has an insertion. Empty if there is no
  alignment.


## Query and target identifiers

`query`
: Query sequence label.

`target`
: Target sequence label. Set to `*` if there is no alignment.


## Sequence lengths

`ql`
: Query sequence length (positive integer).

`tl`
: Target sequence length (positive integer).

`qs`
: Query segment length. Always equal to the query sequence length.

`ts`
: Target segment length. Always equal to the target sequence length.


## Alignment span

The following fields report the first and last aligned positions in the
query (q) or target (t) sequence, using 1-based nucleotide positions.
The `lo`/`hi` variants always span the full alignment including terminal
gaps; the `ilo`/`ihi` variants exclude terminal gaps and report the
span of the actual aligned residues.

`qlo`
: First nucleotide of the query aligned with the target. Always 1 when
  there is an alignment (including terminal gaps at the left). See
  `qilo` to exclude initial gaps.

`qhi`
: Last nucleotide of the query aligned with the target. Always equal to
  the alignment length when there is an alignment (including terminal
  gaps at the right). See `qihi` to exclude terminal gaps.

`qilo`
: First nucleotide of the query aligned with the target, ignoring
  initial terminal gaps. Nucleotide positions use 1-based indexing.

`qihi`
: Last nucleotide of the query aligned with the target, ignoring
  terminal gaps at the right. Nucleotide positions use 1-based indexing.

`tlo`
: First nucleotide of the target aligned with the query. Always 1 when
  there is an alignment (including terminal gaps at the left). See
  `tilo` to exclude initial gaps.

`thi`
: Last nucleotide of the target aligned with the query. Always equal to
  the alignment length when there is an alignment (including terminal
  gaps at the right). See `tihi` to exclude terminal gaps.

`tilo`
: First nucleotide of the target aligned with the query, ignoring
  initial terminal gaps. Nucleotide positions use 1-based indexing.

`tihi`
: Last nucleotide of the target aligned with the query, ignoring
  terminal gaps at the right. Nucleotide positions use 1-based indexing.


## Alignment statistics

`alnlen`
: Length of the pairwise alignment (number of columns, including
  gap-only columns).

`ids`
: Number of matching columns in the alignment (zero or positive
  integer).

`mism`
: Number of mismatching columns in the alignment (zero or positive
  integer).

`gaps`
: Number of gap-containing columns in the alignment (zero or positive
  integer, excluding terminal gaps).

`opens`
: Number of gap-opening columns in the alignment (zero or positive
  integer, excluding terminal gaps).

`exts`
: Number of gap-extension columns in the alignment (zero or positive
  integer).

`pairs`
: Number of columns containing only nucleotides (alignment length minus
  gap-containing columns; zero or positive integer).

`pv`
: Number of positive columns. Equivalent to the number of matches for
  nucleotide sequences.

`pctgaps`
: Number of gap-containing columns expressed as a percentage of the
  alignment length (real value from 0.0 to 100.0).

`pctpv`
: Percentage of positive columns. Equivalent to the percentage of
  matches for nucleotide sequences (real value from 0.0 to 100.0).


## Identity percentages

`id`
: Percentage of identity computed according to the definition selected
  by `--iddef` (default: `id2`; real value from 0.0 to 100.0).

`id0`
: CD-HIT definition: 100 * (matching columns) / (shortest sequence
  length).

`id1`
: Edit distance: 100 * (matching columns) / (alignment length).

`id2`
: Edit distance excluding terminal gaps: 100 * (matching columns) /
  (alignment length - terminal gaps). Default definition for `--id`.

`id3`
: Marine Biological Lab definition, counting each gap opening (internal
  or terminal) as a single mismatch, whether or not the gap was
  extended: 100 * (1.0 - [(mismatches + gap openings) / (longest
  sequence length)]).

`id4`
: BLAST definition, equivalent to `id1` for global pairwise
  alignments. Always equal to `id1`.


## Coverage

`qcov`
: Fraction of the query sequence aligned with the target (real value
  from 0.0 to 100.0). Computed as 100 * (matches + mismatches) /
  (query sequence length). Gap-only columns are not counted.

`tcov`
: Fraction of the target sequence aligned with the query (real value
  from 0.0 to 100.0). Computed as 100 * (matches + mismatches) /
  (target sequence length). Gap-only columns are not counted.


## Score

`raw`
: Raw alignment score (negative, zero, or positive integer). The score
  is the sum of match rewards minus mismatch penalties, gap opening
  penalties, and gap extension penalties, using the parameters set by
  `--match`, `--mismatch`, `--gapopen`, and `--gapext`.

`bits`
: Bit score. Not computed for nucleotide alignments. Always 0.

`evalue`
: E-value. Not computed for nucleotide alignments. Always -1.


## Strand

`qstrand`
: Query strand orientation (`+` or `-` for nucleotide sequences). Empty
  if there is no alignment.

`tstrand`
: Target strand orientation. Always `+`: when a query matches a target
  on the reverse strand, `tstrand` is `+` and `qstrand` is `-`. Empty
  if there is no alignment.


## Unused fields

`qframe`
: Query reading frame (-3 to +3). Only meaningful for coding sequences;
  not computed by vsearch. Always `+0`.

`tframe`
: Target reading frame (-3 to +3). Only meaningful for coding
  sequences; not computed by vsearch. Always `+0`.


# EXAMPLES

Output query label, target label, identity, alignment length, and
number of mismatches:

```sh
vsearch \
    --usearch_global queries.fasta \
    --db db.fasta \
    --id 0.97 \
    --userout results.tsv \
    --userfields query+target+id+alnlen+mism
```

Output all five identity definitions side by side for comparison:

```sh
vsearch \
    --usearch_global queries.fasta \
    --db db.fasta \
    --id 0.0 \
    --userout identity_comparison.tsv \
    --userfields query+target+id0+id1+id2+id3+id4
```

Output query coverage and alignment span excluding terminal gaps:

```sh
vsearch \
    --usearch_global queries.fasta \
    --db db.fasta \
    --id 0.8 \
    --userout coverage.tsv \
    --userfields query+target+qcov+tcov+qilo+qihi+tilo+tihi
```


# SEE ALSO

[`vsearch-pairwise_alignment_parameters(7)`](./vsearch-pairwise_alignment_parameters.7.md),
[`vsearch-allpairs_global(1)`](../commands/vsearch-allpairs_global.1.md),
[`vsearch-cluster_fast(1)`](../commands/vsearch-cluster_fast.1.md),
[`vsearch-cluster_size(1)`](../commands/vsearch-cluster_size.1.md),
[`vsearch-cluster_smallmem(1)`](../commands/vsearch-cluster_smallmem.1.md),
[`vsearch-cluster_unoise(1)`](../commands/vsearch-cluster_unoise.1.md),
[`vsearch-search_exact(1)`](../commands/vsearch-search_exact.1.md),
[`vsearch-usearch_global(1)`](../commands/vsearch-usearch_global.1.md)


#(../commands/fragments/footer.md)
