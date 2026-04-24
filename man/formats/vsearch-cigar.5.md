% vsearch-cigar(5) version 2.31.0 | vsearch file formats
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(../commands/fragments/date.md)

# NAME

cigar --- a compact text-based format for representing pairwise alignments


# DESCRIPTION

CIGAR (Compact Idiosyncratic Gapped Alignment Report) is a text-based
format that encodes a pairwise alignment as a sequence of run-length
encoded operations. Each operation describes a column of the
alignment, or a run of identical consecutive columns, as seen from the
query viewpoint. vsearch emits CIGAR strings in several search,
clustering, and pairwise alignment outputs (see SEE ALSO).

A CIGAR string is an ASCII string consisting of one or more
*operations*. Each operation is written as an optional non-negative
integer *run-length* immediately followed by a single uppercase letter
indicating the *operation type*:

```text
[run-length] operation [run-length] operation ...
```

For example, `3M2I3MD` describes an alignment that starts with three
matches or mismatches, followed by two query insertions, three matches
or mismatches, and a single deletion.


## Operations

vsearch recognizes three operation types, corresponding to the three
possible column types in a gapped pairwise alignment (see
[`vsearch-pairwise_alignment_parameters(7)`](../misc/vsearch-pairwise_alignment_parameters.7.md)):

`M`
: Match or mismatch. The column aligns one query residue with one
  target residue. vsearch does not distinguish matches from mismatches
  in its CIGAR strings.

`I`
: Insertion. A residue is present in the query but absent from the
  target (a gap in the target).

`D`
: Deletion. A residue is present in the target but absent from the
  query (a gap in the query).

Other operation letters defined by the SAM specification (`X`, `=`,
`N`, `S`, `H`, `P`) are not produced by vsearch and are rejected when
reading CIGAR strings.


## Run-length encoding

Consecutive columns of the same operation type are grouped into a
single run. The number of columns in a run is written as a decimal
integer immediately before the operation letter. A run-length of 1 is
implicit and may be omitted: `M` is equivalent to `1M`, and `MID` is
equivalent to `1M1I1D`. Leading zeros are accepted (`03M` is
equivalent to `3M`). A run-length of 0 is accepted and produces a
zero-length operation.

Run-lengths are positive integers at most equal to the largest value
representable by the C `int` type (2,147,483,647 on most platforms,
see `limits.h(0p)`). There is no maximum total length for a CIGAR
string.


## Empty and missing alignments

An empty CIGAR string is valid and represents an empty alignment.

When a command reports that a query has no alignment to a target, the
CIGAR field is replaced with a placeholder whose value depends on the
output format:

- the `aln` and `caln` [userfields](../misc/vsearch-userfields.7.md)
  (`--userout`) produce an empty field;
- the `--blast6out` tabular output produces an empty field;
- the `--samout` (SAM) and `--uc` (UCLUST-like) outputs produce the
  single character `*`.


## Exact-match shorthand

As a vsearch-specific convention, a CIGAR string consisting of the
single character `=` is emitted in `--uc` and `--userout`
(`--userfields caln`) outputs when the query is identical to the
target, ignoring terminal gaps. This shorthand is not part of the SAM
specification and is not produced in `--samout` or `--blast6out`
outputs, which always emit explicit operations.


## Ill-formed strings

vsearch reports a fatal error when reading a CIGAR string that ends
with one or more digits not followed by an operation letter (for
example `12M1`), or that contains any character other than the digits
`0`--`9` and the letters `M`, `I`, and `D`.


# EXAMPLES

The following alignment of an 8-nt query against a 9-nt target

```text
query:   ACGT--TACG
target:  AC-TGGTACG
cigar:   2M1I1M2D4M
```

contains, from left to right:

- 2 columns with matching residues (`AC` / `AC`),
- 1 insertion in the target (`G` present only in the query),
- 1 column with matching residues (`T` / `T`),
- 2 deletions from the query (`GG` present only in the target),
- 4 columns with matching residues (`TACG` / `TACG`).

The same alignment written without run-length compression would be
`MMIMDDMMMM`. Note that `M` does not distinguish matches from
mismatches; a mismatching column is encoded with the same operation
letter.


# SEE ALSO

[`vsearch-allpairs_global(1)`](../commands/vsearch-allpairs_global.1.md),
[`vsearch-cluster_fast(1)`](../commands/vsearch-cluster_fast.1.md),
[`vsearch-cluster_size(1)`](../commands/vsearch-cluster_size.1.md),
[`vsearch-cluster_smallmem(1)`](../commands/vsearch-cluster_smallmem.1.md),
[`vsearch-cluster_unoise(1)`](../commands/vsearch-cluster_unoise.1.md),
[`vsearch-search_exact(1)`](../commands/vsearch-search_exact.1.md),
[`vsearch-usearch_global(1)`](../commands/vsearch-usearch_global.1.md),
[`vsearch-pairwise_alignment_parameters(7)`](../misc/vsearch-pairwise_alignment_parameters.7.md),
[`vsearch-userfields(7)`](../misc/vsearch-userfields.7.md)


#(../commands/fragments/footer.md)
