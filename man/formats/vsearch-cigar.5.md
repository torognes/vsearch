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
query viewpoint. Hence, a cigar string encodes query modifications
needed to equal the target. However, beware that cigar strings in the
SAM format adopt the target viewpoint: they encode target
modifications needed to equal the query (see
[`vsearch-sam(5)`](./vsearch-sam.5.md)). vsearch emits CIGAR strings
in several search, clustering, and pairwise alignment outputs (see SEE
ALSO).

A CIGAR string is an ASCII string consisting of one or more
*operations*. Each operation is written as an optional non-negative
integer *run-length* immediately followed by a single uppercase letter
indicating the *operation type*:

```text
[run-length] operation [run-length] operation ...
```

For example, `3M2I3MD` describes an alignment that starts with three
matches or mismatches, followed by two insertions (residues present
only in the target), three matches or mismatches, and a single
deletion (a residue present only in the query).


## Operations

vsearch recognizes three operation types, corresponding to the three
possible column types in a gapped pairwise alignment (see
[`vsearch-pairwise_alignment_parameters(7)`](../misc/vsearch-pairwise_alignment_parameters.7.md)):

`M`
: Match or mismatch. The column aligns one query residue with one
  target residue. vsearch does not distinguish matches from mismatches
  in its CIGAR strings.

`I`
: Insertion. A residue is present in the target but absent from the
  query (a gap in the query); it must be inserted into the query to
  equal the target.

`D`
: Deletion. A residue is present in the query but absent from the
  target (a gap in the target); it must be deleted from the query to
  equal the target.

Other operation letters defined by the SAM specification (`X`, `=`,
`N`, `S`, `H`, `P`) are never produced by vsearch.


## Run-length encoding

Consecutive columns of the same operation type are grouped into a
single run. The number of columns in a run is written as a decimal
integer immediately before the operation letter, without leading
zeros. A run-length of 1 is implicit and omitted: vsearch writes `M`,
not `1M`, and `MID` rather than `1M1I1D`. Run-lengths are positive and
bounded by the alignment length; there is no maximum total length for
a CIGAR string.


## Empty and missing alignments

An empty CIGAR string is valid and represents an empty alignment.

When a command reports that a query has no alignment to a target, the
CIGAR field is replaced with a placeholder whose value depends on the
output format:

- the `aln` and `caln` [userfields](../misc/vsearch-userfields.7.md)
  (`--userout`) produce an empty field;
- the `--samout` (SAM) and `--uc` (UCLUST-like) outputs produce the
  single character `*`.

(The `--blast6out` tabular output contains no CIGAR field at all.)


## Exact-match shorthand

As a vsearch-specific convention, a CIGAR string consisting of the
single character `=` is emitted in the `--uc` output (and only there)
when the query and the target sequences are strictly identical; with
`--cluster_fast`, terminal gaps are ignored when deciding identity.
This shorthand is not part of the SAM specification; `--samout` and
the `caln` userfield always emit explicit operations.


## Ill-formed strings

vsearch never reads CIGAR strings from its input files: the strings
it writes are always well-formed as described above, and there is no
user-facing path that parses one back.


# EXAMPLES

The following alignment of an 8-nt query against a 9-nt target

```text
query:     ACGT--TACG
target:    AC-TGGTTCG
alignment: MMDMIIMMMM
cigar:     2M1D1M2I4M
```

contains, from left to right:

- 2 columns with matching residues (`AC` / `AC`),
- 1 deletion (`G` present only in the query, a gap in the target),
- 1 column with matching residues (`T` / `T`),
- 2 insertions (`GG` present only in the target, gaps in the query),
- 4 columns with matching residues (`TACG` / `TTCG`).

Without run-length encoding, the alignment reads `MMDMIIMMMM`, after
encoding it reads `2M1D1M2I4M`. Note that `M` does not distinguish
matches from mismatches; a mismatching column is encoded with the same
operation letter.


# SEE ALSO

[`vsearch-allpairs_global(1)`](../commands/vsearch-allpairs_global.1.md),
[`vsearch-cluster_fast(1)`](../commands/vsearch-cluster_fast.1.md),
[`vsearch-cluster_size(1)`](../commands/vsearch-cluster_size.1.md),
[`vsearch-cluster_smallmem(1)`](../commands/vsearch-cluster_smallmem.1.md),
[`vsearch-cluster_unoise(1)`](../commands/vsearch-cluster_unoise.1.md),
[`vsearch-search_exact(1)`](../commands/vsearch-search_exact.1.md),
[`vsearch-usearch_global(1)`](../commands/vsearch-usearch_global.1.md),
[`vsearch-sam(5)`](./vsearch-sam.5.md),
[`vsearch-pairwise_alignment_parameters(7)`](../misc/vsearch-pairwise_alignment_parameters.7.md),
[`vsearch-userfields(7)`](../misc/vsearch-userfields.7.md)


#(../commands/fragments/footer.md)
