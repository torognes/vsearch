% vsearch-pairwise_alignment_parameters(7) version 2.30.6 | vsearch file formats
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(../commands/fragments/date.md)

# NAME

pairwise alignment parameters --- a description of the pairwise
alignment model implemented in vsearch


# DESCRIPTION

vsearch implements an extremely fast Needleman-Wunsch algorithm, making
use of the Streaming SIMD Extensions (SSE2) of post-2003 x86-64 CPUs.
On Power8 CPUs it uses AltiVec/VSX/VMX instructions, and on ARMv8 CPUs
it uses Neon instructions. On other systems it can use the SIMD
Everywhere (simde) library, if available. For comparisons involving
sequences with a length product greater than 25 million (e.g. two
sequences of length 5 kb), vsearch uses a slower alignment method
described by Hirschberg (1975) and Myers and Miller (1988), with much
smaller memory requirements.

The Needleman-Wunsch algorithm performs *global* pairwise alignment: it
aligns two sequences end-to-end, over their full lengths. This differs
from Smith-Waterman, which performs *local* alignment and identifies the
highest-scoring matching subsequence. vsearch always uses global
alignment.

Along the pairwise alignment, each aligned position contributes to the
score, by being either a *match*, a *mismatch*, or a *gap*.

vsearch interprets symbols in DNA/RNA sequences according to the IUPAC
coding system for nucleotides. See
[`vsearch-nucleotides(7)`](./vsearch-nucleotides.7.md) for details.

The alignment score is used solely to find the optimal alignment. The
similarity percentage reported (and compared against `--id`) is
computed from the resulting alignment, not from the score itself. The
identity definition can be changed with `--iddef` (see below).


## Matches and mismatches

When aligning sequences, identical symbols will receive a positive
match score (default +2, see option `--match`). Note that T and U are
considered identical, regardless of their case. If two symbols are not
identical, their alignment results in a negative mismatch score
(default -4, see option `--mismatch`). Aligning a pair of symbols
where at least one of them is an ambiguous symbol (BDHKMNRSVWY) will
always result in a score of zero. Alignment of two identical ambiguous
symbols (for example, R vs R) also receives a score of zero.

Once the optimal pairwise alignment has been found, when computing the
amount of similarity by counting matches and mismatches **after**
alignment, ambiguous nucleotide symbols will count as matching to
other symbols if they have at least one of the nucleotides (ACGTU)
they may represent in common. For example: W will match A and T, but
also any of MRVHDN. When showing alignments (for example with the
output option `--alnout`) matches involving ambiguous symbols will be
shown with a plus character (+) between them while exact matches
between non-ambiguous symbols will be shown with a vertical bar
character (|).


## Gaps

Gaps are further refined into *gap openings* (see option `--gapopen`)
or *gap extensions* (see option `--gapext`). Gaps are asymmetrical: a
gap opening can occur in six different contexts: in the query (Q) or in
the target (T) sequence; inside the sequence (I), or at the left (L) or
right (R) extremity.


### Gap openings

Sequence symbols (Q and T) can be combined with location symbols (L, I,
and R), and numerical values to declare penalties for all possible
contexts: `aQL/bQI/cQR/dTL/eTI/fTR`, where *abcdef* are zero or
positive integers, and `/` is used as a separator.

`--gapopen` *2QL/20QI/2QR/2TL/20TI/2TR*
: Set the six gap opening penalties using a penalty of 20 for opening
  internal gaps and a penalty of 2 for opening terminal gaps, in both
  query and target sequences. This is the default.


To simplify declarations, the symbol (E) can be used to treat both
extremities (L and R) equally, and the symbols Q and T can be omitted
to treat query and target sequences equally.

`--gapopen` *20I/2E*
: Set the six gap opening penalties using a penalty of 20 for opening
  internal gaps and a penalty of 2 for opening terminal gaps, in both
  query and target sequences. This is the default.


If only a numerical value is given, without any sequence or location
symbol, then the penalty applies to all gap openings. For example:

`--gapopen` *20*
: Set the six gap opening penalties using a penalty of 20 for all gap
  openings, internal or terminal, in both query and target sequences.


To forbid gap-opening, an infinite penalty value can be declared with
the symbol `*`.

`--gapopen` *\*I/2E*
: Set the gap opening penalties to an infinite value for internal gap
  openings, in both query and target sequences.


To use vsearch as a semi-global aligner, a null-penalty can be applied
to the left (L) or right (R) gaps.

vsearch always initializes the six gap opening penalties using the
default parameters (20I/2E). The user is then free to declare only the
values they want to modify.

The string is scanned from left to right, accepted symbols are
`0123456789/LIREQT*`, and later values override previous values.

Please note that vsearch, in contrast to usearch, only allows integer
gap penalties. Because the lowest gap penalties are 0.5 by default in
usearch, all default scores and gap penalties in vsearch have been
doubled to maintain equivalent penalties and to produce identical
alignments.


### Gap extensions

Gap extensions follow the same penalty declaration system as gap
openings.

`--gapext` *2I/1E*
: Set the six gap extending penalties using a penalty of 2 for
  extending internal gaps and a penalty of 1 for extending terminal
  gaps, in both query and target sequences. This is the default.


## Identity definitions

The identity percentage computed from a pairwise alignment can be
defined in several ways. The option `--iddef` selects the definition
used when applying the `--id` threshold:

`--iddef` *0*
: CD-HIT definition: (matching columns) / (shortest sequence length).

`--iddef` *1*
: Edit distance: (matching columns) / (alignment length).

`--iddef` *2*
: Edit distance excluding terminal gaps: (matching columns) / (alignment
  length - terminal gaps). This is the default.

`--iddef` *3*
: Marine Biological Lab definition, counting each gap opening (internal
  or terminal) as a single mismatch, whether or not the gap was
  extended: 1.0 - [(mismatches + gap openings) / (longest sequence
  length)].

`--iddef` *4*
: BLAST definition, equivalent to `--iddef 1` for global pairwise
  alignments.

Note that the `--iddef` choice has no effect on the score or selection
of the optimal pairwise alignment. The identity is computed from the
alignment after the fact.


# SEE ALSO

[`vsearch-nucleotides(7)`](./vsearch-nucleotides.7.md),
[`vsearch-allpairs_global(1)`](../commands/vsearch-allpairs_global.1.md),
[`vsearch-cluster_fast(1)`](../commands/vsearch-cluster_fast.1.md),
[`vsearch-cluster_size(1)`](../commands/vsearch-cluster_size.1.md),
[`vsearch-cluster_smallmem(1)`](../commands/vsearch-cluster_smallmem.1.md),
[`vsearch-usearch_global(1)`](../commands/vsearch-usearch_global.1.md)


#(../commands/fragments/footer.md)
