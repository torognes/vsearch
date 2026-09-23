% vsearch-pairwise_alignment_parameters(7) version 2.32.0 | vsearch file formats
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

The alignment score is used to find the optimal alignment. By
default, the similarity percentage reported (and compared against
`--id`) is computed from the columns of the resulting alignment, not
from the score itself. The identity definition can be changed with
`--iddef` (see below); `--iddef 5` computes it from the score.


## Matches and mismatches

When aligning sequences, identical symbols will receive a positive
match score (default +2, see option `--match`). Note that T and U are
considered identical, regardless of their case. If two symbols are not
identical, their alignment results in a negative mismatch score
(default -4, see option `--mismatch`). Aligning a pair of symbols
where at least one of them is an ambiguous symbol (BDHKMNRSVWY) will
always result in a score of zero. Alignment of two identical ambiguous
symbols (for example, R vs R) also receives a score of zero. As an
exception, when `--n_mismatch` is given, any pairing involving an `N`
is scored as a mismatch instead of zero.

Once the optimal pairwise alignment has been found, when computing the
amount of similarity by counting matches and mismatches **after**
alignment, ambiguous nucleotide symbols will count as matching to
other symbols if they have at least one of the nucleotides (ACGTU)
they may represent in common. For example: W will match A and T, but
also any of MRVHDN. When showing alignments (for example with the
output option `--alnout`) matches involving ambiguous symbols will be
shown with a plus character (+) between them while exact matches
between non-ambiguous symbols will be shown with a vertical bar
character (|). Here too, `--n_mismatch` is an exception: a column
holding an N is then counted as a mismatch rather than as a match, and
is shown with a blank character.


## Alignments involving Ns

Since an N stands for any of the four nucleotides, the default rules
above make every column holding an N a free match: it neither helps
nor hurts the score (zero), and it counts as a matching column when
the identity percentage is computed. A query aligned end-to-end over a
run of Ns is therefore reported at 100% identity, even though not a
single known nucleotide was compared. Reference databases assembled
from public repositories do contain such runs, and the resulting hits
are usually unwanted.

Note that the *k*-mer pre-filter is not what lets these hits through:
words containing an ambiguous symbol are never indexed, so a run of Ns
yields no word at all (a sequence made only of Ns cannot be selected as
a target, and vsearch warns about it). The pair is selected on words
shared elsewhere in the target, and it is the global alignment that
then slides the query into the run of Ns, because those columns cost
nothing while gaps do.

The option `--n_mismatch` changes that rule: any column where at least
one of the two symbols is an N, N against N included, is scored as a
mismatch and counted as a mismatch. The other ambiguous symbols
(BDHKMRSVWY) are not affected and keep the default behaviour described
above. There is currently no third option to exclude columns holding
an N from the identity computation altogether.

Ns are not always inherited from the input: `--hardmask` replaces
masked regions with Ns, so masked regions stop being compared and
become free matches. With the default (soft) masking, masked regions
keep their nucleotides, in lowercase, and are compared as usual;
`--n_mismatch` has no effect on them.


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
extremities (L and R) equally --- combining E with L or R in the same
penalty string is a fatal error --- and the symbols Q and T can be
omitted to treat query and target sequences equally.

`--gapopen` *20I/2E*
: Set the six gap opening penalties using a penalty of 20 for opening
  internal gaps and a penalty of 2 for opening terminal gaps, in both
  query and target sequences. This is the default.


If only a numerical value is given, without any sequence or location
symbol, then the penalty applies to all gap openings. For example:

`--gapopen` *20*
: Set the six gap opening penalties using a penalty of 20 for all gap
  openings, internal or terminal, in both query and target sequences.


To forbid a gap opening, an infinite penalty value can be declared with
the symbol `*`: a gap opening whose penalty is infinite is never
allowed. Like a numerical penalty, `*` applies to the selected sequence
(Q or T) and location (L, I, R or E); a bare `*` forbids every gap
opening (internal and terminal, in both sequences).

`--gapopen` *\*I/2E*
: Set the gap opening penalties to an infinite value for internal gap
  openings, in both query and target sequences.

`--gapopen` *\*LQ*
: Forbid the opening of left end-gaps in the query sequence; all other
  gap opening penalties keep their default values.


To use vsearch as a semi-global aligner, a null-penalty can be applied
to the left (L) or right (R) gaps.

vsearch always initializes the six gap opening penalties using the
default parameters (20I/2E). The user is then free to declare only the
values they want to modify.

The string is scanned from left to right, accepted symbols are
`0123456789/LIREQT*`, and later values override previous values.

Each finite penalty must be an integer between 0 and 6553. To declare a
larger, gap-forbidding penalty, use the infinite value `*` (see above).

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


An infinite extension penalty declared with `*` forbids extending a gap
of that class beyond a single position: gaps longer than one are not
allowed, while a gap of length one is still permitted. As with gap
openings, `*` applies to the selected sequence (Q or T) and location (L,
I, R or E).

`--gapext` *\*I*
: Forbid internal gaps longer than one position, in both query and
  target sequences.


### Gap costs

A gap of length *k* costs its opening penalty plus (*k* - 1) extension
penalties: the opening penalty covers the first position of the gap,
and each further position adds one extension penalty. With the default
penalties, an internal gap of length 1 costs 20, and an internal gap of
length 4 costs 20 + 3 × 2 = 26. With `--gapext 0I`, an internal gap
costs 20 whatever its length.

A gap is *terminal* when it is the first or the last operation of the
alignment; any other gap is internal and pays the internal penalties.
When two gaps follow each other at an end of the alignment, one in each
sequence (for example, a gap in the query immediately followed by a
gap in the target, after the last aligned pair), only the outermost of
the two is terminal.

The alignment score (userfield `raw`, see
[`vsearch-userfields(7)`](./vsearch-userfields.7.md)) is the sum of the
scores of the aligned pairs (`--match`, `--mismatch`, or zero for a
pair involving an ambiguous symbol) minus the costs of the gaps.


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

`--iddef` *5*
: Score-based definition: 1.0 - [(`--match` * *L* - score) /
  ((`--match` - `--mismatch`) * *L*)], where score is the alignment
  score and *L* the shortest sequence length, clamped to the range 0.0
  to 1.0. *L* * `--match` is the highest score the pair can reach; the
  amount by which the alignment falls short of it is counted in
  mismatch equivalents, one mismatch costing `--match` - `--mismatch`.
  Each gap therefore weighs its own penalty (see *Gap costs* above),
  and a column holding an ambiguous symbol, which scores zero, weighs
  `--match` / (`--match` - `--mismatch`) of a mismatch. Requires
  `--match` to be greater than `--mismatch`.

Note that the `--iddef` choice has no effect on the score or selection
of the optimal pairwise alignment. The identity is computed from the
alignment after the fact.

What counts as a matching column is the same for definitions 0 to 4,
and follows the rules given above: a column holding an ambiguous symbol
matches whenever the two symbols share at least one of the nucleotides
they represent, unless `--n_mismatch` is given. Definition 5 counts no
columns: it reads the score, in which such a column scores zero.


# SEE ALSO

[`vsearch-nucleotides(7)`](./vsearch-nucleotides.7.md),
[`vsearch-allpairs_global(1)`](../commands/vsearch-allpairs_global.1.md),
[`vsearch-cluster_fast(1)`](../commands/vsearch-cluster_fast.1.md),
[`vsearch-cluster_size(1)`](../commands/vsearch-cluster_size.1.md),
[`vsearch-cluster_smallmem(1)`](../commands/vsearch-cluster_smallmem.1.md),
[`vsearch-usearch_global(1)`](../commands/vsearch-usearch_global.1.md)


#(../commands/fragments/footer.md)
