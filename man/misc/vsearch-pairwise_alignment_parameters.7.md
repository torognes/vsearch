% vsearch-pairwise_alignment_parameters(7) version 2.30.0 | vsearch file formats
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(../commands/fragments/date.md)

# NAME

pairwise alignment parameters --- a description of the pairwise
alignment model implemented in vsearch


# DESCRIPTION

vsearch implements an extremely fast Needleman-Wunsch algorithm
(REFERENCE). There are many ways to align two sequences. This
algorithm finds the pairwise alignment with the best score.

(explain the difference with Smith-Waterman?)

Along the pairwise alignment, each aligned position contributes to the
score, by being either a *match*, a *mismatch*, or a *gap*. 

## match and mismatches

When aligning sequences, identical symbols will receive a positive
match score (default +2, see option `--match`). Note that T and U are
considered identical, regardless of their case. If two symbols are not
identical, their alignment result in a negative mismatch score
(default -4, see option `--mismatch`). Aligning a pair of symbols
where at least one of them is an ambiguous symbol (BDHKMNRSVWY) will
always result in a score of zero. Alignment of two identical ambiguous
symbols (for example, R vs R) also receives a score of zero.

Once the optimal pairwise alignent has been found: When computing the
amount of similarity by counting matches and mismatches **after**
alignment, ambiguous nucleotide symbols will count as matching to
other symbols if they have at least one of the nucleotides (ACGTU)
they may represent in common. For example: W will match A and T, but
also any of MRVHDN. When showing alignments (for example with the
output option `--alnout`) matches involving ambiguous symbols will be
shown with a plus character (+) between them while exact matches
between non-ambiguous symbols will be shown with a vertical bar
character (|).

Note that --iddef choice has no effect on the score and selection of
the optimal pairwise alignment. The similarity is computed after.


## gaps

Gaps are further refined into *gap openings* (see option `--gapopen`)
or *gap extensions* (see option `--gapext`).


### gap openings

Contrary to matches/mismatches, gaps are asymetrical, so a gap opening
can occur in six different contexts: in the query (Q) or in the target
(T) sequence; inside the sequence (I), or at the left (L) or right (R)
extremity. Sequence symbols (Q and T) can be combined with location
symbols (L, I, and R), and numerical values to declare penalties for
all possible contexts: aQL/bQI/cQR/dTL/eTI/fTR, where abcdef are zero
or positive integers, and '/' is used as a separator.

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
openings, internal or terminal, in both query and target
sequences.


To forbid gap-opening, an infinite penalty value can be declared with
the symbol '\*'.

`--gapopen` *\*I/2E*
: Set the gap opening penalties to an infinite value for internal gap
openings, in both query and target sequences.


To use vsearch as a semi-global aligner, a null-penalty can be applied
to the left (L) or right (R) gaps.

vsearch always initializes the six gap opening penalties using the
default parameters (20I/2E). The user is then free to declare only the
values he/she wants to modify.

The string is scanned from left to right, accepted symbols are
(0123456789*/LIREQT), and later values override previous values.

Please note that vsearch, in contrast to usearch, only allows integer
gap penalties. Because the lowest gap penalties are 0.5 by default in
usearch, all default scores and gap penalties in vsearch have been
doubled to maintain equivalent penalties and to produce identical
alignments.


### gap extensions

Gap extensions follow the same penalty declaration system as gap
openings.

`--gapext` *2I/1E*
: Set the six gap extending penalties using a penalty of 2 for
extending internal gaps and a penalty of 1 for extending terminal
gaps, in both query and target sequences. This is the default.


# EXAMPLES

(is this section needed?)


# SEE ALSO

(nothing for now)


#(../commands/fragments/footer.md)
