% vsearch-nucleotides(7) version 2.33.0 | vsearch file formats
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(../commands/fragments/date.md)

# NAME

nucleotides --- a description of the nucleotide symbols accepted by
vsearch


# DESCRIPTION

vsearch interprets symbols in DNA/RNA sequences according to the IUPAC
coding system for nucleotides. This widely available table is
reproduced here in the form of a manpage for ease-of-use:

| Symbol | Description                   | Base represented | Complement |
|--------|-------------------------------|------------------|------------|
| A      | Adenine                       | A                | T          |
| C      | Cytosine                      | C                | G          |
| G      | Guanine                       | G                | C          |
| T      | Thymine                       | T                | A          |
| U      | Uracil                        | U                | A          |
|--------|-------------------------------|------------------|------------|
| B      | not A (B comes after A)       | C or G or T      | V          |
| D      | not C (D comes after C)       | A or G or T      | H          |
| H      | not G (H comes after G)       | A or C or T      | D          |
| K      | bases that are ketones        | G or T           | M          |
| M      | bases with amino groups       | A or C           | K          |
| N      | Nucleic acid (any base)       | A or C or G or T | N          |
| R      | purine                        | A or G           | Y          |
| S      | Strong interaction            | C or G           | S          |
| V      | not T (V comes after T and U) | A or C or G      | B          |
| W      | Weak interaction              | A or T           | W          |
| Y      | pyrimidine                    | C or T           | R          |
| -      | Gap                           |                  |            |

Note that the symbols 'X' (*Masked*) and 'I' (*Inosine*) are **not**
accepted by vsearch. Neither belongs to the IUPAC set above, and both
are stripped from input sequences, with a warning. Inosine is left out
deliberately, and not for want of a mapping: it pairs with all four
bases, most stably with C, so 'I' is read by some tools as a synonym
of 'G' (inosine being guanine without its 2-amino group), and by
others as a universal base equivalent to 'N'. Rather than silently
pick one of the two, vsearch strips the symbol and reports it.

The gap symbol '-' is listed for completeness only: vsearch *writes*
it in alignment outputs, but never accepts it in input sequences
(a fatal error in both fasta and fastq files).

vsearch reads nucleotide sequences only; amino acid sequences are not
supported. A protein input is nevertheless not rejected, because many
one-letter amino acid codes are also valid nucleotide symbols in the
table above: the residues that are not are stripped with a warning,
and the remainder is read as a much shorter nucleotide sequence. A
60-residue protein can thus become a 35-symbol sequence, silently. The
stripping warning carries a reminder to that effect; results computed
from such an input are meaningless.

How these symbols behave when two sequences are compared is a separate
matter, described in
[`vsearch-pairwise_alignment_parameters(7)`](./vsearch-pairwise_alignment_parameters.7.md):
in short, a column holding an ambiguous symbol scores zero, and counts
as a matching column whenever the two symbols share at least one of
the nucleotides they represent. An N therefore matches anything, which
is why a query can align to a run of Ns with 100% identity; the option
`--n_mismatch` counts such columns as mismatches instead.

# SEE ALSO

[`vsearch-pairwise_alignment_parameters(7)`](./vsearch-pairwise_alignment_parameters.7.md),
[`vsearch-fasta(5)`](../formats/vsearch-fasta.5.md),
[`vsearch-fastq(5)`](../formats/vsearch-fastq.5.md)


#(../commands/fragments/footer.md)
