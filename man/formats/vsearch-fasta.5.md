% vsearch-fasta(5) version 2.30.4 | vsearch file formats
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(../commands/fragments/date.md)

# NAME

fasta --- a text-based format for representing nucleotide sequences


# DESCRIPTION

In fasta files, each entry is made of a _header_ and a _sequence_. The
header is defined as the string comprised between the initial '>'
symbol and the first space, tabulation, or new line symbol, unless the
`--notrunclabels` option is in effect, in which case the entire line
is included.

The *header* should contain printable ASCII characters (values ranging
from 33 to 126, see `ascii(7)`). If unprintable ASCII characters are
encountered, vsearch terminates with a fatal error. If non-ASCII
characters are encountered (values ranging from 128 to 255), vsearch
emits a warning.

When using the option `--sizein`, if the header matches the patterns
`>[;]size=integer;identifier`, `>identifier;size=integer;identifier`,
or `>identifier;size=integer[;]`, vsearch will interpret *integer* as
the number of occurrences (or abundance) of the sequence in the
study. This abundance information is used, or created, during chimera
detection, clustering, dereplication, sorting and searching.


(TBD: describe sequences)


#(./fragments/format_sequence.md)


# SEE ALSO

[`vsearch-fastq(5)`](./vsearch-fastq.5.md), [`vsearch-udb(5)`](./vsearch-udb.5.md)


#(../commands/fragments/footer.md)
