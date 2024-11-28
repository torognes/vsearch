% vsearch-fasta(5) version 2.30.0 | vsearch file formats
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

The *header* should contain printable ASCII characters (33-126). The
program will terminate with a fatal error if there are unprintable
ASCII characters (see `ascii(7)`). A warning will be issued if
non-ASCII characters (128-255) are encountered.

When using the option `--sizein`, if the header matches the patterns
`>[;]size=integer;identifier`, `>identifier;size=integer;identifier`,
or `>identifier;size=integer[;]`, vsearch will interpret *integer* as
the number of occurrences (or abundance) of the sequence in the
study. That abundance information is used or created during chimera
detection, clustering, dereplication, sorting and searching.

#(./fragments/format_sequence.md)


# SEE ALSO

[`vsearch-fastq(5)`](./vsearch-fastq.5.md), [`vsearch-udb(5)`](./vsearch-udb.5.md)


#(../commands/fragments/footer.md)
