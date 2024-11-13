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

The _header_ should contain printable ascii characters (33-126). The
program will terminate with a fatal error if there are unprintable
ascii characters (see `ascii(7)`). A warning will be issued if
non-ascii characters (128-255) are encountered.

If the header matches the pattern '>[;]size=integer;label', the
pattern '>label;size=integer;label', or the pattern
'>label;size=integer[;]', vsearch will interpret integer as the number
of occurrences (or abundance) of the sequence in the study. That
abundance information is used or created during chimera detection,
clustering, dereplication, sorting and searching.

The _sequence_ is defined as a string of IUPAC symbols
('ACGTURYSWKMDBHVN'), starting after the end of the header line and
ending before the next header line, or the file's end. vsearch
silently ignores ascii characters 9 to 13, and exits with an error
message if ascii characters 0 to 8, 14 to 31, '.' or '-' are
present. All other ascii or non-ascii characters are stripped and
complained about in a warning message.


# EXAMPLES

(give examples of valid and invalid fasta files)


# SEE ALSO

[`vsearch-fastq(5)`](./vsearch-fastq.5.md), [`vsearch-udb(5)`](./vsearch-udb.5.md)


#(../commands/fragments/footer.md)
