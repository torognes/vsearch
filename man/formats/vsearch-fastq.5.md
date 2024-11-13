% vsearch-fastq(5) version 2.30.0 | vsearch file formats
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(../commands/fragments/date.md)

# NAME

fastq --- a text-based format for representing nucleotide sequences
and their corresponding quality scores


# DESCRIPTION

(WIP)

In fastq files, each entry is made of _sequence header_ starting with
a symbol '@', a nucleotidic _sequence_ (same rules as for fasta
sequences), a _quality header_ starting with a symbol '+', and a
_quality string_ of ASCII characters (offset 33 or 64), each one
encoding the quality value of the corresponding position in the
nucleotidic sequence.

#(./fragments/sequences.md)

In fastq files, each entry is made of a _sequence header_, a
_sequence_, a _quality header_, ... The header is defined as the
string comprised between the initial '>' symbol and the first space,
tabulation, or new line symbol, unless the `--notrunclabels` option is
in effect, in which case the entire line is included.

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



# EXAMPLES

(give examples of valid and invalid fastq files)


# SEE ALSO

[`vsearch-fasta(5)`](./vsearch-fasta.5.md), [`vsearch-udb(5)`](./vsearch-udb.5.md)


#(../commands/fragments/footer.md)
