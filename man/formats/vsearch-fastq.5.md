% vsearch-fastq(5) version 2.30.0 | vsearch file formats
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(../commands/fragments/date.md)

# NAME

fastq --- a text-based format for representing nucleotide sequences
and their corresponding quality scores


# DESCRIPTION

In fastq files, each entry is made of *sequence header* starting with
a symbol '@', a nucleotidic *sequence*, a *quality header* starting
with a symbol '+', and a *quality string* of ASCII characters (offset
33 or 64), each one encoding the quality value of the corresponding
position in the nucleotidic sequence.

The *sequence header* is defined as the string comprised between the
initial '@' symbol and the first space, tabulation, or new line
symbol, unless the `--notrunclabels` option is in effect, in which
case the entire line is included.

The sequence header should contain printable ASCII characters
(33-126). The program will terminate with a fatal error if there are
unprintable ASCII characters (see `ascii(7)`). A warning will be
issued if non-ASCII characters (128-255) are encountered.

If the sequence header contains patterns such as `[@;]size=integer[;]`
or `[@;]ee=float[;]`, vsearch can interpret these annotations and use
them for chimera detection, clustering, dereplication, filtering and
sorting.

#(./fragments/format_sequence.md)

The *quality header* is defined as the string comprised between the
initial '+' symbol and the first space, tabulation, or new line
symbol.

The *quality string* is a string of ASCII characters, starting after
the end of the quality header line and ending before the next header
line, or the file's end. The range of valid ASCII characters can
extend from '!' to '~' when the offset is 33, and from '@' to '~' when
the offset is 64. vsearch silently ignores ASCII characters 9 to 13,
and exits with an error message if ASCII characters 0 to 8, 14 to 31,
‘.’ or ‘-’ are present. All other ASCII or non-ASCII characters are
stripped and complained about in a warning message.


# SEE ALSO

[`vsearch-fasta(5)`](./vsearch-fasta.5.md), [`vsearch-udb(5)`](./vsearch-udb.5.md)


#(../commands/fragments/footer.md)
