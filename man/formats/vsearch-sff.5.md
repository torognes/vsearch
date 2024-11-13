% vsearch-sff(5) version 2.30.0 | vsearch file formats
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(../commands/fragments/date.md)

# NAME

sff --- a binary file format used to encode pyrosequencing results


# DESCRIPTION

Standard flowgram format (sff), used by Roche 454 and Ion Torrent

(complete description of the SFF format, hard to find on the web)

A sff file must start with a 31-byte common header:
- magic value: a 32-bit value spelling '.sff' (0x2e736666),
- version value: a 32-bit value set to 1 (0x00000001)
- index offset: a 64-bit value
- index length: a 32-bit value
- number of reads: a 32-bit value
- header length: a 16-bit value
- key length: a 16-bit value
- flow length: a 16-bit value
- flowgram format: a 8-bit value

sff files are in big endian notation (no known counter-example)


# EXAMPLES

(show how to build sff files?)


# SEE ALSO

[`vsearch-fasta(5)`](./vsearch-fasta.5.md), [`vsearch-fastq(5)`](./vsearch-fastq.5.md)


#(../commands/fragments/footer.md)
