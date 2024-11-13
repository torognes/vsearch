% vsearch-sff(5) version 2.30.0 | vsearch file formats
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(../commands/fragments/date.md)

# NAME

sff --- a binary file format used to encode pyrosequencing results


# DESCRIPTION

used by Roche 454 and ionTorent

(complete description of the SFF format, hard to find on the web)

A sff file must start with a 31-byte common header:
- magic value: 4 bytes spelling '.sff' (0x2e736666),
- version value: a 32-bit value of 1 (0x00000001)



# EXAMPLES

(show how to build sff files?)


# SEE ALSO

[`vsearch-fasta(5)`](./vsearch-fasta.5.md), [`vsearch-fastq(5)`](./vsearch-fastq.5.md)


#(../commands/fragments/footer.md)
