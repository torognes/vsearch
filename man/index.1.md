% vsearch(1) version 2.30.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./commands/fragments/date.md)

# NAME

vsearch --- a versatile open-source tool for metabarcoding and metagenomics


# SYNOPSIS

| **vsearch** \<command\> \[_file_] \[_options_]

(see below for a list of available commands)

WARNING: this is a pre-release of the future revised and updated
vsearch documentation. Users can already refer to this content for the
commands listed below, but be aware that most vsearch commands are not
yet featured.


# DESCRIPTION

vsearch is a versatile open-source tool for microbiome analysis,
including chimera detection, clustering, dereplication and
rereplication, extraction, FASTA/FASTQ/SFF file processing, masking,
orienting, pairwise alignment, restriction site cutting, searching,
shuffling, sorting, subsampling, and taxonomic classification of
amplicon sequences for metabarcoding, metagenomics, genomics, and
population genetics.


# VSEARCH COMMANDS

Each command is described in a dedicated manpage. For example, type
`man vsearch-fasta2fastq` to read about the `--fasta2fastq`
command. Manpages reffering to commands all belong to section 1
(executable programs). Some vsearch manpages belong to section 5 (file
formats and conventions). For example, type `man 'vsearch-fastq(5)'`
to read about the fastq format as-seen by vsearch.

<!---
## Chimera detection:

| **vsearch** (\-\-uchime_denovo | \-\-uchime2_denovo | \-\-uchime3_denovo) _fastafile_ (\-\-chimeras | \-\-nonchimeras | \-\-uchimealns | \-\-uchimeout) _outputfile_ \[_options_]
| **vsearch** \-\-uchime_ref _fastafile_ (\-\-chimeras | \-\-nonchimeras | \-\-uchimealns | \-\-uchimeout) _outputfile_ \-\-db _fastafile_ \[_options_]

## Clustering:

| **vsearch** (\-\-cluster_fast | \-\-cluster_size | \-\-cluster_smallmem | \-\-cluster_unoise) _fastafile_ (\-\-alnout | \-\-biomout | \-\-blast6out | \-\-centroids | \-\-clusters | \-\-mothur_shared_out | \-\-msaout | \-\-otutabout | \-\-profile | \-\-samout | \-\-uc | \-\-userout) _outputfile_ \-\-id _real_ \[_options_]

--->

## FASTA/FASTQ/SFF file processing:

**[`vsearch-fasta2fastq(1)`](./commands/vsearch-fasta2fastq.1.md)**
: Convert a fasta file to a fastq file with fake quality scores.

**[`vsearch-fastq_chars(1)`](./commands/vsearch-fastq_chars.1.md)**
: Analyze fastq files to identify the quality encoding and the range
of quality score values used.

<!---
## Orienting:

**[`vsearch-orient(1)`](./commands/vsearch-orient.1.md)**
: Use a reference database to orient fastq or fasta sequences.

-->

# FILE FORMATS

**[`vsearch-fasta(5)`](./formats/vsearch-fasta.5.md)**
: Specify the fasta format, as used by vsearch.

**[`vsearch-fastq(5)`](./formats/vsearch-fastq.5.md)**
: Specify the fastq format, as used by vsearch.

(also SFF, and UDB).


# SEE ALSO

[swarm](https://github.com/torognes/swarm),
[swipe](https://github.com/torognes/swipe),
[usearch](https://github.com/rcedgar/usearch12)


#(./commands/fragments/footer.md)


# VERSION HISTORY

(inject history here)
