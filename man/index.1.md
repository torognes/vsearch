% vsearch(1) version 2.30.0 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./commands/fragments/date.md)

# NAME

vsearch --- a versatile open-source tool for metabarcoding and metagenomics


# SYNOPSIS

| **vsearch** \<command\> \[_file_] \[_options_]

(see below for a list of available commands)


# DESCRIPTION

vsearch is a versatile open-source tool for microbiome analysis,
including chimera detection, clustering, dereplication and
rereplication, extraction, FASTA/FASTQ/SFF file processing, masking,
orienting, pairwise alignment, restriction site cutting, searching,
shuffling, sorting, subsampling, and taxonomic classification of
amplicon sequences for metabarcoding, metagenomics, genomics, and
population genetics.


# GIT COMMANDS

Chimera detection:

| **vsearch** (\-\-uchime_denovo | \-\-uchime2_denovo | \-\-uchime3_denovo) _fastafile_ (\-\-chimeras | \-\-nonchimeras | \-\-uchimealns | \-\-uchimeout) _outputfile_ \[_options_]
| **vsearch** \-\-uchime_ref _fastafile_ (\-\-chimeras | \-\-nonchimeras | \-\-uchimealns | \-\-uchimeout) _outputfile_ \-\-db _fastafile_ \[_options_]

Clustering:

| **vsearch** (\-\-cluster_fast | \-\-cluster_size | \-\-cluster_smallmem | \-\-cluster_unoise) _fastafile_ (\-\-alnout | \-\-biomout | \-\-blast6out | \-\-centroids | \-\-clusters | \-\-mothur_shared_out | \-\-msaout | \-\-otutabout | \-\-profile | \-\-samout | \-\-uc | \-\-userout) _outputfile_ \-\-id _real_ \[_options_]


# FILE FORMATS

**vsearch-fasta(5)**
: Define fasta format, as used by vsearch.

**vsearch-fastq(5)**
: Define fastq format, as used by vsearch.



# SEE ALSO

swarm, swipe, [`usearch`](https://github.com/rcedgar/usearch12)


#(../commands/fragments/footer.md)


# VERSION HISTORY

(inject history here)
