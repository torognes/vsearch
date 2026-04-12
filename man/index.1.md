% vsearch(1) version 2.30.6 | vsearch manual
% Torbjørn Rognes, Tomás Flouri, and Frédéric Mahé
#(./commands/fragments/date.md)

# NAME

vsearch --- a versatile open-source tool for metabarcoding and metagenomics


# SYNOPSIS

| **vsearch** \<command\> \[_file_] \[_options_]

(see below for a list of all available commands)


# DESCRIPTION

vsearch is a versatile open-source tool for microbiome analysis,
including chimera detection, clustering, dereplication and
rereplication, extraction, FASTA/FASTQ/SFF file processing, masking,
orienting, pairwise alignment, restriction site cutting, searching,
shuffling, sorting, subsampling, and taxonomic classification of
amplicon sequences for metabarcoding, metagenomics, genomics, and
population genetics.

Each command is described in a dedicated manpage. For example, type
`man vsearch-usearch_global` to read about the `--usearch_global`
command. Command manpages belong to section 1 (executable programs).
Format and reference manpages belong to sections 5 and 7 respectively;
for example, type `man 5 vsearch-fastq` to read about the fastq format
as used by vsearch.


# VSEARCH COMMANDS

## General

**[`vsearch-help(1)`](./commands/vsearch-help.1.md)**
: List available commands and options.

**[`vsearch-version(1)`](./commands/vsearch-version.1.md)**
: Write version information, citation, and compression support status.


## Chimera detection

**[`vsearch-uchime_denovo(1)`](./commands/vsearch-uchime_denovo.1.md)**
: Detect chimeras *de novo* using the UCHIME algorithm.

**[`vsearch-uchime2_denovo(1)`](./commands/vsearch-uchime2_denovo.1.md)**
: Detect chimeras *de novo* using the UCHIME2 algorithm.

**[`vsearch-uchime3_denovo(1)`](./commands/vsearch-uchime3_denovo.1.md)**
: Detect chimeras *de novo* using the UCHIME2 algorithm with a stricter
  abundance skew.

**[`vsearch-uchime_ref(1)`](./commands/vsearch-uchime_ref.1.md)**
: Detect chimeras using a reference database.

**[`vsearch-chimeras_denovo(1)`](./commands/vsearch-chimeras_denovo.1.md)**
: Detect chimeras *de novo* in long exact sequences.


## Clustering

**[`vsearch-cluster_fast(1)`](./commands/vsearch-cluster_fast.1.md)**
: Clusterize sequences sorted by decreasing length.

**[`vsearch-cluster_size(1)`](./commands/vsearch-cluster_size.1.md)**
: Clusterize sequences sorted by decreasing abundance.

**[`vsearch-cluster_smallmem(1)`](./commands/vsearch-cluster_smallmem.1.md)**
: Clusterize pre-sorted sequences using minimal memory.

**[`vsearch-cluster_unoise(1)`](./commands/vsearch-cluster_unoise.1.md)**
: Denoise amplicon sequences using the UNOISE3 algorithm.


## Dereplication and rereplication

**[`vsearch-fastx_uniques(1)`](./commands/vsearch-fastx_uniques.1.md)**
: Merge strictly identical fasta or fastq sequences.

**[`vsearch-derep_fulllength(1)`](./commands/vsearch-derep_fulllength.1.md)**
: Merge strictly identical fasta or fastq sequences (fasta-only output).

**[`vsearch-derep_id(1)`](./commands/vsearch-derep_id.1.md)**
: Merge identical fasta or fastq sequences sharing the same label.

**[`vsearch-derep_prefix(1)`](./commands/vsearch-derep_prefix.1.md)**
: Merge fasta or fastq sequences with identical prefixes.

**[`vsearch-derep_smallmem(1)`](./commands/vsearch-derep_smallmem.1.md)**
: Merge strictly identical sequences using minimal memory.

**[`vsearch-rereplicate(1)`](./commands/vsearch-rereplicate.1.md)**
: Use abundance values to rereplicate fasta sequences.


## Extraction of sequences

**[`vsearch-fastx_getseq(1)`](./commands/vsearch-fastx_getseq.1.md)**
: Extract a sequence from a fasta or fastq file by label.

**[`vsearch-fastx_getseqs(1)`](./commands/vsearch-fastx_getseqs.1.md)**
: Extract sequences from a fasta or fastq file by label.

**[`vsearch-fastx_getsubseq(1)`](./commands/vsearch-fastx_getsubseq.1.md)**
: Extract a subsequence from a fasta or fastq file by label.


## FASTA/FASTQ/SFF file processing

**[`vsearch-fasta2fastq(1)`](./commands/vsearch-fasta2fastq.1.md)**
: Convert fasta entries into fastq entries with fake quality scores.

**[`vsearch-fastq_chars(1)`](./commands/vsearch-fastq_chars.1.md)**
: Analyze a fastq file to identify the quality encoding and range of
  quality score values used.

**[`vsearch-fastq_convert(1)`](./commands/vsearch-fastq_convert.1.md)**
: Convert between fastq encoding variants.

**[`vsearch-fastq_eestats(1)`](./commands/vsearch-fastq_eestats.1.md)**
: Report per-position quality and expected error statistics.

**[`vsearch-fastq_eestats2(1)`](./commands/vsearch-fastq_eestats2.1.md)**
: Report read retention across combinations of length and expected
  error cutoffs.

**[`vsearch-fastq_filter(1)`](./commands/vsearch-fastq_filter.1.md)**
: Trim and filter fastq sequences.

**[`vsearch-fastq_join(1)`](./commands/vsearch-fastq_join.1.md)**
: Join paired-end reads into one sequence with a gap.

**[`vsearch-fastq_mergepairs(1)`](./commands/vsearch-fastq_mergepairs.1.md)**
: Merge paired-end reads by aligning overlapping regions.

**[`vsearch-fastq_stats(1)`](./commands/vsearch-fastq_stats.1.md)**
: Analyze fastq sequences and output detailed statistics.

**[`vsearch-fastx_filter(1)`](./commands/vsearch-fastx_filter.1.md)**
: Trim and filter fasta or fastq sequences.

**[`vsearch-fastx_revcomp(1)`](./commands/vsearch-fastx_revcomp.1.md)**
: Reverse-complement fasta or fastq sequences.

**[`vsearch-sff_convert(1)`](./commands/vsearch-sff_convert.1.md)**
: Convert an SFF file to fastq.


## Masking

**[`vsearch-fastx_mask(1)`](./commands/vsearch-fastx_mask.1.md)**
: Mask low-complexity regions in fasta or fastq sequences.

**[`vsearch-maskfasta(1)`](./commands/vsearch-maskfasta.1.md)**
: Mask low-complexity regions in fasta sequences (deprecated;
  use `--fastx_mask`).


## Orienting

**[`vsearch-orient(1)`](./commands/vsearch-orient.1.md)**
: Use a reference database to orient fasta or fastq sequences.


## Pairwise alignment

**[`vsearch-allpairs_global(1)`](./commands/vsearch-allpairs_global.1.md)**
: Perform global pairwise alignments of all sequence pairs.


## Restriction site cutting

**[`vsearch-cut(1)`](./commands/vsearch-cut.1.md)**
: Use a restriction pattern to cut fasta sequences.


## Searching

**[`vsearch-search_exact(1)`](./commands/vsearch-search_exact.1.md)**
: Search for exact full-length matches against a database.

**[`vsearch-usearch_global(1)`](./commands/vsearch-usearch_global.1.md)**
: Search sequences against a reference database using global alignment.


## Shuffling and sorting

**[`vsearch-shuffle(1)`](./commands/vsearch-shuffle.1.md)**
: Randomize the order of fasta or fastq entries.

**[`vsearch-sortbylength(1)`](./commands/vsearch-sortbylength.1.md)**
: Sort fasta or fastq sequences by decreasing length.

**[`vsearch-sortbysize(1)`](./commands/vsearch-sortbysize.1.md)**
: Sort fasta or fastq sequences by decreasing abundance.


## Subsampling

**[`vsearch-fastx_subsample(1)`](./commands/vsearch-fastx_subsample.1.md)**
: Randomly subsample fasta or fastq sequences.


## Taxonomic classification

**[`vsearch-sintax(1)`](./commands/vsearch-sintax.1.md)**
: Classify sequences using the SINTAX algorithm.


## UDB database handling

**[`vsearch-makeudb_usearch(1)`](./commands/vsearch-makeudb_usearch.1.md)**
: Create a UDB database file from a fasta file.

**[`vsearch-udb2fasta(1)`](./commands/vsearch-udb2fasta.1.md)**
: Extract sequences from a UDB database file into a fasta file.

**[`vsearch-udbinfo(1)`](./commands/vsearch-udbinfo.1.md)**
: Display information about a UDB database file.

**[`vsearch-udbstats(1)`](./commands/vsearch-udbstats.1.md)**
: Report statistics about indexed words in a UDB database file.


# FILE FORMATS

**[`vsearch-fasta(5)`](./formats/vsearch-fasta.5.md)**
: The fasta format, as used by vsearch.

**[`vsearch-fastq(5)`](./formats/vsearch-fastq.5.md)**
: The fastq format, as used by vsearch.

**[`vsearch-sff(5)`](./formats/vsearch-sff.5.md)**
: The Standard Flowgram Format (SFF), used by Roche 454 and early Ion
  Torrent PGM sequencing platforms.

**[`vsearch-udb(5)`](./formats/vsearch-udb.5.md)**
: The UDB (USEARCH database) binary format, containing fasta sequences
  and a pre-computed k-mer index.


# REFERENCE PAGES

**[`vsearch-expected_error(7)`](./misc/vsearch-expected_error.7.md)**
: A quality summary metric for fastq sequences.

**[`vsearch-nucleotides(7)`](./misc/vsearch-nucleotides.7.md)**
: The IUPAC nucleotide symbols accepted by vsearch.

**[`vsearch-pairwise_alignment_parameters(7)`](./misc/vsearch-pairwise_alignment_parameters.7.md)**
: The pairwise alignment model implemented in vsearch.

**[`vsearch-userfields(7)`](./misc/vsearch-userfields.7.md)**
: The output fields available with the `--userout` option.


# SEE ALSO

[swarm](https://github.com/torognes/swarm),
[swipe](https://github.com/torognes/swipe),
[usearch](https://github.com/rcedgar/usearch12)


#(./commands/fragments/footer.md)


# VERSION HISTORY

(inject history here)
