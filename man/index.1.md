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

**[`vsearch-cigar(5)`](./formats/vsearch-cigar.5.md)**
: The CIGAR (Compact Idiosyncratic Gapped Alignment Report) format,
  used by vsearch to encode pairwise alignments.

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

New features and important modifications of **vsearch** (short-lived
or minor bug releases may not be mentioned):

**v1.0.0** released November 28th, 2014
:   First public release.

**v1.0.1** released December 1st, 2014
:   Bug fixes (sortbysize, semicolon after size annotation in headers)
    and minor changes (labels as secondary sort key for most sorts,
    treat T and U as identical for dereplication, only output size in
    `--dbmatched` file if `--sizeout` specified).

**v1.0.2** released December 6th, 2014
:   Bug fixes (ssse3/sse4.1 requirement, memory leak).

**v1.0.3** released December 6th, 2014
:   Bug fix (now writes help to stdout instead of stderr).

**v1.0.4** released December 8th, 2014
:   Added `--allpairs_global` option. Reduce memory requirements
    slightly and eliminate memory leaks.

**v1.0.5** released December 9th, 2014
:   Fixes a minor bug with `--allpairs_global` and `--acceptall`
    options.

**v1.0.6** released December 14th, 2014
:   Fixes a memory allocation bug in chimera detection (`--uchime_ref`
    option).

**v1.0.7** released December 19th, 2014
:   Fixes a bug in the output from chimera detection with the
    `--uchimeout` option.

**v1.0.8** released January 22nd, 2015
:   Introduces several changes and bug fixes:

    - a new linear memory aligner for alignment of sequences longer
      than 5,000 nucleotides,
    - a new `--cluster_size` command that sorts sequences by decreasing
      abundance before clustering,
    - meaning of userfields qlo, qhi, tlo, thi changed for
      compatibility with usearch,
    - new userfields qilo, qihi, tilo, tihi give alignment coordinates
      ignoring terminal gaps,
    - in `--uc` output files, a perfect alignment is indicated with a
      `=` sign,
    - the option `--cluster_fast` now sorts sequences by decreasing
      length, then by decreasing abundance and finally by sequence
      identifier,
    - default `--maxseqlength` value set to 50,000 nucleotides,
    - fix for bug in alignment in rare cases,
    - fix for lack of detection of under- or overflow in SIMD aligner.

**v1.0.9** released January 22nd, 2015
:   Fixes a bug in the function sorting sequences by decreasing
    abundance (`--sortbysize`).

**v1.0.10** released January 23rd, 2015
:   Fixes a bug where the `--sizein` option was ignored and always
    treated as on, affecting clustering and dereplication commands.

**v1.0.11** released February 5th, 2015
:   Introduces the possibility to output results in SAM format (for
    clustering, pairwise alignment and searching).

**v1.0.12** released February 6th, 2015
:   Temporarily fixes a problem with long headers in FASTA files.

**v1.0.13** released February 17th, 2015
:   Fix a memory allocation problem when computing multiple sequence
    alignments with the `--msaout` and `--consout` options, as well as
    a memory leak. Also increased line buffer for reading FASTA files
    to 4MB.

**v1.0.14** released February 17th, 2015
:   Fix a bug where the multiple alignment and consensus sequence
    computed after clustering ignored the strand of the sequences. Also
    decreased size of line buffer for reading FASTA files to 1MB again
    due to excessive stack memory usage.

**v1.0.15** released February 18th, 2015
:   Fix bug in calculation of identity metric between sequences when
    using the MBL definition (`--iddef 3`).

**v1.0.16** released February 19th, 2015
:   Integrated patches from Debian for increased compatibility with
    various architectures.

**v1.1.0** released February 20th, 2015
:   Added the `--quiet` option to suppress all output to stdout and
    stderr except for warnings and fatal errors. Added the `--log`
    option to write messages to a log file.

**v1.1.1** released February 20th, 2015
:   Added info about `--log` and `--quiet` options to help text.

**v1.1.2** released March 18th, 2015
:   Fix bug with large datasets. Fix format of help info.

**v1.1.3** released March 18th, 2015
:   Fix more bugs with large datasets.

**v1.2.0--v1.2.19** released July 6th to September 8th, 2015
:   Several new commands and options added. Bugs fixed. Documentation
    updated.

**v1.3.0** released September 9th, 2015
:   Changed to autotools build system.

**v1.3.1** released September 14th, 2015
:   Several new commands and options. Bug fixes.

**v1.3.2** released September 15th, 2015
:   Fixed memory leaks. Added `-h` shortcut for help. Removed extra
    `v` in version number.

**v1.3.3** released September 15th, 2015
:   Fixed bug in hexadecimal digits of MD5 and SHA1 digests. Added
    `--samheader` option.

**v1.3.4** released September 16th, 2015
:   Fixed compilation problems with zlib and bzip2lib.

**v1.3.5** released September 17th, 2015
:   Minor configuration/makefile changes to compile to native CPU and
    simplify makefile.

**v1.4.0** released September 25th, 2015
:   Added `--sizeorder` option.

**v1.4.1** released September 29th, 2015
:   Inserted public domain MD5 and SHA1 code to eliminate dependency
    on crypto and openssl libraries and their licensing issues.

**v1.4.2** released October 2nd, 2015
:   Dynamic loading of libraries for reading gzip and bzip2 compressed
    files if available. Circumvention of missing gzoffset function in
    zlib 1.2.3 and earlier.

**v1.4.3** released October 3rd, 2015
:   Fix a bug with determining amount of memory on some versions of
    Apple OS X.

**v1.4.4** released October 3rd, 2015
:   Remove debug message.

**v1.4.5** released October 6th, 2015
:   Fix memory allocation bug when reading long FASTA sequences.

**v1.4.6** released October 6th, 2015
:   Fix subtle bug in SIMD alignment code that reduced accuracy.

**v1.4.7** released October 7th, 2015
:   Fixes a problem with searching for or clustering sequences with
    repeats. In this new version, vsearch looks at all words occurring
    at least once in the sequences in the initial step. Previously only
    words occurring exactly once were considered. In addition, vsearch
    now requires at least 10 words to be shared by the sequences,
    previously only 6 were required. If the query contains less than 10
    words, all words must be present for a match. This change seems to
    lead to slightly reduced recall, but somewhat increased precision,
    ending up with slightly improved overall accuracy.

**v1.5.0** released October 7th, 2015
:   Introduces the new option `--minwordmatches` that allows the user
    to specify the minimum number of matching unique words before a
    sequence is considered further. New default values for different
    word lengths are also set. The minimum word length is increased
    to 7.

**v1.6.0** released October 9th, 2015
:   Adds the relabeling options (`--relabel`, `--relabel_md5` and
    `--relabel_sha1`) to the shuffle command. Also adds the `--xsize`
    option to the clustering, dereplication, shuffling and sorting
    commands.

**v1.6.1** released October 14th, 2015
:   Fix bugs and update manual and help text regarding relabelling. Add
    all relabelling options to the subsampling command. Add the
    `--xsize` option to chimera detection, dereplication and fastq
    filtering commands. Refactoring of code.

**v1.7.0** released October 14th, 2015
:   Add `--relabel_keep` option.

**v1.8.0** released October 19th, 2015
:   Added `--search_exact`, `--fastx_mask` and `--fastq_convert`
    commands. Changed most commands to read FASTQ input files as well
    as FASTA files. Modified `--fastx_revcomp` and `--fastx_subsample`
    to write FASTQ files.

**v1.8.1** released November 2nd, 2015
:   Fixes for compatibility with QIIME and older OS X versions.

**v1.9.0** released November 12th, 2015
:   Added the `--fastq_mergepairs` command and associated options. This
    command has not been tested well yet. Included additional files to
    avoid dependency of autoconf for compilation. Fixed an error where
    identifiers in fasta headers were not truncated at tabs, just
    spaces. Fixed a bug in detection of the file format (FASTA/FASTQ)
    of a gzip compressed input file.

**v1.9.1** released November 13th, 2015
:   Fixed memory leak and a bug in score computation in
    `--fastq_mergepairs`, and improved speed.

**v1.9.2** released November 17th, 2015
:   Fixed a bug in the computation of some values with
    `--fastq_stats`.

**v1.9.3** released November 19th, 2015
:   Workaround for missing x86intrin.h with old compilers.

**v1.9.4** released December 3rd, 2015
:   Fixed incrementation of counter when relabeling dereplicated
    sequences.

**v1.9.5** released December 3rd, 2015
:   Fixed bug resulting in inferior chimera detection performance.

**v1.9.6** released January 8th, 2016
:   Fixed bug in aligned sequences produced with `--fastapairs` and
    `--userout` (qrow, trow) options.

**v1.9.7** released January 12th, 2016
:   Masking behaviour is changed somewhat to keep the letter case of
    the input sequences unchanged when no masking is performed. Masking
    is now performed also during chimera detection. Documentation
    updated.

**v1.9.8** released January 22nd, 2016
:   Fixed bug causing segfault when chimera detection is performed on
    extremely short sequences.

**v1.9.9** released January 22nd, 2016
:   Adjusted default minimum number of word matches during searches for
    improved performance.

**v1.9.10** released January 25th, 2016
:   Fixed bug related to masking and lower case database sequences.

**v1.10.0** released February 11th, 2016
:   Parallelized and improved merging of paired-end reads and adjusted
    some defaults. Removed progress indicator when stderr is not a
    terminal. Added `--fasta_score` option to report chimera scores in
    FASTA files. Added `--rereplicate` and `--fastq_eestats` commands.
    Fixed typos. Added relabelling to files produced with `--consout`
    and `--profile` options.

**v1.10.1** released February 23rd, 2016
:   Fixed a bug affecting the `--fastq_mergepairs` command causing
    FASTQ headers to be truncated at first space (despite the bug fix
    release 1.9.0 of November 12th, 2015). Full headers are now
    included in the output (no matter if `--notrunclabels` is in effect
    or not).

**v1.10.2** released March 18th, 2016
:   Fixed a bug causing a segmentation fault when running
    `--usearch_global` with an empty query sequence. Also fixed a bug
    causing imperfect alignments to be reported with an alignment
    string of `=` in uc output files. Fixed typos in man file. Fixed
    fasta/fastq processing code regarding presence or absence of
    compression library header files.

**v1.11.1** released April 13th, 2016
:   Added strand information in UC file for `--derep_fulllength` and
    `--derep_prefix`. Added expected errors (ee) to header of FASTA
    files specified with `--fastaout` and `--fastaout_discarded` when
    `--eeout` or `--fastq_eeout` option is in effect for fastq_filter
    and fastq_mergepairs. The options `--eeout` and `--fastq_eeout` are
    now equivalent.

**v1.11.2** released June 21st, 2016
:   Two bugs were fixed. The first issue was related to the
    `--query_cov` option that used a different coverage definition than
    the qcov userfield. The coverage is now defined as the fraction of
    the whole query sequence length that is aligned with matching or
    mismatching residues in the target. All gaps are ignored. The other
    issue was related to the consensus sequences produced during
    clustering when only N's were present in some positions. Previously
    these would be converted to A's in the consensus. The behaviour is
    changed so that N's are produced in the consensus, and it should
    now be more compatible with usearch.

**v2.0.0** released June 24th, 2016
:   This major new version supports reading from pipes. Two new options
    are added: `--gzip_decompress` and `--bzip2_decompress`. One of
    these options must be specified if reading compressed input from a
    pipe, but are not required when reading from ordinary files. The
    vsearch header that was previously written to stdout is now written
    to stderr. This enables piping of results for further processing.
    The file name `-` now represents standard input (/dev/stdin) or
    standard output (/dev/stdout) when reading or writing files,
    respectively. Code for reading FASTA and FASTQ files has been
    refactored.

**v2.0.1** released June 30th, 2016
:   Avoid segmentation fault when masking very long sequences.

**v2.0.2** released July 5th, 2016
:   Avoid warnings when compiling with GCC 6.

**v2.0.3** released August 2nd, 2016
:   Fixed bad compiler options resulting in Illegal instruction errors
    when running precompiled binaries.

**v2.0.4** released September 1st, 2016
:   Improved error message for bad FASTQ quality values. Improved
    manual.

**v2.0.5** released September 9th, 2016
:   Add options `--fastaout_discarded` and `--fastqout_discarded` to
    output discarded sequences from subsampling to separate files.
    Updated manual.

**v2.1.0** released September 16th, 2016
:   New command: `--fastx_filter`. New options: `--fastq_maxlen`,
    `--fastq_truncee`. Allow `--minwordmatches` down to 3.

**v2.1.1** released September 23rd, 2016
:   Fixed bugs in output to UC-files. Improved help text and manual.

**v2.1.2** released September 28th, 2016
:   Fixed incorrect abundance output from fastx_filter and fastq_filter
    when relabelling.

**v2.2.0** released October 7th, 2016
:   Added OTU table generation options `--biomout`,
    `--mothur_shared_out` and `--otutabout` to the clustering and
    searching commands.

**v2.3.0** released October 10th, 2016
:   Allowed zero-length sequences in FASTA and FASTQ files. Added
    `--fastq_trunclen_keep` option. Fixed bug with output of OTU tables
    to pipes.

**v2.3.1** released November 16th, 2016
:   Fixed bug where `--minwordmatches 0` was interpreted as the default
    minimum word matches for the given word length instead of zero.
    When used in combination with `--maxaccepts 0` and `--maxrejects 0`
    it will allow complete bypass of kmer-based heuristics.

**v2.3.2** released November 18th, 2016
:   Fixed bug where vsearch reported the ordinal number of the target
    sequence instead of the cluster number in column 2 on H-lines in
    the uc output file after clustering. For search and alignment
    commands both usearch and vsearch reports the target sequence
    number here.

**v2.3.3** released December 5th, 2016
:   A minor speed improvement.

**v2.3.4** released December 9th, 2016
:   Fixed bug in output of sequence profiles and updated documentation.

**v2.4.0** released February 8th, 2017
:   Added support for Linux on Power8 systems (ppc64le) and Windows on
    x86_64. Improved detection of pipes when reading FASTA and FASTQ
    files. Corrected option for specifying output from fastq_eestats
    command in help text.

**v2.4.1** released March 1st, 2017
:   Fixed an overflow bug in fastq_stats and fastq_eestats affecting
    analysis of very large FASTQ files. Fixed maximum memory usage
    reporting on Windows.

**v2.4.2** released March 10th, 2017
:   Default value for fastq_minovlen increased to 16 in accordance
    with help text and for compatibility with usearch. Minor changes
    for improved accuracy of paired-end read merging.

**v2.4.3** released April 6th, 2017
:   Fixed bug with progress bar for shuffling. Fixed missing N-lines in
    UC files with usearch_global, search_exact and allpairs_global when
    the output_no_hits option was not specified.

**v2.4.4** released August 28th, 2017
:   Fixed a few minor bugs, improved error messages and updated
    documentation.

**v2.5.0** released October 5th, 2017
:   Support for UDB database files. New commands: fastq_stripright,
    fastq_eestats2, makeudb_usearch, udb2fasta, udbinfo, and udbstats.
    New general option: no_progress. New options minsize and maxsize to
    fastx_filter. Minor bug fixes, error message improvements and
    documentation updates.

**v2.5.1** released October 25th, 2017
:   Fixed bug with bad default value of 1 instead of 32 for
    minseqlength when using the makeudb_usearch command.

**v2.5.2** released October 30th, 2017
:   Fixed bug where `-` as an argument to the fastq_eestats2 option
    was treated literally instead of equivalent to stdin.

**v2.6.0** released November 10th, 2017
:   Rewritten paired-end reads merger with improved accuracy. Decreased
    default value for fastq_minovlen option from 16 to 10. The default
    value for the fastq_maxdiffs option is increased from 5 to 10.
    There are now other more important restrictions that will avoid
    merging reads that cannot be reliably aligned.

**v2.6.1** released December 8th, 2017
:   Improved parallelisation of paired end reads merging.

**v2.6.2** released December 18th, 2017
:   Fixed option xsize that was partially inactive for commands
    uchime_denovo, uchime_ref, and fastx_filter.

**v2.7.0** released February 13th, 2018
:   Added commands cluster_unoise, uchime2_denovo and uchime3_denovo
    contributed by Davide Albanese based on Robert Edgar's papers.
    Refactored fasta and fastq print functions as well as code for
    extraction of abundance and other attributes from the headers.

**v2.7.1** released February 16th, 2018
:   Fix several bugs on Windows related to large files, use of `-` as a
    file name to mean stdin or stdout, alignment errors, missed kmers
    and corrupted UDB files. Added documentation of UDB-related
    commands.

**v2.7.2** released April 20th, 2018
:   Added the sintax command for taxonomic classification. Fixed a bug
    with incorrect FASTA headers of consensus sequences after
    clustering.

**v2.8.0** released April 24th, 2018
:   Added the fastq_maxdiffpct option to the fastq_mergepairs command.

**v2.8.1** released June 22nd, 2018
:   Fixes for compilation warnings with GCC 8.

**v2.8.2** released August 21st, 2018
:   Fix for wrong placement of semicolons in header lines in some cases
    when using the sizeout or xsize options. Reduced memory requirements
    for full-length dereplication in cases with many duplicate
    sequences. Improved wording of fastq_mergepairs report. Updated
    manual regarding use of sizein and sizeout with dereplication.
    Changed a compiler option.

**v2.8.3** released August 31st, 2018
:   Fix for segmentation fault for `--derep_fulllength` with `--uc`.

**v2.8.4** released September 3rd, 2018
:   Further reduce memory requirements for dereplication when not using
    the uc option. Fix output during subsampling when quiet or log
    options are in effect.

**v2.8.5** released September 26th, 2018
:   Fixed a bug in fastq_eestats2 that caused the values for large
    lengths to be much too high when the input sequences had varying
    lengths.

**v2.8.6** released October 9th, 2018
:   Fixed a bug introduced in version 2.8.2 that caused
    derep_fulllength to include the full FASTA header in its output
    instead of stopping at the first space (unless the notrunclabels
    option is in effect).

**v2.9.0** released October 10th, 2018
:   Added the fastq_join command.

**v2.9.1** released October 29th, 2018
:   Changed compiler options that select the target cpu and tuning to
    allow the software to run on any 64-bit x86 system, while tuning
    for more modern variants. Avoid illegal instruction error on some
    architectures. Update documentation of rereplicate command.

**v2.10.0** released December 6th, 2018
:   Added the sff_convert command to convert SFF files to FASTQ. Added
    some additional option argument checks. Fixed segmentation fault bug
    after some fatal errors when a log file was specified.

**v2.10.1** released December 7th, 2018
:   Improved sff_convert command. It will now read several variants of
    the SFF format. It is also able to read from a pipe. Warnings are
    given if there are minor problems. Error messages have been
    improved. Minor speed and memory usage improvements.

**v2.10.2** released December 10th, 2018
:   Fixed bug in sintax with reversed order of domain and kingdom.

**v2.10.3** released December 19th, 2018
:   Ported to Linux on ARMv8 (aarch64). Fixed compilation warning with
    gcc version 8.1.0 and 8.2.0.

**v2.10.4** released January 4th, 2019
:   Fixed serious bug in x86_64 SIMD alignment code introduced in
    version 2.10.3. Added link to BioConda in README. Fixed bug in
    fastq_stats with sequence length 1. Fixed use of equals symbol in
    UC files for identical sequences with cluster_fast.

**v2.11.0** released February 13th, 2019
:   Added ability to trim and filter paired-end reads using the reverse
    option with the fastx_filter and fastq_filter commands. Added
    `--xee` option to remove ee attributes from FASTA headers. Minor
    invisible improvement to the progress indicator.

**v2.11.1** released February 28th, 2019
:   Minor change to the handling of the weak_id and id options when
    using cluster_unoise.

**v2.12.0** released March 19th, 2019
:   Take sequence abundance into account when computing consensus
    sequences or profiles after clustering. Warn when rereplicating
    sequences without abundance info. Guess offset 33 in more cases
    with fastq_chars. Stricter checking of option arguments and option
    combinations.

**v2.13.0** released April 11th, 2019
:   Added the `--fastx_getseq`, `--fastx_getseqs` and
    `--fastx_getsubseq` commands to extract sequences from a FASTA or
    FASTQ file based on their labels. Improved handling of ambiguous
    nucleotide symbols. Corrected behaviour of `--uchime_ref` command
    with options `--self` and `--selfid`. Strict detection of illegal
    options for each command.

**v2.13.1** released April 26th, 2019
:   Minor changes to the allowed options for each command. All commands
    now allow the log, quiet and threads options. If more than 1 thread
    is specified for commands that are not multi-threaded, a warning
    will be issued. Minor changes to the manual.

**v2.13.2** released April 30th, 2019
:   Fixed bug related to improper handling of newlines on Windows.
    Allowed option strand plus to uchime_ref for compatibility.

**v2.13.3** released April 30th, 2019
:   Fixed bug in FASTQ parsing introduced in version 2.13.2.

**v2.13.4** released May 10th, 2019
:   Added information about support for gzip- and bzip2-compressed
    input files to the output of the version command. Adapted source
    code for compilation on FreeBSD and NetBSD systems.

**v2.13.5** released July 2nd, 2019
:   Added cut command to fragment sequences at restriction sites.
    Silenced output from the fastq_stats command if quiet option was
    given. Updated manual.

**v2.13.6** released July 2nd, 2019
:   Added info about cut command to output of help command.

**v2.13.7** released September 2nd, 2019
:   Fixed bug in consensus sequence introduced in version 2.13.0.

**v2.14.0** released September 11th, 2019
:   Added relabel_self option. Made fasta_width, sizein, sizeout and
    relabelling options valid for certain commands.

**v2.14.1** released September 18th, 2019
:   Fixed bug with sequences written to file specified with fastaout_rev
    for commands fastx_filter and fastq_filter.

**v2.14.2** released January 28th, 2020
:   Fixed some issues with the cut, fastx_revcomp, fastq_convert,
    fastq_mergepairs, and makeudb_usearch commands. Updated manual.

**v2.15.0** released June 19th, 2020
:   Update manual and documentation. Turn on notrunclabels option for
    sintax command by default. Change maxhits 0 to mean unlimited hits,
    like the default. Allow non-ascii characters in headers, with a
    warning. Sort centroids and uc too when clusterout_sort specified.
    Add cluster id to centroids output when clusterout_id specified.
    Improve error messages when parsing FASTQ files. Add missing
    fastq_qminout option and fix label_suffix option for
    fastq_mergepairs. Add derep_id command that dereplicates based on
    both label and sequence. Remove compilation warnings.

**v2.15.1** released October 28th, 2020
:   Fix for dereplication when including reverse complement sequences
    and headers. Make some extra checks when loading compression
    libraries and add more diagnostic output about them to the output
    of the version command. Report an error when fastx_filter is used
    with FASTA input and options that require FASTQ input. Update
    manual.

**v2.15.2** released January 26th, 2021
:   No real functional changes, but some code and compilation changes.
    Compiles successfully on macOS running on Apple Silicon (ARMv8).
    Binaries available. Code updated for C++11. Minor adaptations for
    Windows compatibility, including the use of the C++ standard library
    for regular expressions. Minor changes for compatibility with
    Power8. Switch to C++ header files.

**v2.16.0** released March 22nd, 2021
:   Adds the orient command. Handles empty input files properly.
    Documentation has been updated.

**v2.17.0** released March 29th, 2021
:   The fastq_mergepairs command has been changed. It now allows
    merging of sequences with overlaps as short as 5 bp if the
    `--fastq_minovlen` option has been adjusted down from the default
    10. In addition, much fewer pairs of reads should now be rejected
    with the reason 'multiple potential alignments' as the algorithm
    for detecting those have been changed.

**v2.17.1** released June 14th, 2021
:   Modernized code. Minor changes to help info.

**v2.18.0** released August 27th, 2021
:   Added the fasta2fastq command. Fixed search bug on ppc64le. Fixed
    bug with removal of size and ee info in uc files. Fixed compilation
    errors in some cases. Made some general code improvements. Updated
    manual.

**v2.19.0** released December 21st, 2021
:   Added the lcaout and lca_cutoff options to enable the output of
    last common ancestor (LCA) information about hits when searching.
    The randseed option was added as a valid option to the sintax
    command. Code improvements.

**v2.20.0** released January 10th, 2022
:   Added the fastx_uniques command and the fastq_qout_max option for
    dereplication of FASTQ files. Some code cleaning.

**v2.20.1** released January 11th, 2022
:   Fixes a bug in fastq_mergepair that caused an occasional hang at
    the end when using multiple threads.

**v2.21.0** released January 12th, 2022
:   Adds the sample, qsegout and tsegout options. Enables the use of
    UDB databases with uchime_ref.

**v2.21.1** released January 18th, 2022
:   Fix a problem with dereplication of empty input files. Update
    Altivec code on ppc64le for improved compiler compatibility
    (vector->__vector).

**v2.21.2** released September 12th, 2022
:   Fix problems with the lcaout option when using maxaccepts above 1
    and either lca_cutoff below 1 or with top_hits_only enabled. Update
    documentation. Update code to avoid compiler warnings.

**v2.22.0** released September 19th, 2022
:   Add the derep_smallmem command for dereplication using little
    memory.

**v2.22.1** released September 19th, 2022
:   Fix compiler warning.

**v2.23.0** released July 7th, 2023
:   Update documentation. Add citation file. Modernize and improve
    code. Fix several minor bugs. Fix compilation with GCC 13. Print
    stats after fastq_mergepairs to log file instead of stderr. Handle
    sizein option correctly with dbmatched option for usearch_global.
    Allow maxseqlength option for makeudb_usearch. Fix memory allocation
    problem with chimera detection. Add lengthout and xlength options.
    Increase precision for eeout option. Add warning about sintax
    algorithm, random seed and multiple threads. Refactor chimera
    detection code. Add undocumented experimental long_chimeras_denovo
    command. Fix segfault with clustering. Add more references.

**v2.24.0** released October 26th, 2023
:   Update documentation. Improve code. Allow up to 20 parents for the
    undocumented and experimental chimeras_denovo command. Fix
    compilation warnings for sha1.c. Compile for release (not debug) by
    default.

**v2.25.0** released November 10th, 2023
:   Allow a given percentage of mismatches between chimeras and parents
    for the experimental chimeras_denovo command.

**v2.26.0** released November 24th, 2023
:   Enable the maxseqlength and minseqlength options for the chimera
    detection commands. When the usearch_global or search_exact commands
    are used, OTU tables will include samples and OTUs with no matches.

**v2.26.1** released November 25th, 2023
:   No real changes, but the previous version was released without
    proper updates to the source code.

**v2.27.0** released January 19th, 2024
:   The usearch_global and search_exact commands now support FASTQ files
    as well as FASTA files as input. This version of vsearch includes
    clarifications and updates to the manual. Some code has been
    refactored. Generic Dockerfiles for major Linux distributions have
    been included. Some warnings from compilers and other tools have
    been eliminated. The release for Windows will also include DLLs for
    the two compression libraries.

**v2.27.1** released April 6th, 2024
:   Fixes the weak_id option and makes searches report weak hits in
    some cases. Updates the names of the compression libraries to
    libz.so.1 and libbz2.so.1 on Linux to make them work on common
    Linux distributions without installing additional packages.
    README.md has been updated with information about compression
    libraries on Windows.

**v2.28.0** released April 26th, 2024
:   The sintax command has been improved in several ways in this
    version of vsearch. Please note that several details of this
    algorithm are not clearly described in the preprint, and the
    implementation in vsearch differs from that in usearch. The former
    vsearch version did not always choose the most common taxonomic
    entity over the 100 bootstraps among the database sequences with
    the highest amount of word similarity to the query. Instead, if
    several sequences had an equal similarity with the query, the
    sequence encountered in the earliest bootstrap was chosen. The
    confidence level was calculated based on this sequence compared to
    the selected sequences from the other 99 bootstraps. This could
    lead to a suboptimal choice with a low confidence. In the new
    version, the most common of the sequences with the highest amount
    of word similarity across the 100 bootstraps will be selected, and
    ties will be broken randomly. Another problem with the old
    implementation was that if several sequences had the same amount of
    word similarity, the shortest one in the reference database would
    be chosen, and if they were equally long, the earliest in the
    database file would be chosen. A new option called sintax_random
    has now been introduced. This option will randomly select one of
    the sequences with the highest number of shared words with the
    query, without considering their length or position. This avoids a
    bias towards shorter reference sequences. This option is strongly
    recommended and will probably soon be the default. Furthermore, a
    ninth taxonomic rank, strain (letter t), is now recognized. The
    speed of the sintax command has also been significantly improved at
    least in some cases. Run vsearch with the randseed option and 1
    thread to ensure reproducibility of the random choices in the
    algorithm.

**v2.28.1** released April 26th, 2024
:   Fix a segmentation fault that could occur with the blast6out and
    output_no_hits options.

**v2.29.0** released September 26th, 2024
:   Fixes seven bugs (see changelog below), adds initial support for
    RISC-V architectures, and improves code quality and code testing
    (1,210 new tests):

    - add: experimental support for RISCV64 and other 64-bit
      little-endian architectures, thanks to Michael R. Crusoe and his
      fellow Debian developers (issue #566),
    - add: official support for clang-19 and gcc 14,
    - add: beta support for clang-20,
    - remove: unused `--output` option for command `--fastq_stats`
      (issue #572),
    - fix: bug in `--sintax` when selecting the best lineage (only low
      confidence values below 0.5 were affected) (issue #573),
    - fix: out-of-bounds error in `--fastq_stats` when processing empty
      reads (issue #571),
    - fix: bug in `--cut`, patterns with multiple cutting sites were not
      detected (commit 4c4f9fa),
    - fix: memory error (segmentation fault) when using `--derep_id`
      and `--strand` (issue #565),
    - fix: `--fastq_join` now obeys `--quiet` and `--log` options
      (commit 87f968b),
    - fix: `--fastq_join` quality padding is now also set to Q40 when
      quality offset is 64 (commit be0bf9b),
    - fix: (partial) `--fastq_join`'s handling of abundance annotations
      (commit f2bbcb4),
    - improve: additional safeguards to validate input values and to
      make sure that they are within acceptable limits. Changes concern
      options `--abskew` (commit a530dd8) and `--fastq_maxdiffs`
      (commit 4b254db),
    - improve: code quality (1.3k+ commits, 6k+ clang-tidy warnings
      eliminated),
    - improve: documentation and help messages (issue #568),
    - improve: complete refactoring and modernization of a subset of
      commands (`--sortbylength`, `--sortbysize`, `--shuffle`,
      `--rereplicate`, `--cut`, `--fastq_join`, `--fasta2fastq`,
      `--fastq_chars`),
    - improve: code coverage of our test-suite for the above-mentioned
      commands (1,210 new tests, 4,753 in total).

**v2.29.1** released October 24th, 2024
:   Fix a segmentation fault that could occur during alignment in
    version 2.29.0, for example with `--uchime_ref`. Some improvements
    to code and documentation.

**v2.29.2** released December 20th, 2024
:   Fix a segmentation fault during clustering when the set of clusters
    is empty. Initial documentation in markdown format available on
    GitHub Pages.

**v2.29.3** released February 3rd, 2025
:   Released to mitigate a bug that occurs when compiling the
    `align_simd.cc` file on x86_64 systems with the GNU C++ compiler
    version 9 or later with the `-O3` optimization option. It results
    in incorrect code that may cause bad alignments in some
    circumstances. We are investigating this issue further, but for now
    we recommend compiling with the `-O2` flag. The README.md file and
    the Dockerfiles have been updated to reflect this. The binaries
    released with this version will include this fix.

**v2.29.4** released February 14th, 2025
:   Adjust the window size used for chimera detection down from 64 to
    32. The window size was by accident increased from 32 to 64 in
    version 2.23.0, leading to somewhat fewer chimeras being predicted.
    In addition, a compiler pragma has been included in align_simd.cc
    to further protect the compiler from generating wrong code.

**v2.30.0** released February 27th, 2025
:   Add options `--n_mismatch`, `--fastq_minqual`, and
    `--fastq_truncee_rate`. The `--n_mismatch` option will count N's as
    mismatches in alignments, which may be useful to get sensible
    alignments for sequences with lots of N's. By default N's are
    counted as matches. Both the scoring and the counting of matches
    are affected. The new `--fastq_minqual` option for the
    `--fastq_filter` and `--fastx_filter` commands will discard
    sequences with any bases with a quality score below the given
    value. The new `--fastq_truncee_rate` option for the same commands
    will truncate sequences at the first position where the number of
    expected errors per base is above the given value.

**v2.30.1** released October 3rd, 2025
:   Incorporates many code improvements, more extensive testing, better
    documentation and some minor bug fixes:

    - fix: use-after-free introduced in commit de6c1d8 (Jun 13, 2024),
    - fix: (harmless) out-of-bounds memory issue in `--derep_prefix`
      (commit 8a0a508b),
    - fix: (harmless) memory leak in `--fastx_getseqs --label_field`
      (commit a9c42713),
    - fix: (harmless) memory leak when using option `--userfields`
      (`--allpairs_global` commit 03b95bcf; `--cluster_*` commit
      2fde5472; `--search_exact` commit 45cd56d6; `--usearch_global`
      commit d83bfee9),
    - fix: (harmless) valgrind error, use of uninitialized values
      (commit 8bab2444), also eliminates a pesky compilation warning,
    - change: passing a negative value to `--fastq_truncee_rate` is now
      an error (commit a120f371),
    - change: passing a negative value to `--fastq_minqual` is now an
      error (commit ff5b0c99),
    - change: passing a non-ASCII symbol to `--join_padgap` or
      `--join_padgapq` is now an error (commit a708f5b3),
    - change: when using `--gapopen "*"` to forbid gap opening, the
      penalty is now set to INT_MAX (rather than 1,000). This might
      change alignment results for users who relied on the old behavior
      (thanks to Denis Filloux, issue #602, commit 96e9cf9e),
    - change: when using command `--chimeras_denovo`, `--tabbedout` or
      `--alnout` can be the only output files specified (commit
      d51f0a4),
    - change: when using command `--chimeras_denovo`, option
      `--lengthout` is now accepted (commit 3f55cc6b),
    - change: when using command `--chimeras_denovo`, option `--xlength`
      is now accepted (commit 0fd346cd),
    - add: compilation option GLIBCXX_DEBUG when compiling for
      debugging (commit 5cf4a6c1),
    - add: official support for clang 20,
    - add: initial support for clang 21, initial support for GCC 15,
    - add: experimental support for clang 22,
    - improve: more accurate line number when reporting illegal
      characters in fastq headers (commit 539084e9),
    - improve: more accurate line number when reporting non-ASCII
      characters in fastq headers (commit 98a851ed),
    - improve: remove checks for unneeded libraries during compilation
      (commit 249bb5d5),
    - improve: code quality (8,536 clang-tidy warnings eliminated),
    - improve: documentation and help messages (issue #604),
    - improve: complete refactoring and modernization of the command
      `--fastq_stats`,
    - improve: command `--fastq_stats` is now up to twice faster
      (tested on x86-64),
    - improve: extensive test-suites for `--fastq_stats` and
      `--sff_convert`,
    - improve: code coverage of our test-suite.

**v2.30.2** released December 12th, 2025
:   Fixes two minor issues. Allow a UDB file to be written to a pipe
    or stdout (issue #599). Correct computation of median cluster size
    after dereplication (issue #611).

**v2.30.3** released January 12th, 2026
:   Fixes memory allocation bugs in the chimera detection code that
    caused a segmentation fault in rare cases (issue #615).

**v2.30.4** released January 19th, 2026
:   Fixes issue #617. Due to a bug, incorrect scores were reported in
    the FASTA file headers of some non-chimeric sequences during chimera
    detection. The non-chimeric status was correct.

**v2.30.5** released March 10th, 2026
:   Includes the following changes:

    - fix: out-of-bound look-ups when printing alignment rows (commit
      2f3d1ef). Thanks to user @gbbio for reporting this issue (#618).
      This is a regression, introduced in v2.30.1,
    - fix: options `--xlength` and `--xee` were not removing length and
      ee attributes as expected (commit 69fccdbc). This is a
      regression, introduced in v2.30.1,
    - fix: range checking for the option `--max_unmasked_pct` (commit
      e1a5611505),
    - add: individual manpages for the commands `--shuffle` and `--cut`,
    - improve: eliminated 50 clang-tidy warnings,
    - improve: code coverage of our test-suite.

**v2.30.6** released March 27th, 2026
:   Includes the following changes:

    - fix: out-of-bound look-ups when printing empty last common
      ancestor (LCA) results with the `--lcaout` option. Thanks to user
      AntoninLCH for reporting this issue,
    - improve: code refactoring,
    - improve: code coverage of our test-suite.
