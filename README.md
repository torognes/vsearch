# VSEARCH

## Introduction

The aim of the project is to create an alternative to the [USEARCH](http://www.drive5.com/usearch/) tool developed by Robert C. Edgar (2010). The new tool should have:

* open source code with an appropriate open source license
* 64-bit design that handles very large databases and more than 4GB of memory

A tool called VSEARCH has been implemented. Exactly the same option names as USEARCH has been used in order to make it possible to make VSEARCH almost a drop-in replacement.
The basic `usearch_global` algorithm for global alignments using nucleotide sequences is implemented, as well as the `derep_fulllength` and the `sortbysize` and `sortbylength` commands.
At this stage it does not support amino acid sequences, local alignments, clustering, chimera detection etc.
Searches have been parallelized using threads and SIMD. VSEARCH includes a SIMD vectorized full global aligner, while USEARCH uses heuristic aligner.

In the example below, VSEARCH will identify sequences in database.fsa at least 90% identical on the plus strand to the query sequences in queries.fsa and write the results to alnout.txt.

`./vsearch-0.0.12-linux-x86_64 --usearch_global queries.fsa --db database.fsa --strand plus --id 0.9 --alnout alnout.txt`


## Implementation details

**Algorithm:** VSEARCH indexes the unique kmers in the database in a way similar USEARCH, but is currently limited to continuous words (non-spaced seeds). It samples every unique kmer from each query sequence and identifies the number of kmer matches in each database sequence. It then examines the database sequences in order of decreasing number of kmer matches. A full global alignment is computed and those database sequences that satisfy all accept options are accepted while the others are rejected. The `--maxrejects` and `--maxaccepts` options are supported in this process, indicating the maximum number of non-matching and matching databases considered, respectively. Please see the USEARCH paper and supplementary for details.

**Kmer selection:** How many and which kmers USEARCH chooses from the query sequence is not well documented. It is also not known exactly which database sequences are examined, and in which order. We have therefore experimented with various strategies in order to obtain good performance. Our procedure seems to give results equal to or better than USEARCH. 
We have chosen to select all unique kmers from the query. At least 6 of these kmers must be present in the database sequence before it will be considered. Also, at least 1 out of 16 query kmers need to be present in the database sequence. Furthermore, if several database sequences have the same number of kmer matches, they will be examined in order of decreasing sequence length.

**Output formats:** All the output options of usearch are supported. The output can be written in different formats specified with the `--alnout`, `--blast6out`, `--fastapairs`, `--uc` and the flexible format specified with the `--userout` and `--userfields`. Also the `--matched`, `--notmatched`, `--dbmatched` and `--dbnotmatched` FASTA output files are supported.

**Alignment:** VSEARCH currently uses a 8-way SIMD vectorized full dynamic programming algorithm (Needleman-Wunsch) for the global alignments instead of the less sensitive default procedure employed by USEARCH involving seeding, extension and banded dynamic programming. If the `--fulldp` option is specified to USEARCH it will also use a full dynamic programming approach, but USEARCH is then considerably slower.

**Accuracy:** The accuracy of VSEARCH has been assessed and compared to USEARCH. The Rfam 11.0 database was used for the assessment, as described on the [USEARCH website](http://drive5.com/usearch/benchmark_rfam.html). A similar procedure was described in the USEARCH paper using the Rfam 9.1 database.
The database was initially shuffled. Then the first sequence from each of the 2085 Rfam families with at least two members was selected as queries while the rest was used as the database. The ability of VSEARCH and USEARCH to identify another member of the same family as the top hit was measured. Recall and precision was calculated. In most cases VSEARCH had slightly better recall and precision than USEARCH.
The recall of VSEARCH was usually about 92.3-93.5% and the precision was usually 93.0-94.1%. When run without the `--fulldp` option the recall of USEARCH was usually about 83.0-85.3% while precision was 98.5-99.0%. When run with the `--fulldp` option the recall of USEARCH was usually about 92.0-92.8% and the precision was about 92.2-93.0%.
Please see the files in the `eval` folder for the scripts used for this assessment.

**Speed:** The speed of VSEARCH appears about equal to USEARCH when USEARCH is run without the `--fulldp` option. When USEARCH is run with the `--fulldp` option, VSEARCH may be considerable faster, but it depends on the options and sequences used.
For the accuracy assessment searches in Rfam 11.0, VSEARCH took 71 seconds for 100 replicates of the same query sequences, whereas USEARCH without the `--fulldp` option needed 63 seconds and USEARCH with `--fulldp` needed 70 seconds. This includes time for loading and indexing the database (about 3 secs for VSEARCH, 6 secs for USEARCH). The measurements were made on a Apple MacBook Pro Retina 2013 with four 2.3GHz Intel Core i7 cores (8 virtual cores) using the default number of threads (8).
The dereplication and sorting commands seems to be considerably faster in VSEARCH than in USEARCH.

**Command line options:** The options currently supported by VSEARCH is indicated below. Please run VSEARCH with the `--help` option to see more information about the options.

**Extensions:** A shuffle command has been added. By specifying a FASTA file using the `--shuffle` option, and an output file with the `--output` option, VSEARCH will shuffle the sequences in a pseudo-random order. A positive integer may be specified as the seed with the `--seed` option to generate the same shuffling several times. By default, or when `--seed 0` is specified, the pseudo-random number generator will be initialized with pseudo-random data from the machine to give different numbers each time it is run.
Another extension implemented is that dereplication will honor the `--sizein` option and add together the abundances of the sequences that are merged.
The width of FASTA formatted output files may be specified with the `--fasta_width` option.


## Command line options supported

General options:

* `--help`
* `--version`
* `--fasta_width <int>` (Default 80)
* `--maxseqlength <int>` (Default 50000)
* `--minseqlength <int>` (Default 1 for sort/shuffle or 32 for search/dereplicate)
* `--notrunclabels`
* `--strand <plus|both>` (Required for `search_global`)
* `--threads <int>` (Default 0=available cores)
* `--uc <filename>`
* `--uc_allhits`

Search options:

* `--usearch_global <filename>`
* `--alnout <filename>`
* `--blast6out <filename>`
* `--db <filename>` (Required)
* `--dbmatched <filename>`
* `--dbnotmatched <filename>`
* `--fastapairs <filename>`
* `--fulldp` (VSEARCH always computes full dynamic programming alignments)
* `--gapext <string>` (Default 2I/1E)
* `--gapopen <string>` (Default 20I/2E)
* `--id <real>` (Required)
* `--idprefix <int>`
* `--idsuffix <int>`
* `--leftjust`
* `--match <int>` (Default 2)
* `--matched <filename>`
* `--maxaccepts <int>` (Default 1)
* `--maxdiffs <int>`
* `--maxgaps <int>`
* `--maxhits`
* `--maxid <real>`
* `--maxqsize <int>`
* `--maxqt <real>`
* `--maxrejects <int>` (Default 16)
* `--maxsizeratio <real>`
* `--maxsl <real>`
* `--maxsubs <int>`
* `--mid <real>`
* `--mincols <int>`
* `--minqt <real>`
* `--minsizeratio <real>`
* `--minsl <real>`
* `--mintsize <int>`
* `--mismatch <int>` (Default -4)
* `--notmatched <filename>`
* `--output_no_hits`
* `--query_cov <real>`
* `--rightjust`
* `--rowlen <int>` (Default 64)
* `--self`
* `--selfid`
* `--target_cov <real>`
* `--top_hits_only`
* `--userfields <string>`
* `--userout <filename>`
* `--weak_id <real>`
* `--wordlength <int>` (Default 8)

Dereplication, sorting and shuffling options:

* `--derep_fulllength <filename>`
* `--shuffle <filename>`
* `--sortbylength <filename>`
* `--sortbysize <filename>`
* `--maxsize <int>` (Default inf.)
* `--minsize <int>` (Default 0)
* `--minuniquesize <int>`
* `--output <filename>` (Required for `shuffle`, `sortbylength` and `sortbysize`)
* `--relabel`
* `--seed <int>` (Default 0=randomize)
* `--sizein`
* `--sizeout`
* `--strand <plus|both>`
* `--topn <int>`


## Main limitations

* **Commands:** No clustering or chimera checking, yet.
* **Masking:** Currently, VSEARCH does not mask the sequences while USEARCH performs masking by default.
* **Indexing options:** Only continuous seeds are supported.


## License

The code is currently licensed under the GNU Affero General Public License version 3.

VSEARCH includes code from Google's [CityHash project](http://code.google.com/p/cityhash/) by Geoff Pike and Jyrki Alakuijala, providing some excellent hash functions available under a MIT license.

VSEARCH binaries may include code from the [zlib](http://www.zlib.net) library copyright Jean-loup Gailly and Mark Adler.

VSEARCH binaries may include code from the [bzip2](http://www.bzip.org) library copyright Julian R. Seward.


## Code

The code is written in C++ but most of it is actually C with some C++ syntax conventions. The files are:

* **vsearch.h** - C header file for entire project
* **arch.cc** - Architecture specific code (Mac/Linux)
* **db.cc** - Handles the database file read, access etc
* **dbindex.cc** - Indexes the database by identifying unique kmers in the sequences and make a database hash table
* **derep.cc** - Code for dereplication
* **maps.cc** - Various character mapping arrays
* **nw.cc** - New Needleman-Wunsch global alignment, serial. Only for testing.
* **query.cc** - reads the fasta file containing the query sequences.
* **results.cc** - Output results in various formats (alnout, userout, blast6, uc)
* **search.cc** - Search database
* **searchsimd.cc** - SIMD parallel global alignment of 1 query with 8 database sequences
* **showalign.cc** - Output an alignment in a human-readable way given a CIGAR-string and the sequences
* **shuffle.cc** - Shuffle sequences
* **sortbylength.cc** - Code for sorting by length
* **sortbysize.cc** - Code for sorting by size (abundance)
* **unique.cc** - Find unique kmers in a sequence
* **userfields.cc** - Code for parsing the userfields option argument
* **util.cc** - Various common utility functions
* **vsearch.cc** - Main program file, general initialization, reads arguments and parses options, writes info.

VSEARCH may be compiled with zip or bzip2 integration that allows it to read compressed FASTA files. The [zlib](http://www.zlib.net/) and the [bzip2](http://www.bzip.org/) libraries are needed for this.


## Bugs

VSEARCH has not been tested comprehensively yet. All bug reports are highly appreciated.


## Future work

Some issues to work on:

* testing and debugging
* masking
* clustering
* chimera filtering


## The VSEARCH Team

The following people have contributed to VSEARCH:

* Torbjørn Rognes
* Tomás Flouri
* Frédéric Mahé
* Christopher Quince
* Umer Zeeshan Ijaz


## References

Edgar, Robert C. (2010) Search and clustering orders of magnitude faster than BLAST. Bioinformatics, 26 (19): 2460-2461. doi:[10.1093/bioinformatics/btq461](http://dx.doi.org/10.1093/bioinformatics/btq461)
