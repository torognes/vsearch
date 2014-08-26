# VSEARCH

## Introduction

The aim of the project is to create an alternative to the [USEARCH](http://www.drive5.com/usearch/) tool developed by Robert C. Edgar (2010). The new tool should have:

* open source code with an appropriate open source license
* 64-bit design that handles very large databases and more than 4GB of memory

A tool called VSEARCH has been implemented. Exactly the same option names as USEARCH has been used in order to make it possible to make VSEARCH almost a drop-in replacement.
The basic `usearch_global` algorithm for global alignments using nucleotide sequences is implemented, as well as the `derep_fulllength` and the `sortbysize` and `sortbylength` commands.
At this stage it does not support amino acid sequences, local alignments, clustering, chimera detection etc. Searches have been parallelized with threads. SIMD-parallelization of global alignments is planned, but not yet implemented.

In the example below, VSEARCH will identify sequences in database.fsa at least 90% identical on the plus strand to the query sequences in queries.fsa and write the results to alnout.txt.

`./vsearch-0.0.10-linux-x86_64 --usearch_global queries.fsa --db database.fsa --strand plus --id 0.9 --alnout alnout.txt`


## Implementation details

**Algorithm:** VSEARCH indexes the unique kmers in the database in a way probably similar USEARCH, but is currently limited to continuous words (non-spaced seeds). It samples every unique kmer from each query sequence and identifies the database sequences with the largest number of kmer matches. It then examines the database sequences in order of decreasing number of kmer matches. A full global alignment is computed and those database sequences that satisfy all accept options are accepted while the others are rejected. The `--maxrejects` and `--maxaccepts` options are supported in this process, indicating the maximum number of non-matching and matching databases considered, respectively. Please see the USEARCH paper and supplementary for details.

**Kmer selection:** How many and which kmers USEARCH chooses from the query sequence is not well documented. It is also not known exactly which database sequences are examined, and in which order. We have therefore experimented with various possibilities in order to obtain good performance. Our procedure seems to give results equal to or better than USEARCH. 
We have choosen to select all unique kmers from the query. At least 7 of these kmers must be present in the database sequence before it will be considered. Also, at least 1 of every 16 query kmers need to be present in the database sequence. Furthermore, if several database sequences have the same number of kmer matches, they will be examined in order of decreasing length.

**Output formats:** All the output options of usearch are supported. The output can be written in five different formats specified with the `--alnout`, `--blast6out`, `--fastapairs`, `--uc` and the flexible format specified with the `--userout` and `--userfields`. Also the `--matched`, `--notmatched`, `--dbmatched` and `--dbnotmatched` FASTA output files are supported.

**Alignment:** VSEARCH currently uses a relatively slow non-vectorized full dynamic programming Needleman-Wunsch algorithm for the global alignments (similar to USEARCH with `--fulldp`) instead of the quicker (but possible less senstive) procedure involving seeding, extension and banded dynamic programming employed by USEARCH. This could be replaced by a vectorized variant using SIMD in VSEARCH to get up to comparable speed. A banded SIMD variant could also be valuable.

**Speed:** The speed of VSEARCH appears to be roughly equal to USEARCH when USEARCH is run with the `--fulldp` option. When USEARCH is run without the `--fulldp` option, it may be considerable faster, but it is quite variable and depends on several parameters. Running USEARCH without the `--fulldp` option also seems to results in reduced recall.
The dereplication and sorting commands seems to be considerably faster in VSEARCH than in USEARCH.

**Accuracy:** The accuracy of VSEARCH has been assessed and compared to USEARCH. The Rfam 11.0 database was used for the assessment, as described on the [USEARCH website](http://drive5.com/usearch/benchmark_rfam.html). A similar procedure was described in the USEARCH paper with the Rfam 9.1 database.
The database was initially shuffled. Then the first sequence from each of 2085 Rfam families with more than one member was selected. The ability of VSEARCH and USEARCH to select another member of the same familiy as the top hit was measured. Recall and precision was calculated. In most cases VSEARCH had somewhat better recall and precision than USEARCH. Please see the files in the `eval` folder for details.

**Command line options:** The options currently supported by VSEARCH is indicated below. Please run vsearch with the `--help` option to see more information about the options.

**Extensions:** A shuffle command has been added. By specifiying a FASTA file using the `--shuffle` option, and an output file with the `--output` option, VSEARCH will shuffle the sequences pseudo-randomly. A seed by be specified with the `--seed` option to generate the same shuffling several times. By default, or when `--seed 0` is specified, the pseudo-random generator will be initialized with pseudo-random data from the local machine.
Another extension implemented is that dereplication will honor the `--sizein` option and add together the abundances of the sequences that are merged.
The width of FASTA formatted output files may be specified with the `--fasta_width` option.


## Command line options supported

General options:

* `--help`
* `--version`
* `--maxseqlength <int>` (Default 50000)
* `--minseqlength <int>` (Default 1 for sort/shuffle or 32 for search/dereplicate)
* `--notrunclabels`
* `--strand <plus|both>`
* `--threads <int>` (Default 0=available cores)
* `--uc <filename>`
* `--uc_allhits`

Search options:

* `--usearch_global <filename>`
* `--alnout <filename>`
* `--blast6out <filename>`
* `--db <filename>`
* `--dbmatched <filename>`
* `--dbnotmatched <filename>`
* `--fasta_width <int>` (Default 80)
* `--fastapairs <filename>`
* `--fulldp` (Always full dynamic programming alignments)
* `--gapext <string>` (Default 2I/1E)
* `--gapopen <string>` (Default 20I/2E)
* `--id <real>`
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
* `--output <filename>`
* `--relabel`
* `--seed <int>` (Default 0=randomize)
* `--sizein`
* `--sizeout`
* `--strand <plus|both>`
* `--topn <int>`


## Main limitations

* **Commands:** No clustering or chimera checking, yet.
* **Speed:** Only non-vectorized full alignment implemented. A vectorized version (SIMD) will be made. Otherwize VSEARCH is generally faster than USEARCH.
* **Masking:** Currently, VSEARCH does not mask the sequences while USEARCH performs masking by default.
* **Indexing options:** Only continuous seeds are supported.


## License

The code is currently licensed under the GNU Affero General Public License version 3.


## Code

The code is written in C++ but most of it is actually C with some C++ syntax conventions. The files are:

* **vsearch.h** - C header file for entire project
* **arch.cc** - Architecture specific code (Mac/Linux)
* **bzquery.cc** - reads compressed fasta query files. Not currently used.
* **db.cc** - Handles the database file read, access etc
* **dbindex.cc** - Indexes the database by identifying unique kmers in the sequences and make a database hash table
* **derep.cc** - Code for dereplication
* **maps.cc** - Various character mapping arrays
* **nw.cc** - New Needleman-Wunsch global alignment, serial
* **nws.cc** - Old Needleman-Wunsch-Sellers global alignment, serial. Not used.
* **query.cc** - reads the fasta file containing the query sequences.
* **results.cc** - Output results in various formats (alnout, userout, blast6, uc)
* **search.cc** - search database
* **showalign.cc** - Output an alignment in a human-readable way given a CIGAR-string and the sequences
* **shuffle.cc** - Shuffle sequences
* **sortbylength.cc** - Code for sorting by length
* **sortbysize.cc** - Code for sorting by size (abundance)
* **unique.cc* - Find unique kmers in a sequence
* **userfields.cc** - Code for parsing the userfields option argument
* **util.cc** - Various common utility functions
* **vsearch.cc** - Main program file, general initialization, reads arguments and parses options, writes info.


## Bugs

VSEARCH has not been tested comprehensively yet. All bug reports are highly appreciated.


## Future work

Some issues to work on:

* testing and debugging
* parallelisation with SIMD-based global alignment
* masking
* clustering
* chimera filtering


## References

Edgar, Robert C. (2010) Search and clustering orders of magnitude faster than BLAST. Bioinformatics, 26 (19): 2460-2461. doi:[10.1093/bioinformatics/btq461](http://dx.doi.org/10.1093/bioinformatics/btq461)
