# VSEARCH

## Introduction

The aim of this project is to create an alternative to the [USEARCH](http://www.drive5.com/usearch/) tool developed by Robert C. Edgar (2010).

The new tool should:

* have open source code with an appropriate open source license
* be free of charge, gratis
* have a 64-bit design that handles very large databases and much more than 4GB of memory
* be as accurate or more accurate than usearch
* be as fast or faster than usearch

We have implemented a tool called VSEARCH which supports the following commands: `usearch_global`, `cluster_smallmem`, `cluster_fast`, `derep_fulllength`, `sortbysize`, `sortbylength` and `maskfasta`, as well as almost all their options.

VSEARCH stands for vectorized search, as the tool takes advantage of parallelism in the form of SIMD vectorization to perform accurate alignments at high speed.

Searches have been parallelized using threads and SIMD. VSEARCH includes a SIMD vectorized full global aligner, while USEARCH uses a heuristic seed and extend aligner.

VSEARCH does not support amino acid sequences, local alignments, or chimera detection. These features may be added in the future.

Exactly the same option names as USEARCH version 7 has been used in order to make VSEARCH an almost drop-in replacement.


## Example

In the example below, VSEARCH will identify sequences in the file database.fsa that are at least 90% identical on the plus strand to the query sequences in the file queries.fsa and write the results to the file alnout.txt.

`./vsearch-0.1.0-linux-x86_64 --usearch_global queries.fsa --db database.fsa --strand plus --id 0.9 --alnout alnout.txt`


## Details

**Algorithm:** VSEARCH indexes the unique kmers in the database in a way similar USEARCH, but is currently limited to continuous words (non-spaced seeds). It samples every unique kmer from each query sequence and identifies the number of matching kmers in each database sequence. It then examines the database sequences in order of decreasing number of kmer matches. A full global alignment is computed and those database sequences that satisfy all accept options are accepted while the others are rejected. The `--maxrejects` and `--maxaccepts` options are supported in this process, indicating the maximum number of non-matching and matching databases considered, respectively. Please see the USEARCH paper and supplementary for details.

**Kmer selection:** How many and which kmers USEARCH chooses from the query sequence is not well documented. It is also not known exactly which database sequences are examined, and in which order. We have therefore experimented with various strategies in order to obtain good performance. Our procedure seems to give results equal to or better than USEARCH.

We have chosen to select all unique kmers from the query. At least 6 of these kmers must be present in the database sequence before it will be considered. Also, at least 1 out of 16 query kmers need to be present in the database sequence. Furthermore, if several database sequences have the same number of kmer matches, they will be examined in order of decreasing sequence length.

It appears that there are differences in usearch between the searches performed by the `usearch_global` command and the clustering commands. Notably, it appears like `usearch_global` simply ignores the options `wordlength`, `slots` and `pattern`, while the clustering commands takes them into account.

**Output formats:** Almost all output options of usearch are supported. The output can be written in different formats specified with the `--alnout`, `--blast6out`, `--fastapairs`, `--uc` and the flexible format specified with the `--userout` and `--userfields`. Also the `--matched`, `--notmatched`, `--dbmatched` and `--dbnotmatched` FASTA output files are supported. The only exceptions are the `--consout`, `--construncate` and `--msaout` clustering output options which are not supported.

**Alignment:** VSEARCH uses a 8-way 16-bit SIMD vectorized full dynamic programming algorithm (Needleman-Wunsch) for the global alignments instead of the less sensitive default procedure employed by USEARCH involving seeding, extension and banded dynamic programming. If the `--fulldp` option is specified to USEARCH it will also use a full dynamic programming approach, but USEARCH is then considerably slower.

**Accuracy:** The accuracy of VSEARCH searches has been assessed and compared to USEARCH version 7.0.1090. The Rfam 11.0 database was used for the assessment, as described on the [USEARCH website](http://drive5.com/usearch/benchmark_rfam.html). A similar procedure was described in the USEARCH paper using the Rfam 9.1 database.

The database was initially shuffled. Then the first sequence from each of the 2085 Rfam families with at least two members was selected as queries while the rest was used as the database. The ability of VSEARCH and USEARCH to identify another member of the same family as the top hit was measured, and then recall and precision was calculated.

When USEARCH was run without the `--fulldp` option, VSEARCH had much better recall than USEARCH, but the precision was lower. The [F<sub>1</sub>-score](http://en.wikipedia.org/wiki/F1_score) was considerably higher for VSEARCH. When USEARCH was run with `--fulldp`, VSEARCH had slightly better recall, precision and F-score than USEARCH.

The recall of VSEARCH was usually about 92.3-93.5% and the precision was usually 93.0-94.1%. When run without the `--fulldp` option the recall of USEARCH was usually about 83.0-85.3% while precision was 98.5-99.0%. When run with the `--fulldp` option the recall of USEARCH was usually about 92.0-92.8% and the precision was about 92.2-93.0%.

Please see the files in the `eval` folder for the scripts used for this assessment.

The increased sensitivity of VSEARCH often leads to larger and fewer clusters than USEARCH.
Further improvements in speed may be obtained by intra-sequence SIMD parallelization of the alignments, as well as parallelization of the clustering algorithms.

**Speed:** The speed of VSEARCH searches appears to slightly faster than USEARCH when USEARCH is run without the `--fulldp` option. When USEARCH is run with the `--fulldp` option, VSEARCH may be considerable faster, but it depends on the options and sequences used.

For the accuracy assessment searches in Rfam 11.0 with 100 replicates of the query sequences, VSEARCH needed 51 seconds, whereas USEARCH needed 60 seconds without the `--fulldp` option and 70 seconds with `--fulldp`. This includes time for loading, masking and indexing the database (about 2 secs for VSEARCH, 5 secs for USEARCH). The measurements were made on a Apple MacBook Pro Retina 2013 with four 2.3GHz Intel Core i7 cores (8 virtual cores) using the default number of threads (8).

The speed of clustering with VSEARCH relative to USEARCH depends on how many threads are used. Running with a single thread VSEARCH currently seems to be 2-4 times slower than with USEARCH, depending on parameters. Clustering has been parallelized with threads in VSEARCH, but clustering does not seem to be parallelized in USEARCH (despite what the name and documentation for `cluster_fast` seems to indicate). Clustering with VSEARCH using 4-8 threads is often faster than USEARCH. The speed of VSEARCH might be further improved with an intra-sequence SIMD-vectorized aligner.

The dereplication and sorting commands seems to be considerably faster in VSEARCH than in USEARCH.

**Clustering:** The clustering commands `cluster_smallmem` and `cluster_fast` have been implemented with all options except `--consout`, `--msaout` and `--construncate`.

**Masking:** VSEARCH by default uses an optimzed reimplementation of the well-known DUST algorithm by Tatusov and Lipman to mask simple repeats and low-complexity regions in the sequences. USEARCH by default uses an undocumented rapid masking method called "fastnucleo" that seems to mask fewer and smaller regions. USEARCH may also be run with the DUST masking method, but the masking then takes something like 30 times longer.

**Extensions:** A shuffle command has been added. By specifying a FASTA file using the `--shuffle` option, and an output file with the `--output` option, VSEARCH will shuffle the sequences in a pseudo-random order. A positive integer may be specified as the seed with the `--seed` option to generate the same shuffling several times. By default, or when `--seed 0` is specified, the pseudo-random number generator will be initialized with pseudo-random data from the machine to give different numbers each time it is run.

Another extension implemented is that dereplication will honor the `--sizein` option and add together the abundances of the sequences that are merged.

The width of FASTA formatted output files may be specified with the `--fasta_width` option. Both the `--fasta_width` option and the `--rowlen` option will not wrap sequence and alignments when a value of zero is specified.

VSEARCH implements the old USEARCH option `--iddef` to specify the definition of identity used to rank the hits. Values accepted are 0 (CD-HIT definition using shortest sequence as numerator), 1 (edit distance), 2 (edit distance excluding terminal gaps, default), 3 (Marine Biological Lab definition where entire gaps are considered a single difference) or 4 (BLAST, same as 2). See the [USEARCH User Guide 4.1](http://drive5.com/usearch/UsearchUserGuide4.1.pdf) page 42-44 for details. Also id0 to id4 are accepted as arguments to the `--userfields` option.

**Command line options:** The options currently supported by VSEARCH are indicated below. Please run VSEARCH with the `--help` option to see more information about the options.


## Command line options supported

General options:

* `--help`
* `--version`
* `--iddef <int>` (Default 2)
* `--fasta_width <int>` (Default 80)
* `--maxseqlength <int>` (Default 50000)
* `--minseqlength <int>` (Default 1 for sort/shuffle or 32 for search/dereplicate)
* `--notrunclabels`
* `--strand <plus|both>`
* `--threads <int>` (Default 0=available cores)
* `--uc <filename>`
* `--uc_allhits`

Clustering and searching options:

* `--cluster_fast <filename>`
* `--cluster_smallmem <filename>`
* `--usearch_global <filename>`
* `--alnout <filename>`
* `--blast6out <filename>`
* `--centroids <filename>`
* `--clusters <prefix>`
* `--consout <filename>`
* `--construncate` (Not implemented)
* `--db <filename>` (Required)
* `--dbmask dust|none|soft` (Default dust)
* `--dbmatched <filename>`
* `--dbnotmatched <filename>`
* `--fastapairs <filename>`
* `--fulldp` (Ignored; VSEARCH always computes full dynamic programming alignments)
* `--gapext <string>` (Default 2I/1E)
* `--gapopen <string>` (Default 20I/2E)
* `--hardmask`
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
* `--maxrejects <int>` (Default 32)
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
* `--msaout <filename>` (Not implemented)
* `--notmatched <filename>`
* `--output_no_hits`
* `--qmask dust|none|soft` (Default dust)
* `--query_cov <real>`
* `--rightjust`
* `--rowlen <int>` (Default 64)
* `--self`
* `--selfid`
* `--target_cov <real>`
* `--top_hits_only`
* `--userfields <string>`
* `--userout <filename>`
* `--usersort`
* `--weak_id <real>`
* `--wordlength <int>` (Default 8)

Dereplication, masking, shuffling and sorting options:

* `--derep_fulllength <filename>`
* `--maskfasta <filename>`
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

* VSEARCH cannot perform chimera detection.


## VSEARCH license and third party licenses

The code is currently licensed under the GNU Affero General Public License version 3.

VSEARCH includes code from Google's [CityHash project](http://code.google.com/p/cityhash/) by Geoff Pike and Jyrki Alakuijala, providing some excellent hash functions available under a MIT license.

VSEARCH includes code derived from Tatusov and Lipman's DUST program that is in the public domain.

VSEARCH binaries may include code from the [zlib](http://www.zlib.net) library copyright Jean-loup Gailly and Mark Adler.

VSEARCH binaries may include code from the [bzip2](http://www.bzip.org) library copyright Julian R. Seward.


## Code

The code is written in C++ but most of it is actually C with some C++ syntax conventions.

    File     | Description
-------------|------
**align.cc** | New Needleman-Wunsch global alignment, serial. Only for testing.
**alignsimd.cc** | SIMD parallel global alignment of 1 query with 8 database sequences
**arch.cc** | Architecture specific code (Mac/Linux)
**bitmap.cc** | Implementation of bitmaps
**cluster.cc** | Clustering (cluster\_fast and cluster\_smallmem)
**db.cc** | Handles the database file read, access etc
**dbindex.cc** | Indexes the database by identifying unique kmers in the sequences and make a database hash table
**derep.cc** | Dereplication
**maps.cc** | Various character mapping arrays
**mask.cc** | Masking (DUST)
**minheap.cc** | A minheap implementation for the list of top kmer matches
**query.cc** | Reads the fasta file containing the query sequences.
**results.cc** | Output results in various formats (alnout, userout, blast6, uc)
**search.cc** | Implements search using global alignment
**searchcore.cc** | Core search functions for searching and clustering
**showalign.cc** | Output an alignment in a human-readable way given a CIGAR-string and the sequences
**shuffle.cc** | Shuffle sequences
**sortbylength.cc** | Code for sorting by length
**sortbysize.cc** | Code for sorting by size (abundance)
**unique.cc** | Find unique kmers in a sequence
**userfields.cc** | Code for parsing the userfields option argument
**util.cc** | Various common utility functions
**vsearch.cc** | Main program file, general initialization, reads arguments and parses options, writes info.

VSEARCH may be compiled with zlib or bzip2 integration that allows it to read compressed FASTA files. The [zlib](http://www.zlib.net/) and the [bzip2](http://www.bzip.org/) libraries are needed for this.


## Bugs

VSEARCH has not been tested comprehensively yet. All bug reports are highly appreciated.


## Future work

Some issues to work on:

* testing and debugging
* intra-sequence SIMD parallelization using the striped approach (Farrar 2007) or the plain vertical approach (Rognes & Seeberg 2000)
* multiple sequence alignment and related options for consensus sequences from clustering
* chimera filtering


## The VSEARCH team

The following people have contributed to VSEARCH:

* Torbj&oslash;rn Rognes
* Tom&aacute;&scaron; Flouri
* Fr&eacute;d&eacute;ric Mah&eacute;
* Christopher Quince
* Umer Zeeshan Ijaz


## References

Edgar, Robert C. (2010) Search and clustering orders of magnitude faster than BLAST. Bioinformatics, 26 (19): 2460-2461. doi:[10.1093/bioinformatics/btq461](http://dx.doi.org/10.1093/bioinformatics/btq461)
