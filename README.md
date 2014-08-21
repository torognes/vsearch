# VSEARCH


## Introduction

The aim of the project is to create an alternative to the [USEARCH](http://www.drive5.com/usearch/) tool developed by Robert C. Edgar (2010). The new tool should have:

* open source code with an appropriate open source license
* 64-bit design that handles very large databases and more than 4GB of memory

A tool called VSEARCH has been implemented. Exactly the same option names as USEARCH has been used in order to make it possible to make VSEARCH almost a drop-in replacement.
The basic `usearch_global` algorithm for global alignments using nucleotide sequences is implemented, as well as the `derep_fulllength` and the `sortbysize` and `sortbylength` commands.
At this stage it does not support amino acid sequences, local alignments, clustering, chimera detection etc. It currently does not use any parallellisation (neither pthreads nor SIMD).

In the example below, VSEARCH will identify sequences in database.fsa at least 90% identical on the plus strand to the query sequences in queries.fsa and write the results to alnout.txt.

`./vsearch-0.0.8-linux-x86_64 --usearch_global queries.fsa --db database.fsa --strand plus --id 0.9 --alnout alnout.txt`


## Implementation details

**Algorithm:** VSEARCH indexes the unique kmers in the database in the same way as USEARCH, but is currently limited to continuous words (non-spaced seeds). It evenly samples a certain number of unique kmers from each query sequence and identifies the database sequences with the largest number of kmer matches. It then examines the database sequences in order of decreasing number of kmer matches. A full global alignment is computed and those database sequences that satisfy all accept options are accepted while the others are rejected. The `--maxrejects` and `--maxaccepts` options are supported in this process, indicating the maximum number of non-matching and matching databases considered, respectively. Please see the USEARCH paper and supplementary for details.

**Kmer selection:** How many and which unique kmers USEARCH chooses from the query sequence is not well documented, but our procedure seems to give results approximately equal to USEARCH. According to the documentation, USEARCH is supposed to sample *x* of the unique kmers in the query uniformly along the query sequence, where *x* is chosen so as to expect 8 matches with the targets satisfying the accept criteria (specified with the `--id` option). The exact procedure is unknown. The choice of *x* should be looked further into. Currently, VSEARCH samples up to 255 unique kmers from the query.

**Output formats:** All the output options of usearch are supported. The output can be written in five different formats specified with the `--alnout`, `--blast6out`, `--fastapairs`, `--uc` and the flexible format specified with the `--userout` and `--userfields`. Also the `--matched`, `--notmatched`, `--dbmatched` and `--dbnotmatched` FASTA output files are supported.

**Alignment:** VSEARCH currently uses a relatively slow non-vectorized full dynamic programming Needleman-Wunsch algorithm for the global alignments (similar to USEARCH with `--fulldp`) instead of the quicker (but possible less senstive) procedure involving seeding, extension and banded dynamic programming employed by USEARCH. This could be replaced by a vectorized variant using SIMD in VSEARCH to get up to comparable speed. A banded SIMD variant could also be valuable.

**Performance:** Based on very limited testing, the VSEARCH speed appears to be about half that of USEARCH when USEARCH is run with the `--fulldp` option, but it depends a lot on the exact parameters and sequences. With `--id 1.0` the are almost equal in speed, but with lower thresholds, e.g. `--id 0.85` USEARCH is considerably faster, but also much less senstive. The accuracy generally seems comparable to USEARCH. But all of this needs to be looked more into. The dereplication and sorting commands seems to be considerably faster in VSEARCH than in USEARCH.

**Command line options:** A list of all the options supported is available in the doc folder. Please see the file usearch\_options.md. The options currently supported by VSEARCH is indicated below and also in the options_supported.md


## Command line options supported

General options:

* `--help`
* `--version`
* `--threads <int>` (Default 1) (Currently only one thread supported)
* `--maxseqlength <int>` (Default 50000)
* `--notrunclabels`

Search options:

* `--usearch_global <filename>`
* `--db <filename>`
* `--strand <plus|both>`
* `--maxhits`
* `--top_hits_only`
* `--output_no_hits`
* `--alnout <filename>`
* `--rowlen <int>` (Default 64)
* `--userout <filename>`
* `--userfields <string>`
* `--uc <filename>`
* `--uc_allhits`
* `--fastapairs <filename>`
* `--blast6out <filename>`
* `--matched <filename>`
* `--notmatched <filename>`
* `--dbmatched <filename>`
* `--sizeout`
* `--dbnotmatched <filename>`
* `--id <real>`
* `--weak_id <real>`
* `--query_cov <real>`
* `--target_cov <real>`
* `--idprefix <int>`
* `--idsuffix <int>`
* `--minqt <real>`
* `--maxqt <real>`
* `--minsl <real>`
* `--maxsl <real>`
* `--leftjust`
* `--rightjust`
* `--self`
* `--selfid`
* `--maxid <real>`
* `--minsizeratio <real>`
* `--maxsizeratio <real>`
* `--maxdiffs <int>`
* `--maxsubs <int>`
* `--maxgaps <int>`
* `--mincols <int>`
* `--maxqsize <int>`
* `--mintsize <int>`
* `--mid <real>`
* `--maxaccepts <int>` (Default 1)
* `--maxrejects <int>` (Default 16)
* `--wordlength <int>` (Default 8)
* `--match <int>` (Default 2)
* `--mismatch <int>` (Default -4)
* `--gapopen <string>` (Default 20I/2E)
* `--gapext <string>` (Default 2I/1E)
* `--fulldp` (Always full dynamic programming alignments)
* `--minseqlength <int>` (Default 1)

Dereplication options:

* `--derep_fulllength <filename>`
* `--minseqlength <int>` (Default 32)
* `--output <filename>`
* `--uc <filename>`
* `--uc_allhits`
* `--sizein`
* `--sizeout`
* `--minuniquesize <int>`
* `--strand <plus|both>`
* `--topn <int>`

Abundance sort options:

* `--sortbysize <filename>`
* `--output <filename>`
* `--minseqlength <int>` (Default 32)
* `--minsize <int>` (Default 0)
* `--maxsize <int>` (Default inf.)
* `--relabel`
* `--sizeout`

Length sort options:

* `--sortbylength <filename>`
* `--output <filename>`
* `--minseqlength <int>` (Default 1)
* `--relabel`
* `--sizeout`


## Main limitations

* **Commands:** No clustering or chimera checking, yet.
* **Speed:** Only non-vectorized full alignment, no parallelization. Threads are currently not supported.
* **Masking:** Currently, VSEARCH does not mask the sequences while USEARCH performs masking by default.
* **Indexing options:** Only continuous seeds are supported.

## License

The code is currently licensed under the GNU Affero General Public License version 3. The ordinary GNU GPL is one alternative. We should discuss this.


## Code

The code is written in C++ but most of it is actually C with some C++ syntax conventions. The files are:

* **vsearch.h** - C header file for entire project
* **vsearch.cc** - Main program file, general initialization, reads arguments and parses options, writes info.
* **query.cc** - reads the fasta file containing the query sequences.
* **bzquery.cc** - reads compressed fasta query files. Not currently used.
* **db.cc** - Handles the database file read, access etc
* **dbindex.cc** - Indexes the database by identifying unique kmers in the sequences and make a database hash table
* **nw.cc** - New Needleman-Wunsch global alignment, serial
* **nws.cc** - Old Needleman-Wunsch-Sellers global alignment, serial. Not used.
* **showalign.cc** - Output an alignment in a human-readable way given a CIGAR-string and the sequences
* **search.cc** - search database
* **util.cc** - Various common utility functions
* **results.cc** - Output results in various formats (alnout, userout, blast6, uc)
* **derep.cc** - Code for dereplication
* **userfields.cc** - Code for parsing the userfields option argument
* **arch.cc** - Architecture specific code (Mac/Linux)
* **maps.cc** - Various character mapping arrays
* **sortbylength.cc** - Code for sorting by length
* **sortbysize.cc** - Code for sorting by size (abundance)

## Bugs

VSEARCH has not been tested comprehensively. All bug reports are appreciated.


## Future work

Some issues to work on:

* testing and debugging
* precision and recall comparison with USEARCH
* performance (speed) comparison with USEARCH
* optimize selection of unique kmers
* optimize order of database sequence comparison
* parallelisation with pthreads
* parallelisation with SIMD-based global alignment
* clustering
* chimera filtering
* masking


## References

Edgar, Robert C. (2010) Search and clustering orders of magnitude faster than BLAST. Bioinformatics, 26 (19): 2460-2461. doi:[10.1093/bioinformatics/btq461](http://dx.doi.org/10.1093/bioinformatics/btq461)
