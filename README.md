# VSEARCH


## Introduction

The aim of the project is to create an alternative to the USEARCH tool. The new tool should have:

* open source code with an appropriate open source license
* 64-bit design that handles very large databases and more than 4GB of memory

A tool called VSEARCH has been implemented. Exactly the same option names as USEARCH has been used in order to make it possible to make VSEARCH almost a drop-in replacement. The basic usearch\_global algorithm for global alignments using nucleotide sequences is implemented. At this stage it does not support amino acid sequences, local alignments, clustering, chimera detection etc. It currently does not use pthreads nor SIMD.

In the example below, VSEARCH will identify sequences in database.fsa at least 90% identical to the query sequences in queries.fsa and write the results to alnout.txt.

`./vsearch --usearch_global queries.fsa --db database.fsa --alnout alnout.txt --id 0.9 --strand plus --threads 1 --fulldp`


## Implementation details

**Algorithm:** VSEARCH indexes the unique kmers in the database in the same way as USEARCH, but is currently limited to continuous words (non-spaced seeds). It evenly samples a certain number of unique kmers from each query sequence and identifies the database sequences with the largest number of kmer matches. It then examines the database sequences in order of decreasing number of kmer matches. A full global alignment is computed and those database sequences that satisfy the sequence identity fraction criteria specified using the `--id` and `--self` options are accepted. The `--maxrejects` and `--maxaccepts` options are supported in this process, indicating the maximum number of non-matching and matching databases considered, respectively. Please see the USEARCH paper and supplementary for details.

**Kmer selection:** How many and which unique kmers USEARCH chooses from the query sequence is not well documented, but our procedure seems to give results approximately equal to USEARCH. In VSEARCH, *x* of the unique kmers in the query are sampled uniformly along the query sequence, where *x* is chosen so as to expect 8 matches with the targets satisfying the accept criteria (specified with the `--id` option). USEARCH seems to do a similar thing but the exact procedure is unknown. The choice of *x* should be looked further into.

**Output formats:** The output can be written in four different formats specified with the `--alnout`, `--blast6out`, `--uc` and the flexible format specified with the `--userout` and `--userfields`.

**Alignment:** VSEARCH currently uses a relatively slow non-vectorized full dynamic programming Needleman-Wunsch-Sellers algorithm for the alignments (similar to USEARCH with `--fulldp`) instead of the quicker (but possible less senstive) procedure involving seeding, extension and banded dynamic programming employed by USEARCH. This could be replaced by a vectorized variant using SIMD in VSEARCH to get up to comparable speed. A banded SIMD variant could also be valuable.

**Performance:** Based on very limited testing, the VSEARCH speed appears to be about half that of USEARCH when USEARCH is run with the `--fulldp` option. The accuracy also seems comparable but very variable relative to USEARCH. This needs to be looked more into.

**Command line options:** A list of all the options supported is available in the doc folder. Please see the file usearch\_options.md. The options currently supported by VSEARCH is indicated below and also in the options_supported.md


## Command line options supported

Required options:

* `--usearch_global <filename>`
* `--db <filename>`
* `--id <real>`
* `--strand <plus|both>` (Currently only plus strand supported)

Output options (at least one must be specified):
* `--alnout <filename>`
* `--blast6out <filename>`
* `--uc <filename>`
* `--userout <filename>` and `--userfields <list of fields separated by +>`

Optional options:

* `--maxaccepts <int>` (Default 1)
* `--maxrejects <int>` (Default 16)
* `--self`
* `--wordlength <int>` (Default 8)
* `--match <int>` (Default 1)
* `--mismatch <int>` (Default -2)
* `--gapopen <int>` (Default 10) (Only one common gap opening penalty supported)
* `--gapext <int>` (Default 1) (Only one common gap extension penalty supported)
* `--rowlen <int>` (Default 64)
* `--fulldp` (Default)
* `--threads <int>` (Default 1) (Currently only one thread supported)


## Main limitations

* **Commands:** Only usearch_global is supported.
* **Masking:** Currently, VSEARCH does not mask the sequences while USEARCH performs masking by default.
* **Strands:** Only the plus strand is searched.
* **Accept/reject options**: Only the `--id` and `--self` options is supported.
* **Indexing options:** Only continuous seeds are supported.
* **Gap penalties:** Only standard gap open and extension penalties are supported. Specific left/interior/right/end/query/target gap penalties are not supported.
* **Speed:** Only non-vectorized full alignment, no parallelization. Threads are currently not supported.


## License

The code is currently licensed under the GNU Affero General Public License version 3. The ordinary GNU GPL is one alternative. We should discuss this.


## Code

The code is written in C++ but most of it is actually C with some C++ syntax conventions. The files are:

* **vsearch.h** - C header file for entire project
* **vsearch.cc** - Main program file, general initialization, reads arguments and parses options, writes info.
* **query.cc** - reads the fasta file containing the query sequences.
* **db.cc** - Handles the database file read, access etc
* **dbindex.cc** - Indexes the database by identifying unique kmers in the sequences and make a database hash table
* **nw.cc** - Needleman-Wunsch-Sellers global alignment, serial
* **showalign.cc** - Output an alignment in a human-readable way given a CIGAR-string and the sequences
* **search.cc** - search database
* **util.cc** - Various common utility functions
* **results.cc** - Output results in various formats (alnout, userout, blast6, uc)
* **derepl.cc** - Code for dereplication (very rudimentary at the moment)
* **userfields.cc** - Code for parsing the userfields option argument


## Bugs

At the moment VSEARCH is not well tested. How well it works on really large databases is not checked. There might be bugs related to crossing the 4GB memory space.


## Future work

Some issues to work on:

* special gap penalties
* search both strands
* precision and recall comparison with USEARCH
* performance comparison with USEARCH
* testing and debugging
* understand how many and which unique kmers from the query USEARCH chooses
* improving the kmer match counting and target prioritization during searching
* parallelisation with pthreads
* parallelisation with SIMD-based global alignment
* clustering (i.e. uclust)
* chimera filtering
* more accept options
* masking
