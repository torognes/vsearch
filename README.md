# VSEARCH


## Introduction

The aim of the project is to create an alternative to the USEARCH tool. The new tool should have:

* open source code with an appropriate open source license
* 64-bit design that handles very large databases and more than 4GB of memory

A tool called VSEARCH has been implemented. Exactly the same option names as USEARCH has been used in order to make it possible to make VSEARCH almost a drop-in replacement. The basic usearch_global algorithm for global alignments using nucleotide sequences is implemented. At this stage it does not support amino acid sequences, local alignments, clustering etc. No uclust/uchime/uparse/ublast. It currently does not use pthreads.

In the example below, VSEARCH will identify sequences in database.fsa at least 90% identical to the query sequences in queries.fsa and write the results to alnout.txt.

`./vsearch --usearch_global queries.fsa --db database.fsa --alnout alnout.txt --id 0.9 --strand plus --threads 1 --fulldp`


## Implementation details

**Algorithm:** VSEARCH indexes the unique kmers in the database in the same way as USEARCH, but is currently limited to continuous words (non-spaced seeds). It currently pseudo-randomly picks a certain number of unique kmers from each query sequence and identifies the database sequences with the largest number of kmer matches.  It then examines the database sequences in order of decreasing number of kmer matches. A full global alignment is computed and those that satisfy the sequence identity fraction criteria specified using the `--id` option are accepted. The `--maxrejects` and `--maxaccepts` options are supported in this process, indicating the maximum number of non-matching and matching databases considered, respectively. Please see the USEARCH paper and supplementary for details.

**Kmer selection:** How many and which unique kmers USEARCH chooses from the query sequence is not well documented, but our procedure seems to give results approximately equal to USEARCH. In VSEARCH, *x* of the unique kmers in the query are chosen pseudo-randomly; where *x* is chosen so as to expect 8 or more matches with a the targets satisfying the accept criteria. It appears that USEARCH chooses a similar number of unique kmers by sampling at equal distances along the query sequence. We should probably switch to sampling at equal distances too. The choice of *x* must be looked further into.

**Alignment:** The VSEARCH tool currently uses a serial full dynamic programming Needleman-Wunsch-Sellers algorithm for the alignments (similar to USEARCH with `--fulldp`) instead of the procedure involving seeding, extension and banded dynamic programming. This could be replaced by a parallel variant using SIMD. A banded SIMD variant could also be valuable. Pre-screening the matches using a kmer vector approach (as in SWARM) is performed and has a small positive impact in some cases.

**Performance:** The speed appears comparable to USEARCH when USEARCH is run with the `--fulldp` option, but sometimes considerably slower and sometimes considerably faster. The accuracy also seems comparable but very variable relative to USEARCH. More testing remains.

**Command line options:** I have made a list of all the options that usearch_global supports. I have indicated the options currently supported by VSEARCH in bold. Please see the file usearch_options.md


## Main limitations

* **Commands:** Only usearch_global is supported.
* **Masking:** Currently, VSEARCH does not mask the sequences while USEARCH performs masking by default.
* **Strands:** Only the plus strand is searched.
* **Threads:** Threads are currently not supported.
* **Output:** Only the `--alnout` type of output (human readable alignments) is currently supported.
* **Accept options**: Only the `--id` accept option is supported.
* **Indexing options:** Only continuous seeds are supported.
* **Gap penalties:** Only standard gap open and extension penalties are supported. Specific left/interior/right/end/query/target gap penalties are not supported.


## License

I would like the code to be licensed under the GNU Affero General Public License version 3. But we can of course discuss this.


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
* **kmervector.cc** - functions for creating and comparing kmer (qgram) vectors for sequences
* **popcount.cc** - population count implementations for various architectures


## Bugs

At the moment VSEARCH is not well tested. How well it works on really large databases is not checked. There might be bugs related to crossing the 4GB memory space.


## Future work

Some issues to work on:

* testing
* performance comparison
* improving identifying and selection of unique kmers
* improving the kmer match counting and target prioritization during searching
* parallelisation with pthreads
* parallelisation with SIMD-based global alignment
* clustering (i.e. uclust)
* more accept options
* search both strands
* masking
