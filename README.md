# VSEARCH

## Introduction

The aim of this project is to create an alternative to the [USEARCH](http://www.drive5.com/usearch/) tool developed by Robert C. Edgar (2010). The new tool should:

* have open source code with an appropriate open source license
* be free of charge, gratis
* have a 64-bit design that handles very large databases and much more than 4GB of memory
* be as accurate or more accurate than usearch
* be as fast or faster than usearch

We have implemented a tool called VSEARCH which supports searching, clustering, chimera detection, dereplication, sorting and masking (commands `--usearch_global`, `--cluster_smallmem`, `--cluster_fast`, `--uchime_ref`, `--uchime_denovo`, `--derep_fulllength`, `--sortbysize`, `--sortbylength`, `--maskfasta` and `--allpairs_global`, as well as almost all their options).

VSEARCH stands for vectorized search, as the tool takes advantage of parallelism in the form of SIMD vectorization as well as multiple threads to perform accurate alignments at high speed. VSEARCH uses an optimal global aligner (full dynamic programming Needleman-Wunsch), in contrast to USEARCH which by default uses a heuristic seed and extend aligner. This results in more accurate alignments and overall improved sensitivity (recall) with VSEARCH, especially for alignments with gaps.

The same option names as in USEARCH version 7 has been used in order to make VSEARCH an almost drop-in replacement.

VSEARCH binaries are provided for x86-64 systems running GNU/Linux or OS X (10.7 or higher).

When compiled with the zlib and bzip2 libraries (as in the supplied binaries), VSEARCH can directly read input query and database files that are compressed (.gz and .bz2).

VSEARCH does not support amino acid sequences or local alignments. These features may be added in the future.

## Getting Help

If you can't find an answer in the VSEARCH documentation, please visit the [VSEARCH Web Forum](https://groups.google.com/forum/#!forum/vsearch-forum) to post a question or start a discussion.

## Example

In the example below, VSEARCH will identify sequences in the file database.fsa that are at least 90% identical on the plus strand to the query sequences in the file queries.fsa and write the results to the file alnout.txt.

`./vsearch-1.1.3-linux-x86_64 --usearch_global queries.fsa --db database.fsa --id 0.9 --alnout alnout.txt`

## Download and install

**binaries** The latest VSEARCH binaries are available [here](https://github.com/torognes/vsearch/releases) for [GNU/Linux on x86-64 systems](https://github.com/torognes/vsearch/blob/master/bin/vsearch-1.1.3-linux-x86_64) and [Apple Mac OS X on x86-64 systems](https://github.com/torognes/vsearch/blob/master/bin/vsearch-1.1.3-osx-x86_64). These executables include support for input files compressed by zlib and bzip2 (with files usually ending in `.gz` or `.bz2`). Download the appropriate executable and make a symbolic link to the vsearch binary in a folder included in your `$PATH`. You may use the following commands, assuming `~/bin` is in your `$PATH` (substitute `linux` with `osx` in those lines if you're on a Mac):

```sh
cd ~
mkdir -p bin
cd bin
wget https://github.com/torognes/vsearch/releases/download/v1.1.3/vsearch-1.1.3-linux-x86_64
chmod a+x vsearch-1.1.3-linux-x86_64
ln -s vsearch-1.1.3-linux-x86_64 vsearch
```

**source code** Use the following commands to clone the entire repository and build the executable:

```sh
git clone https://github.com/torognes/vsearch.git
cd vsearch/src
make -f Makefile
cd ../bin/
```

The alternative makefiles `Makefile.ZLIB`, `Makefile.BZLIB` and `Makefile.static` may be used to include support for compressed input files using zlib, bzip2 or both. The first two alternatives uses dynamic linking to the compression libraries, while the third uses static linking. The compression libraries [zlib](http://www.zlib.net) and/or [bzip2](http://www.bzip.org) must be downloaded and installed in folders called `zlib` and/or `bzip2`, respectively, below the main `vsearch` folder.

**documentation** The VSEARCH user's manual is available in the `doc` folder in the form of a [man page](https://github.com/torognes/vsearch/blob/master/doc/vsearch.1) and a [pdf ](https://github.com/torognes/vsearch/blob/master/doc/vsearch_manual.pdf). To install the manpage, copy `vsearch.1` file or a create a symbolic link to `vsearch.1` in a folder included in your `$MANPATH`.

**Galaxy wrapper** Thanks to the work of the [Intergalactic Utilities Commission](https://wiki.galaxyproject.org/IUC) members, vsearch is now part of the [Galaxy ToolShed](https://toolshed.g2.bx.psu.edu/view/iuc/vsearch/).

**Debian package** Thanks to the [Debian Med](https://www.debian.org/devel/debian-med/) team, there is now a [vsearch](https://packages.debian.org/sid/vsearch) package in [Debian](https://www.debian.org/). The `example/test` data is available in a separate [vsearch-data](https://packages.debian.org/sid/science/vsearch-data) package.

**Homebrew package** Thanks to [Torsten Seeman](https://github.com/tseemann), a vsearch package for [Homebrew](http://brew.sh/) [has been made](https://github.com/Homebrew/homebrew-science/pull/2409).

## Implementation details and initial assessment

**Search algorithm:** VSEARCH indexes the unique kmers in the database in a way similar to USEARCH, but is currently limited to continuous words (non-spaced seeds) of length 3-15. It samples every unique kmer from each query sequence and identifies the number of matching kmers in each database sequence. It then examines the database sequences in order of decreasing number of kmer matches. A full global alignment is computed and those target sequences that satisfy all accept options are retained while the others are rejected. The `--maxrejects` and `--maxaccepts` options are supported in this process, indicating the maximum number of non-matching and matching target sequences considered, respectively. Please see the USEARCH paper and supplementary for details.

**Kmer selection:** How many and which kmers USEARCH chooses from the query sequence is not well documented. It is also not known exactly which database sequences are examined, and in which order. We have therefore experimented with various strategies in order to obtain good performance. Our procedure seems to give results equal to or better than USEARCH.

We have chosen to select all unique kmers from the query. At least 6 of these kmers must be present in the database sequence before it will be considered. Also, at least 1 out of 16 query kmers need to be present in the database sequence. Furthermore, if several database sequences have the same number of kmer matches, they will be examined in order of decreasing sequence length.

It appears that there are differences in usearch between the searches performed by the `--usearch_global` command and the clustering commands. Notably, it appears like `--usearch_global` simply ignores the options `--wordlength`, `--slots` and `--pattern`, while the clustering commands takes them into account. VSEARCH supports the `--wordlength` option for kmer lengths from 3 to 15, but the options `--slots` and `--pattern` are ignored.

**Alignment:** VSEARCH uses a 8-way 16-bit SIMD vectorized implementation of the full dynamic programming algorithm (Needleman-Wunsch) for global sequence alignment. It is an adaptation of the method described by Rognes (2011). Due to the extreme memory requirements of this method when aligning two long sequences (e.g. more than 5000bp long), an alternative algorithm described by Hirschberg (1975) and Myers and Miller (1988) is used when aligning a pair of long sequences. This alternative algorithm uses only a linear amount of memory but is much slower. USEARCH by default uses a heuristic procedure involving seeding, extension and banded dynamic programming. If the `--fulldp` option is specified to USEARCH it will also use a full dynamic programming approach, but USEARCH is then considerably slower.

**Search Accuracy:** The accuracy of VSEARCH searches has been assessed and compared to USEARCH version 7.0.1090. The Rfam 11.0 database was used for the assessment, as described on the [USEARCH website](http://drive5.com/usearch/benchmark_rfam.html). A similar procedure was described in the USEARCH paper using the Rfam 9.1 database.

The database was initially shuffled. Then the first sequence from each of the 2085 Rfam families with at least two members was selected as queries while the rest was used as the database. The ability of VSEARCH and USEARCH to identify another member of the same family as the top hit was measured, and then recall and precision was calculated.

When USEARCH was run without the `--fulldp` option, VSEARCH had much better recall than USEARCH, but the precision was lower. The [F<sub>1</sub>-score](http://en.wikipedia.org/wiki/F1_score) was considerably higher for VSEARCH. When USEARCH was run with `--fulldp`, VSEARCH had slightly better recall, precision and F-score than USEARCH.

The recall of VSEARCH was usually about 92.3-93.5% and the precision was usually 93.0-94.1%. When run without the `--fulldp` option the recall of USEARCH was usually about 83.0-85.3% while precision was 98.5-99.0%. When run with the `--fulldp` option the recall of USEARCH was usually about 92.0-92.8% and the precision was about 92.2-93.0%.

Please see the files in the `eval` folder for the scripts used for this assessment.

**Search speed:** The speed of VSEARCH searches appears to be somewhat faster than USEARCH when USEARCH is run without the `--fulldp` option. When USEARCH is run with the `--fulldp` option, VSEARCH may be considerable faster, but it depends on the options and sequences used.

For the accuracy assessment searches in Rfam 11.0 with 100 replicates of the query sequences, VSEARCH needed 46 seconds, whereas USEARCH needed 60 seconds without the `--fulldp` option and 70 seconds with `--fulldp`. This includes time for loading, masking and indexing the database (about 2 seconds for VSEARCH, 5 seconds for USEARCH). The measurements were made on a Apple MacBook Pro Retina 2013 with four 2.3GHz Intel Core i7 cores (8 virtual cores) using the default number of threads (8).

**Memory:** VSEARCH is a 64-bit program and supports very large databases if you have enough memory. Search and clustering might use a lot of memory, especially if run with many threads. Memory usage has not been compared with USEARCH yet.

**Clustering:** The clustering commands `--cluster_smallmem` and `--cluster_fast` have been implemented. These commands support multiple threads. The only difference between `--cluster_smallmem` and `--cluster_fast` is that `--cluster_fast` will sort the sequences by length before clustering, while `--cluster_smallmem` require the sequences to be in length-sorted order unless the `--usersort` option is specified. An additional clustering command called `--cluster_size` has been added that will sort your sequences by abundance before clustering. There is no significant difference in speed or memory usage between these commands. The increased sensitivity of VSEARCH often leads to larger and fewer clusters than USEARCH.

The speed of clustering with VSEARCH relative to USEARCH depends on how many threads are used. Running with a single thread VSEARCH currently seems to be 2-4 times slower than with USEARCH, depending on parameters. Clustering has been parallelized with threads in VSEARCH, but clustering does not seem to be parallelized in USEARCH (despite what the name and documentation for `--cluster_fast` seems to indicate). Clustering with VSEARCH using 4-8 threads is often faster than USEARCH. The speed of VSEARCH might be further improved with an intra-sequence SIMD-vectorized aligner.

**Chimera detection:** Chimera detection using the algorithm described by Edgar *et al.* (2011) has been implemented in VSEARCH. Both the `--uchime_ref` and `--uchime_denovo` commands and all their options are supported.

A preliminary assessment of the accuracy of VSEARCH on chimera detection has been performed using the SIMM dataset described in the UCHIME paper. See the `eval/chimeval.sh` script and the results in `eval/chimeval.txt` for details. On the datasets with 1-5% substitutions, VSEARCH is generally on par with the original UCHIME implementation (version 4.2.40), and a bit more accurate than the implementation in USEARCH (version 7.0.1090). On the datasets with 1-5% indels, VSEARCH is clearly more accurate than both UCHIME and USEARCH.

VSEARCH is about 40% faster than USEARCH on *de novo* chimera detection and about 30% faster on detection against a reference database. In VSEARCH `uchime_ref` is multithreaded, while `uchime_denovo` is not.

**Dereplication and sorting:** The dereplication and sorting commands seems to be considerably faster in VSEARCH than in USEARCH.

**Masking:** VSEARCH by default uses an optimized multithreaded re-implementation of the well-known DUST algorithm by Tatusov and Lipman (source: ftp://ftp.ncbi.nlm.nih.gov/pub/tatusov/dust/version1/src/) to mask simple repeats and low-complexity regions in the sequences before searching and clustering. USEARCH by default uses an undocumented rapid masking method called "fastnucleo" that seems to mask fewer and smaller regions than dust. USEARCH may also be run with the DUST masking method, but the masking then takes something like 30 times longer.

**Extensions:** A shuffle command has been added. By specifying a FASTA file using the `--shuffle` option, and an output file with the `--output` option, VSEARCH will shuffle the sequences in a pseudo-random order. An integer may be specified as the seed with the `--seed` option to generate the same shuffling several times. By default, or when `--seed 0` is specified, the pseudo-random number generator will be initialized with pseudo-random data from the machine to give different numbers each time it is run.

Another extension implemented is that `--derep_fulllength` and `--cluster_fast` will honour the `--sizein` option and add together the abundances for the sequences that are clustered.

An additional clustering command called `--cluster_size` has been added that will sort sequences by abundance before clustering. 

The commands `--sortbylength` and `--sortbysize` supports the `--topn` option to output no more than the given number of sequences.

The width of FASTA formatted output files may be specified with the `--fasta_width` option and the width of alignments produced with the `--alnout` and `--uchimealn` options may be specified with the `--rowlen` and `--alignwidth` options, respectively. When an argument of zero (0) is specified for these options, sequences and alignments will not be wrapped.

VSEARCH implements the old USEARCH option `--iddef` to specify the definition of identity used to rank the hits. Values accepted are 0 (CD-HIT definition using shortest sequence as numerator), 1 (edit distance), 2 (edit distance excluding terminal gaps, default), 3 (Marine Biological Lab definition where entire gaps are considered a single difference) or 4 (BLAST, same as 2). See the [USEARCH User Guide 4.1](http://drive5.com/usearch/UsearchUserGuide4.1.pdf) page 42-44 for details. Also `id0`, `id1`, `id2`, `id3` and `id4`  are accepted as arguments to the `--userfields` option.

**Command line options:** The options currently supported by VSEARCH are indicated below. Please run VSEARCH with the `--help` option to see more information about the options.


## Command line options supported

General options:

* `--help`
* `--version`
* `--fasta_width <int>` (Default 80)
* `--maxseqlength <int>` (Default 50000)
* `--minseqlength <int>` (Default 1 for sort/shuffle or 32 for search/dereplicate)
* `--notrunclabels`
* `--threads <int>` (Default 0 means all available cores)
* `--quiet`
* `--log <filename>`

Chimera detection options:

* `--abskew <real>` (Default 2.0)
* `--alignwidth <int>` (Default 60)
* `--chimeras <filename>`
* `--db <filename>`
* `--dn <real>` (Default 1.4)
* `--mindiffs <int>` (Default 3)
* `--mindiv <real>` (Default 0.8)
* `--minh <real>` (Default 0.28)
* `--nonchimeras <filename>`
* `--self`
* `--selfid`
* `--uchime_denovo <filename>`
* `--uchime_ref <filename>`
* `--uchimealns <filename>`
* `--uchimeout <filename>`
* `--uchimeout5`
* `--xn <real>` (Default 8.0)

Clustering options (most searching options also apply):

* `--centroids <filename>`
* `--cluster_fast <filename>`
* `--cluster_size <filename>`
* `--cluster_smallmem <filename>`
* `--clusters <prefix>`
* `--consout <filename>`
* `--cons_truncate` (Ignored - Not implemented yet)
* `--id <real>` (Required)
* `--iddef <int>` (Default 2)
* `--msaout <filename>`
* `--qmask dust|none|soft` (Default dust)
* `--sizein`
* `--sizeout`
* `--strand <plus|both>` (Default plus)
* `--uc <filename>`
* `--usersort`

Dereplication options:

* `--derep_fulllength <filename>`
* `--sortbylength <filename>`
* `--sortbysize <filename>`
* `--maxsize <int>` (Default inf.)
* `--minsize <int>` (Default 0)
* `--minuniquesize <int>` (Default 1)
* `--output <filename>`
* `--relabel`
* `--sizein`
* `--sizeout`
* `--topn <int>` (Default all)

Masking options:

* `--hardmask`
* `--maskfasta <filename>`
* `--output_no_hits`
* `--qmask dust|none|soft` (Default dust)

Pairwise alignment options (most searching options also apply):

*  `--allpairs_global <filename>`
*  `--acceptall`

Searching options:

* `--alnout <filename>`
* `--blast6out <filename>`
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
* `--iddef <int>` (Default 2)
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
* `--notmatched <filename>`
* `--output_no_hits`
* `--pattern <string>` (Ignored)
* `--qmask dust|none|soft` (Default dust)
* `--query_cov <real>`
* `--rightjust`
* `--rowlen <int>` (Default 60)
* `--samout <filename>`
* `--self`
* `--selfid`
* `--sizeout`
* `--slots <int>` (Ignored)
* `--strand <plus|both>` (Default plus)
* `--target_cov <real>`
* `--top_hits_only`
* `--uc <filename>`
* `--uc_allhits`
* `--usearch_global <filename>`
* `--userfields <string>`
* `--userout <filename>`
* `--weak_id <real>`
* `--wordlength <int>` (Default 8)

Shuffling options:
* `--output <filename>`
* `--shuffle <filename>`
* `--seed <int>` (Default 0=randomize)
* `--topn <int>` (Default all)


## VSEARCH license and third party licenses

The VSEARCH code is licensed under the GNU Affero General Public License version 3.

VSEARCH includes code from Google's [CityHash project](http://code.google.com/p/cityhash/) by Geoff Pike and Jyrki Alakuijala, providing some excellent hash functions available under a MIT license.

VSEARCH includes code derived from Tatusov and Lipman's DUST program that is in the public domain.

VSEARCH binaries may include code from the [zlib](http://www.zlib.net) library copyright Jean-loup Gailly and Mark Adler, distributed under the [zlib license](http://www.zlib.net/zlib_license.html).

VSEARCH binaries may include code from the [bzip2](http://www.bzip.org) library copyright Julian R. Seward, distributed under a BSD-style license.


## Code

The code is written in C++ but most of it is actually mostly C with some C++ syntax conventions.

File | Description
---|---
**align.cc** | New Needleman-Wunsch global alignment, serial. Only for testing.
**align_simd.cc** | SIMD parallel global alignment of 1 query with 8 database sequences
**allpairs.cc** | All-vs-all optimal global pairwise alignment (no heuristics)
**arch.cc** | Architecture specific code (Mac/Linux)
**bitmap.cc** | Implementation of bitmaps
**chimera.cc** | Chimera detection
**cluster.cc** | Clustering (cluster\_fast and cluster\_smallmem)
**cpu.cc** | Code dependent on specific cpu features (e.g. ssse3)
**db.cc** | Handles the database file read, access etc
**dbindex.cc** | Indexes the database by identifying unique kmers in the sequences
**derep.cc** | Dereplication
**linmemalign.cc** | Linear memory global sequence aligner
**maps.cc** | Various character mapping arrays
**mask.cc** | Masking (DUST)
**minheap.cc** | A minheap implementation for the list of top kmer matches
**msa.cc** | Simple multiple sequence alignment and consensus sequence computation for clusters
**query.cc** | Reads the fasta file containing the query sequences.
**results.cc** | Output results in various formats (alnout, userout, blast6, uc)
**search.cc** | Implements search using global alignment
**searchcore.cc** | Core search functions for searching, clustering and chimera detection
**showalign.cc** | Output an alignment in a human-readable way given a CIGAR-string and the sequences
**shuffle.cc** | Shuffle sequences
**sortbylength.cc** | Code for sorting by length
**sortbysize.cc** | Code for sorting by size (abundance)
**string.h** | Code for a simple string class
**unique.cc** | Find unique kmers in a sequence
**userfields.cc** | Code for parsing the userfields option argument
**util.cc** | Various common utility functions
**vsearch.cc** | Main program file, general initialization, reads arguments and parses options, writes info.

VSEARCH may be compiled with zlib or bzip2 integration that allows it to read compressed FASTA files. The [zlib](http://www.zlib.net/) and the [bzip2](http://www.bzip.org/) libraries are needed for this.


## Bugs

VSEARCH has not been tested comprehensively yet. All bug reports are highly appreciated.
You may submit a bug report here on GitHub as an [issue](https://github.com/torognes/vsearch/issues),
you could post a message on the [VSEARCH Web Forum](https://groups.google.com/forum/#!forum/vsearch-forum)
or you could send an email to [torognes@ifi.uio.no](mailto:torognes@ifi.uio.no?subject=bug_in_vsearch).


## Limitations

VSEARCH is designed for rather short sequences, and will be slow when sequences are longer than about 5,000 bp. This is because it always performs optimal global alignment on selected sequences.


## Future work

Some issues to work on:

* testing and debugging
* performance evaluation
* heuristics for alignment of long sequences (e.g. banded alignment around selected diagonals)?
* intra-sequence SIMD parallelization (using the striped approach (Farrar 2007) or the plain vertical approach (Rognes & Seeberg 2000))


## The VSEARCH team

The main contributors to VSEARCH:

* Tom&aacute;&scaron; Flouri <tomas.flouri@h-its.org> (Coding, testing)
* Umer Zeeshan Ijaz <umer.ijaz@glasgow.ac.uk> (Feature suggestions)
* Fr&eacute;d&eacute;ric Mah&eacute; <mahe@rhrk.uni-kl.de> (Documentation, testing, feature suggestions)
* Ben Nichols <b.nichols.1@research.gla.ac.uk> (Evaluation)
* Christopher Quince <c.quince@warwick.ac.uk> (Initiator, feature suggestions, evaluation)
* Torbj&oslash;rn Rognes <torognes@ifi.uio.no> (Coding, testing, documentation, evaluation)


## Acknowledgements

Thanks to the following for patches and other suggestions for improvements:

* Jeff Epler <jepler@unpythonic.net>
* Andreas Tille <tille@debian.org>


## Citing VSEARCH

No papers about VSEARCH have been published yet, but a manuscript is in preparation.
For now, please cite the [VSEARCH GitHub repository](https://github.com/torognes/vsearch).
Release 1.0.16 has doi [10.5281/zenodo.15524](http://dx.doi.org/10.5281/zenodo.15524).


## Test datasets

Test datasets (found in the `data` folder) were obtained from
the [BioMarks project](http://biomarks.eu/) (Logares et al. 2014),
the [TARA OCEANS project](http://oceans.taraexpeditions.org/) (Karsenti et al. 2011) and
the [Protist Ribosomal Database](http://ssu-rrna.org/) (Guillou et al. 2012).


## References

* Edgar RC (2010)
**Search and clustering orders of magnitude faster than BLAST.**
*Bioinformatics*, 26 (19): 2460-2461.
doi:[10.1093/bioinformatics/btq461](http://dx.doi.org/10.1093/bioinformatics/btq461)

* Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R (2011)
**UCHIME improves sensitivity and speed of chimera detection.**
*Bioinformatics*, 27 (16): 2194-2200.
doi:[10.1093/bioinformatics/btr381](http://dx.doi.org/10.1093/bioinformatics/btr381)

* Farrar M (2007)
**Striped Smith-Waterman speeds database searches six times over other SIMD implementations.**
*Bioinformatics* (2007) 23 (2): 156-161.
doi:[10.1093/bioinformatics/btl582](http://dx.doi.org/10.1093/bioinformatics/btl582)

* Guillou L., Bachar D., Audic S., Bass D., Berney C., Bittner L., Boutte C., Burgaud G., de Vargas C., Decelle J., del Campo J., Dolan J., Dunthorn M., Edvardsen B., Holzmann M., Kooistra W., Lara E., Lebescot N., Logares R., Mahé F., Massana R., Montresor M., Morard R., Not F., Pawlowski J., Probert I., Sauvadet A.-L., Siano R., Stoeck T., Vaulot D., Zimmermann P. & Christen R. (2013)
**The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote Small Sub-Unit rRNA sequences with curated taxonomy.**
*Nucleic Acids Research*, 41 (D1), D597-D604.
doi:[10.1093/nar/gks1160](http://dx.doi.org/10.1093/nar/gks1160)

* Hirschberg D.S (1975) **A linear space algorithm for computing maximal common subsequences.** *Comm ACM*, 18(6), 341-343. doi:[10.1145/360825.360861](http://dx.doi.org/10.1145/360825.360861)

* Karsenti E., González Acinas S., Bork P., Bowler C., de Vargas C., Raes J., Sullivan M. B., Arendt D., Benzoni F., Claverie J.-M., Follows M., Jaillon O., Gorsky G., Hingamp P., Iudicone D., Kandels-Lewis S., Krzic U., Not F., Ogata H., Pesant S., Reynaud E. G., Sardet C., Sieracki M. E., Speich S., Velayoudon D., Weissenbach J., Wincker P. & the Tara Oceans Consortium (2011)
**A holistic approach to marine eco-systems biology.**
*PLoS Biology*, 9(10), e1001177.
doi:[10.1371/journal.pbio.1001177](http://dx.doi.org/10.1371/journal.pbio.1001177)

* Logares R., Audic S., Bass D., Bittner L., Boutte C., Christen R., Claverie J.-M., Decelle J., Dolan J. R., Dunthorn M., Edvardsen B., Gobet A., Kooistra W. H. C. F., Mahé F., Not F., Ogata H., Pawlowski J., Pernice M. C., Romac S., Shalchian-Tabrizi K., Simon N., Stoeck T., Santini S., Siano R., Wincker P., Zingone A., Richards T., de Vargas C. & Massana R. (2014) **The patterning of rare and abundant community assemblages in coastal marine-planktonic microbial eukaryotes.**
*Current Biology*, 24(8), 813-821.
doi:[10.1016/j.cub.2014.02.050](http://dx.doi.org/10.1016/j.cub.2014.02.050)

* Myers E.W., & Miller W. (1988) **Optimal alignments in linear space.**
*Comput Appl Biosci*, 4(1), 11-17.
doi:[10.1093/bioinformatics/4.1.11](http://dx.doi.org/10.1093/bioinformatics/4.1.11)

* Rognes T (2011)
**Faster Smith-Waterman database searches by inter-sequence SIMD parallelisation.**
*BMC Bioinformatics*, 12: 221.
doi:[10.1186/1471-2105-12-221](http://dx.doi.org/10.1186/1471-2105-12-221)
