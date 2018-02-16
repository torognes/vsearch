# VSEARCH

## Introduction

The aim of this project is to create an alternative to the [USEARCH](http://www.drive5.com/usearch/) tool developed by Robert C. Edgar (2010). The new tool should:

* have open source code with an appropriate open source license
* be free of charge, gratis
* have a 64-bit design that handles very large databases and much more than 4GB of memory
* be as accurate or more accurate than usearch
* be as fast or faster than usearch

We have implemented a tool called VSEARCH which supports *de novo* and reference based chimera detection, clustering, full-length and prefix dereplication, rereplication, reverse complementation, masking, all-vs-all pairwise global alignment, exact and global alignment searching, shuffling, subsampling and sorting. It also supports FASTQ file analysis, filtering, conversion and merging of paired-end reads.

VSEARCH stands for vectorized search, as the tool takes advantage of parallelism in the form of SIMD vectorization as well as multiple threads to perform accurate alignments at high speed. VSEARCH uses an optimal global aligner (full dynamic programming Needleman-Wunsch), in contrast to USEARCH which by default uses a heuristic seed and extend aligner. This usually results in more accurate alignments and overall improved sensitivity (recall) with VSEARCH, especially for alignments with gaps.

VSEARCH binaries are provided for x86-64 systems running GNU/Linux, macOS (version 10.7 or higher) and Windows (64-bit, version 7 or higher), as well as ppc64le systems running GNU/Linux.

VSEARCH can directly read input query and database files that are compressed using gzip and bzip2 (.gz and .bz2) if the zlib and bzip2 libraries are available.

Most of the nucleotide based commands and options in USEARCH version 7 are supported, as well as some in version 8. The same option names as in USEARCH version 7 has been used in order to make VSEARCH an almost drop-in replacement. VSEARCH does not support amino acid sequences or local alignments. These features may be added in the future.

## Getting Help

If you can't find an answer in the [VSEARCH documentation](https://github.com/torognes/vsearch/releases/download/v2.7.1/vsearch_manual.pdf), please visit the [VSEARCH Web Forum](https://groups.google.com/forum/#!forum/vsearch-forum) to post a question or start a discussion.

## Example

In the example below, VSEARCH will identify sequences in the file database.fsa that are at least 90% identical on the plus strand to the query sequences in the file queries.fsa and write the results to the file alnout.txt.

`./vsearch --usearch_global queries.fsa --db database.fsa --id 0.9 --alnout alnout.txt`

## Download and install

**Source distribution** To download the source distribution from a [release](https://github.com/torognes/vsearch/releases) and build the executable and the documentation, use the following commands:

```
wget https://github.com/torognes/vsearch/archive/v2.7.1.tar.gz
tar xzf v2.7.1.tar.gz
cd vsearch-2.7.1
./autogen.sh
./configure
make
make install  # as root or sudo make install
```

You may customize the installation directory using the `--prefix=DIR` option to `configure`. If the compression libraries [zlib](http://www.zlib.net) and/or [bzip2](http://www.bzip.org) are installed on the system, they will be detected automatically and support for compressed files will be included in vsearch. Support for compressed files may be disabled using the `--disable-zlib` and `--disable-bzip2` options to `configure`. A PDF version of the manual will be created from the `vsearch.1` manual file if `ps2pdf` is available, unless disabled using the `--disable-pdfman` option to `configure`. Other  options may also be applied to `configure`, please run `configure -h` to see them all. GNU autotools (version 2.63 or later) and the gcc compiler is required to build vsearch.

The IBM XL C++ compiler is recommended on ppc64le systems.

The Windows binary was compiled using the [Mingw-w64](https://mingw-w64.org/) C++ cross-compiler.

**Cloning the repo** Instead of downloading the source distribution as a compressed archive, you could clone the repo and build it as shown below. The options to `configure` as described above are still valid.

```
git clone https://github.com/torognes/vsearch.git
cd vsearch
./autogen.sh
./configure
make
make install  # as root or sudo make install
```

**Binary distribution** Starting with version 1.4.0, binary distribution files containing pre-compiled binaries as well as the documentation will be made available as part of each [release](https://github.com/torognes/vsearch/releases). The included executables include support for input files compressed by zlib and bzip2 (with files usually ending in `.gz` or `.bz2`).

Binary distributions are provided for x86-64 systems running GNU/Linux, macOS (version 10.7 or higher) and Windows (64-bit, version 7 or higher), as well as ppc64le systems running GNU/Linux.

Download the appropriate executable for your system using the following commands if you are using a Linux x86_64 system:

```sh
wget https://github.com/torognes/vsearch/releases/download/v2.7.1/vsearch-2.7.1-linux-x86_64.tar.gz
tar xzf vsearch-2.7.1-linux-x86_64.tar.gz
```

Or these commands if you are using a Linux ppc64le system:

```sh
wget https://github.com/torognes/vsearch/releases/download/v2.7.1/vsearch-2.7.1-linux-ppc64le.tar.gz
tar xzf vsearch-2.7.1-linux-ppc64le.tar.gz
```

Or these commands if you are using a Mac:

```sh
wget https://github.com/torognes/vsearch/releases/download/v2.7.1/vsearch-2.7.1-macos-x86_64.tar.gz
tar xzf vsearch-2.7.1-macos-x86_64.tar.gz
```

Or if you are using Windows, download and extract (unzip) the contents of this file:

```
https://github.com/torognes/vsearch/releases/download/v2.7.1/vsearch-2.7.1-win-x86_64.zip
```

Linux and Mac: You will now have the binary distribution in a folder called `vsearch-2.7.1-linux-x86_64` or `vsearch-2.7.1-macos-x86_64` in which you will find three subfolders `bin`, `man` and `doc`. We recommend making a copy or a symbolic link to the vsearch binary `bin/vsearch` in a folder included in your `$PATH`, and a copy or a symbolic link to the vsearch man page `man/vsearch.1` in a folder included in your `$MANPATH`. The PDF version of the manual is available in `doc/vsearch_manual.pdf`.

Windows: You will now have the binary distribution in a folder called `vsearch-2.7.1-win-x86_64`. The vsearch executable is called `vsearch.exe`. The manual in PDF format is called `vsearch_manual.pdf`.


**Documentation** The VSEARCH user's manual is available in the `man` folder in the form of a [man page](https://github.com/torognes/vsearch/blob/master/man/vsearch.1). A pdf version ([vsearch_manual.pdf](https://github.com/torognes/vsearch/releases/download/v2.7.1/vsearch_manual.pdf)) will be generated by `make`. To install the manpage manually, copy the `vsearch.1` file or a create a symbolic link to `vsearch.1` in a folder included in your `$MANPATH`. The manual in both formats is also available with the binary distribution. The manual in PDF form ([vsearch_manual.pdf](https://github.com/torognes/vsearch/releases/download/v2.7.1/vsearch_manual.pdf)) is also attached to the latest [release](https://github.com/torognes/vsearch/releases).


## Plugins, packages, and wrappers

**QIIME 2 plugin** Thanks to the [QIIME 2](https://github.com/qiime2) team, there is now a plugin called [q2-vsearch](https://github.com/qiime2/q2-vsearch) for [QIIME 2](https://qiime2.org).

**Homebrew package** Thanks to [Torsten Seeman](https://github.com/tseemann), a [vsearch package](https://github.com/Homebrew/homebrew-science/pull/2409) for [Homebrew](http://brew.sh/) has been made.

**Debian package** Thanks to the [Debian Med](https://www.debian.org/devel/debian-med/) team, there is now a [vsearch](https://packages.debian.org/sid/vsearch) package in [Debian](https://www.debian.org/).

**Galaxy wrapper** Thanks to the work of the [Intergalactic Utilities Commission](https://wiki.galaxyproject.org/IUC) members, vsearch is now part of the [Galaxy ToolShed](https://toolshed.g2.bx.psu.edu/view/iuc/vsearch/).


## Converting output to a biom file for use in QIIME and other software

With the `from-uc`command in [biom](http://biom-format.org/) 2.1.5 or later, it is possible to convert data in a `.uc` file produced by vsearch into a biom file that can be read by QIIME and other software. It is described [here](https://gist.github.com/gregcaporaso/f3c042e5eb806349fa18).

Please note that vsearch version 2.2.0 and later are able to directly output OTU tables in biom 1.0 format as well as the classic and mothur formats.


## Implementation details and initial assessment

Please see the paper for details:

Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584
doi: [10.7717/peerj.2584](https://doi.org/10.7717/peerj.2584)


## Dependencies

When compiling VSEARCH the header files for the following two optional libraries are required if support for gzip and bzip2 compressed FASTA and FASTQ input files is needed:

* libz (zlib library) (zlib.h header file) (optional)
* libbz2 (bzip2lib library) (bzlib.h header file) (optional)

On Windows these libraries are called zlib1.dll and bz2.dll.

VSEARCH will automatically check whether these libraries are available and load them dynamically.


## VSEARCH license and third party licenses

The VSEARCH code is dual-licensed either under the GNU General Public License version 3 or under the BSD 2-clause license. Please see LICENSE.txt for details.

VSEARCH includes code from several other projects. We thank the authors for making their source code available.

VSEARCH includes code from Google's [CityHash project](http://code.google.com/p/cityhash/) by Geoff Pike and Jyrki Alakuijala, providing some excellent hash functions available under a MIT license.

VSEARCH includes code derived from Tatusov and Lipman's DUST program that is in the public domain.

VSEARCH includes public domain code written by Alexander Peslyak for the MD5 message digest algorithm.

VSEARCH includes public domain code written by Steve Reid and others for the SHA1 message digest algorithm.

The VSEARCH distribution includes code from GNU Autoconf which normally is available under the GNU General Public License, but may be distributed with the special autoconf configure script exception.

VSEARCH may include code from the [zlib](http://www.zlib.net) library copyright Jean-loup Gailly and Mark Adler, distributed under the [zlib license](http://www.zlib.net/zlib_license.html).

VSEARCH may include code from the [bzip2](http://www.bzip.org) library copyright Julian R. Seward, distributed under a BSD-style license.


## Code

The code is written in C++ but most of it is actually mostly C with some C++ syntax conventions.

File | Description
---|---
**abundance.cc** | Code for extracting and printing abundance information from FASTA headers
**align.cc** | New Needleman-Wunsch global alignment, serial. Only for testing.
**align_simd.cc** | SIMD parallel global alignment of 1 query with 8 database sequences
**allpairs.cc** | All-vs-all optimal global pairwise alignment (no heuristics)
**arch.cc** | Architecture specific code (Mac/Linux)
**bitmap.cc** | Implementation of bitmaps
**chimera.cc** | Chimera detection
**city.cc** | CityHash code
**cluster.cc** | Clustering (cluster\_fast and cluster\_smallmem)
**cpu.cc** | Code dependent on specific cpu features (e.g. ssse3)
**db.cc** | Handles the database file read, access etc
**dbhash.cc** | Database hashing for exact searches
**dbindex.cc** | Indexes the database by identifying unique kmers in the sequences
**derep.cc** | Dereplication
**dynlibs.cc** | Dynamic loading of compression libraries
**eestats.cc** | Produce statistics for fastq_eestats command
**fasta.cc** | FASTA file parser
**fastq.cc** | FASTQ file parser
**fastqops.cc** | FASTQ file statistics etc
**fastx.cc** | Detection of FASTA and FASTQ files, wrapper for FASTA and FASTQ parsers
**kmerhash.cc** | Hash for kmers used by paired-end read merger
**linmemalign.cc** | Linear memory global sequence aligner
**maps.cc** | Various character mapping arrays
**mask.cc** | Masking (DUST)
**md5.c** | MD5 message digest
**mergepairs.cc** | Paired-end read merging
**minheap.cc** | A minheap implementation for the list of top kmer matches
**msa.cc** | Simple multiple sequence alignment and consensus sequence computation for clusters
**otutable.cc** | Generate OTU tables in various formats
**rerep.cc** | Rereplication
**results.cc** | Output results in various formats (alnout, userout, blast6, uc)
**search.cc** | Implements search using global alignment
**searchcore.cc** | Core search functions for searching, clustering and chimera detection
**searchexact.cc** | Exact search functions
**sha1.c** | SHA1 message digest
**showalign.cc** | Output an alignment in a human-readable way given a CIGAR-string and the sequences
**shuffle.cc** | Shuffle sequences
**sortbylength.cc** | Code for sorting by length
**sortbysize.cc** | Code for sorting by size (abundance)
**subsample.cc** | Subsampling reads from a FASTA file
**udb.cc** | UDB database file handling
**unique.cc** | Find unique kmers in a sequence
**userfields.cc** | Code for parsing the userfields option argument
**util.cc** | Various common utility functions
**vsearch.cc** | Main program file, general initialization, reads arguments and parses options, writes info.
**xstring.h** | Code for a simple string class

VSEARCH may be compiled with zlib or bzip2 integration that allows it to read compressed FASTA files. The [zlib](http://www.zlib.net/) and the [bzip2](http://www.bzip.org/) libraries are needed for this.


## Bugs

All bug reports are highly appreciated.
You may submit a bug report here on GitHub as an [issue](https://github.com/torognes/vsearch/issues),
you could post a message on the [VSEARCH Web Forum](https://groups.google.com/forum/#!forum/vsearch-forum)
or you could send an email to [torognes@ifi.uio.no](mailto:torognes@ifi.uio.no?subject=bug_in_vsearch).


## Limitations

VSEARCH is designed for rather short sequences, and will be slow when sequences are longer than about 5,000 bp. This is because it always performs optimal global alignment on selected sequences.


## Future work

Some issues to work on:

* testing and debugging
* heuristics for alignment of long sequences (e.g. banded alignment around selected diagonals)?


## The VSEARCH team

The main contributors to VSEARCH:

* Torbj&oslash;rn Rognes <torognes@ifi.uio.no> (Coding, testing, documentation, evaluation)
* Fr&eacute;d&eacute;ric Mah&eacute; <mahe@rhrk.uni-kl.de> (Documentation, testing, feature suggestions)
* Tom&aacute;&scaron; Flouri <tomas.flouri@h-its.org> (Coding, testing)
* Christopher Quince <c.quince@warwick.ac.uk> (Initiator, feature suggestions, evaluation)
* Ben Nichols <b.nichols.1@research.gla.ac.uk> (Evaluation)


## Acknowledgements

Special thanks to the following people for patches, suggestions, computer access etc:

* Davide Albanese
* Colin Brislawn
* Jeff Epler
* Christopher M. Sullivan
* Andreas Tille
* Sarah Westcott

## Citing VSEARCH

Please cite the following publication if you use VSEARCH:

Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584.
doi: [10.7717/peerj.2584](https://doi.org/10.7717/peerj.2584)

Please note that citing any of the underlying algorithms, e.g. UCHIME, may also be appropriate.

## Test datasets

Test datasets (found in the separate vsearch-data repository) were
obtained from
the [BioMarks project](http://biomarks.eu/) (Logares et al. 2014),
the [TARA OCEANS project](http://oceans.taraexpeditions.org/) (Karsenti et al. 2011)
and the [Protist Ribosomal Database](http://ssu-rrna.org/) (Guillou et al. 2012).


## References

* Edgar RC (2010)
**Search and clustering orders of magnitude faster than BLAST.**
*Bioinformatics*, 26 (19): 2460-2461.
doi:[10.1093/bioinformatics/btq461](http://dx.doi.org/10.1093/bioinformatics/btq461)

* Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R (2011)
**UCHIME improves sensitivity and speed of chimera detection.**
*Bioinformatics*, 27 (16): 2194-2200.
doi:[10.1093/bioinformatics/btr381](http://dx.doi.org/10.1093/bioinformatics/btr381)

* Guillou L, Bachar D, Audic S, Bass D, Berney C, Bittner L, Boutte C, Burgaud G, de Vargas C, Decelle J, del Campo J, Dolan J, Dunthorn M, Edvardsen B, Holzmann M, Kooistra W, Lara E, Lebescot N, Logares R, Mahé F, Massana R, Montresor M, Morard R, Not F, Pawlowski J, Probert I, Sauvadet A-L, Siano R, Stoeck T, Vaulot D, Zimmermann P & Christen R (2013)
**The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote Small Sub-Unit rRNA sequences with curated taxonomy.**
*Nucleic Acids Research*, 41 (D1), D597-D604.
doi:[10.1093/nar/gks1160](http://dx.doi.org/10.1093/nar/gks1160)

* Karsenti E, González Acinas S, Bork P, Bowler C, de Vargas C, Raes J, Sullivan M B, Arendt D, Benzoni F, Claverie J-M, Follows M, Jaillon O, Gorsky G, Hingamp P, Iudicone D, Kandels-Lewis S, Krzic U, Not F, Ogata H, Pesant S, Reynaud E G, Sardet C, Sieracki M E, Speich S, Velayoudon D, Weissenbach J, Wincker P & the Tara Oceans Consortium (2011)
**A holistic approach to marine eco-systems biology.**
*PLoS Biology*, 9(10), e1001177.
doi:[10.1371/journal.pbio.1001177](http://dx.doi.org/10.1371/journal.pbio.1001177)

* Logares R, Audic S, Bass D, Bittner L, Boutte C, Christen R, Claverie J-M, Decelle J, Dolan J R, Dunthorn M, Edvardsen B, Gobet A, Kooistra W H C F, Mahé F, Not F, Ogata H, Pawlowski J, Pernice M C, Romac S, Shalchian-Tabrizi K, Simon N, Stoeck T, Santini S, Siano R, Wincker P, Zingone A, Richards T, de Vargas C & Massana R (2014) **The patterning of rare and abundant community assemblages in coastal marine-planktonic microbial eukaryotes.**
*Current Biology*, 24(8), 813-821.
doi:[10.1016/j.cub.2014.02.050](http://dx.doi.org/10.1016/j.cub.2014.02.050)

* Rognes T (2011)
**Faster Smith-Waterman database searches by inter-sequence SIMD parallelisation.**
*BMC Bioinformatics*, 12: 221.
doi:[10.1186/1471-2105-12-221](http://dx.doi.org/10.1186/1471-2105-12-221)
