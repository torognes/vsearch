[![Build Status](https://app.travis-ci.com/torognes/vsearch.svg?branch=master)](https://app.travis-ci.com/torognes/vsearch)

# VSEARCH

## Introduction

The aim of this project is to create an alternative to the [USEARCH](https://www.drive5.com/usearch/) tool developed by Robert C. Edgar (2010). The new tool should:

* have open source code with an appropriate open source license
* be free of charge, gratis
* have a 64-bit design that handles very large databases and much more than 4GB of memory
* be as accurate or more accurate than usearch
* be as fast or faster than usearch

We have implemented a tool called VSEARCH which supports *de novo* and reference based chimera detection, clustering, full-length and prefix dereplication, rereplication, reverse complementation, masking, all-vs-all pairwise global alignment, exact and global alignment searching, shuffling, subsampling and sorting. It also supports FASTQ file analysis, filtering, conversion and merging of paired-end reads.

VSEARCH stands for vectorized search, as the tool takes advantage of parallelism in the form of SIMD vectorization as well as multiple threads to perform accurate alignments at high speed. VSEARCH uses an optimal global aligner (full dynamic programming Needleman-Wunsch), in contrast to USEARCH which by default uses a heuristic seed and extend aligner. This usually results in more accurate alignments and overall improved sensitivity (recall) with VSEARCH, especially for alignments with gaps.

[VSEARCH binaries](https://github.com/torognes/vsearch/releases/latest) are provided for GNU/Linux on five 64-bit processor architectures: x86_64, POWER8 (ppc64le), ARMv8 (aarch64), little-endian 64-bit RISC-V (riscv64), and little-endian 64-bit MIPS (mips64el). Binaries are also provided for macOS (version 10.9 Mavericks or later) on Intel (x86_64) and Apple Silicon (ARMv8), as well as Windows (64-bit, version 7 or higher, on x86_64). VSEARCH contains native SIMD code for three processor architectures (SSE2/SSSE3, AltiVec/VMX/VSX, Neon). In addition, VSEARCH uses the SIMD Everywhere (SIMDe) library to enable building on riscv64, mips64el, and other little-endian architectures, but the performance may be lower than a native implementation.

| CPU \ OS      | GNU/Linux     | macOS  | Windows   |
| ------------- | :-----------: | :----: | :-------: |
| x86_64        |  ✔            |  ✔     |  ✔        |
| ARMv8         |  ✔            |  ✔     |           |
| POWER8        |  ✔            |        |           |
| RISC-V 64 LE  |  ✔            |        |           |
| MIPS 64 LE    |  not tested   |        |           |

Various packages, plugins and wrappers for VSEARCH are also available from other sources - see [below](https://github.com/torognes/vsearch#packages-plugins-and-wrappers).

The source code compiles correctly with `gcc` (versions 4.8.5 to 14.0)
and `llvm-clang` (3.8 to 19.0). The source code should also compile on
[FreeBSD](https://www.freebsd.org/) and
[NetBSD](https://www.netbsd.org/) systems.

VSEARCH can directly read input query and database files that are compressed using gzip (.gz) and bzip2 (.bz2) if the zlib and bzip2 libraries are available.

Most of the nucleotide based commands and options in USEARCH version 7 are supported, as well as some in version 8. The same option names as in USEARCH version 7 has been used in order to make VSEARCH an almost drop-in replacement. VSEARCH does not support amino acid sequences or local alignments. These features may be added in the future.

## Getting Help

If you can't find an answer in [online documentation](https://torognes.github.io/vsearch/), or in the [manpage](https://github.com/torognes/vsearch/releases/download/v2.30.1/vsearch_manual.pdf), please visit the [VSEARCH Web Forum](https://groups.google.com/forum/#!forum/vsearch-forum) to post a question or start a discussion.

## Example

In the example below, VSEARCH will identify sequences in the file database.fsa that are at least 90% identical on the plus strand to the query sequences in the file queries.fsa and write the results to the file alnout.txt.

`./vsearch --usearch_global queries.fsa --db database.fsa --id 0.9 --alnout alnout.txt`

## Download and install

**Source distribution** To download the source distribution from a [release](https://github.com/torognes/vsearch/releases) and build the executable and the documentation, use the following commands:

```
wget https://github.com/torognes/vsearch/archive/v2.30.1.tar.gz
tar xzf v2.30.1.tar.gz
cd vsearch-2.30.1
./autogen.sh
./configure CFLAGS="-O2" CXXFLAGS="-O2"
make ARFLAGS="cr"
sudo make install
```

You may customize the installation directory using the `--prefix=DIR` option to `configure`. If the compression libraries [zlib](https://www.zlib.net) and/or [bzip2](https://www.sourceware.org/bzip2/) are installed on the system, they will be detected automatically and support for compressed files will be included in vsearch (see section **Dependencies** below). Support for compressed files may be disabled using the `--disable-zlib` and `--disable-bzip2` options to `configure`. A PDF version of the manual will be created from the `vsearch.1` manual file if `ps2pdf` is available, unless disabled using the `--disable-pdfman` option to `configure`. It is recommended to run configure with the options `CFLAGS="-O2"` and `CXXFLAGS="-O2"`. Other  options may also be applied to `configure`, please run `configure -h` to see them all. GNU autoconf (version 2.63 or later), automake and the GCC C++ (`g++`) compiler is required to build vsearch. Version 3.82 or later of `make` may be required on Linux, while version 3.81 is sufficient on macOS.

Warning: Compiling the `align_simd.cc` file on x86_64 systems using the GNU C++ compiler version 9 or later with the `-O3` optimization option on may result in incorrect code that may cause bad alignments in some circumstances. This was due to the `-ftree-partial-pre` optimization enabled by `-O3`. A compiler pragma has been inserted in the code to specifically turn off this optimization for the affected code. Using `-O3` should be safe.

To build VSEARCH on Debian and similar Linux distributions (Ubuntu etc) you'll need the following packages: autoconf, automake, g++, ghostscript, groff, libbz2-dev, make, zlib1g-dev. Include libsimde-dev to build on riscv64 or mips64el.

To build VSEARCH on Fedora and similar Linux distributions (RHEL, Centos etc) you'll need the following packages: autoconf, automake, bzip2-devel, gcc-c++, ghostscript, groff-base, make, zlib-devel.

Instead of downloading the source distribution as a compressed archive, you could clone the repo and build it as shown below. The options to `configure` as described above are still valid.

```
git clone https://github.com/torognes/vsearch.git
cd vsearch
./autogen.sh
./configure CFLAGS="-O2" CXXFLAGS="-O2"
make ARFLAGS="cr"
sudo make install
```

**Binary distribution**: Starting with version 1.4.0, binary distribution files containing pre-compiled binaries as well as the documentation will be made available as part of each [release](https://github.com/torognes/vsearch/releases). The included executables include support for input files compressed by zlib and bzip2 (with files usually ending in `.gz` or `.bz2`).

Binary distributions are provided for x86-64 systems running GNU/Linux, macOS (version 10.7 or higher) or Windows (64-bit, version 7 or higher), 64-bit AMDv8 (aarch64) systems running GNU/Linux or macOS, as well as POWER8 (ppc64le), 64-bit little-endian RISC-V (risv64), and 64-bit little endian MIPS (mips64el) systems running GNU/Linux. A universal macOS binary is also provided. In addition, an x86_64 binary built for the discontinued RHEL 7 and CentOS 7 linux distributions is provided. The other Linux binaries are built on Debian 11 (oldstable, Bullseye). Static binaries are available for all Linux architectures except x86_64, these can be used on systems that do not have all the necessary libraries installed. The Windows binary was built with cross compilation using [Mingw-w64](http://mingw-w64.org/).

Download the appropriate executable for your system using the following commands if you are using a Linux or macOS system:

```sh
wget https://github.com/torognes/vsearch/releases/download/v{VERSION}/vsearch-{VERSION}-{OS}-{ARCH}.tar.gz
tar xzf vsearch-{VERSION}-{OS}-{ARCH}.tar.gz
```

Replace `{VERSION}` with the VSEARCH version number (e.g. `2.30.1`), `{OS}` with the target operating system (`linux` or `macos`), and `{ARCH}` with the architecture (`x86_64`, `aarch64`, `ppc64le`, `riscv64`, or `mips64el`). You could add `-static` after `{ARCH}` to get a statically compiled version for Linux (except x86_64). The name of the binary for the RHEL 7 and CentOS 7 Linux distributions ends in `-ubi7`.

Or, if you are using Windows, download and extract (unzip) the contents of this file:

```
https://github.com/torognes/vsearch/releases/download/v{VERSION}/vsearch-{VERSION}-win-x86_64.zip
```

**Linux and Mac**: You will now have the binary distribution in a folder called `vsearch-{VERSION}-{OS}-{ARCH}` in which you will find three subfolders `bin`, `man` and `doc`. We recommend making a copy or a symbolic link to the vsearch binary `bin/vsearch` in a folder included in your `$PATH`, and a copy or a symbolic link to the vsearch man page `man/vsearch.1` in a folder included in your `$MANPATH`. The PDF version of the manual is available in `doc/vsearch_manual.pdf`.

**Windows**: You will now have the binary distribution in a folder
called `vsearch-{VERSION}-win-x86_64`. The vsearch executable is called
`vsearch.exe`. The manual in PDF format is called
`vsearch_manual.pdf`. If you want to be able to call `vsearch.exe`
from any command prompt window, you can put the VSEARCH executable in
a folder (for instance `C:\Users\<yourname>\bin`), and add the new
folder to the user `Path`: open the `Environment Variables` window by
searching for it in the Start menu, `Edit` user variables, add
`;C:\Users\<yourname>\bin` to the end of the `Path` variable, and save
your changes. The windows distribution also includes the `libbz2.dll`
and `zlib1.dll` files required for reading compressed input
files. These DLL's have been obtained for mingw-w64 from the MSYS2
platform.

**Documentation:** The VSEARCH user's manual is available in the `man` folder in the form of a [man page](https://github.com/torognes/vsearch/blob/master/man/vsearch.1). A pdf version ([vsearch_manual.pdf](https://github.com/torognes/vsearch/releases/download/v2.30.1/vsearch_manual.pdf)) will be generated by `make`. To install the manpage manually, copy the `vsearch.1` file or a create a symbolic link to `vsearch.1` in a folder included in your `$MANPATH`. The manual in both formats is also available with the binary distribution. The manual in PDF form ([vsearch_manual.pdf](https://github.com/torognes/vsearch/releases/download/v2.30.1/vsearch_manual.pdf)) is also attached to the latest [release](https://github.com/torognes/vsearch/releases).


## Packages, plugins, and wrappers

**Conda package** Thanks to the [BioConda](https://bioconda.github.io/) team, there is now a [vsearch package](https://anaconda.org/bioconda/vsearch) in [Conda](https://conda.io/).

**Debian package** Thanks to the [Debian Med](https://www.debian.org/devel/debian-med/) team, there is now a [vsearch](https://packages.debian.org/sid/vsearch) package in [Debian](https://www.debian.org/).

**FreeBSD ports package** Thanks to [Jason Bacon](https://github.com/outpaddling), a [vsearch](https://www.freebsd.org/cgi/ports.cgi?query=vsearch&stype=all) [FreeBSD ports](https://www.freebsd.org/ports/) package is available. Install the binary package with `pkg install vsearch`, or build from source with additional optimizations.

**Galaxy wrapper** Thanks to the work of the [Intergalactic Utilities Commission](https://wiki.galaxyproject.org/IUC) members, VSEARCH is now part of the [Galaxy ToolShed](https://toolshed.g2.bx.psu.edu/view/iuc/vsearch/).

**Homebrew package** Thanks to [Torsten Seeman](https://github.com/tseemann), a [vsearch package](https://formulae.brew.sh/formula/vsearch) for [Homebrew](http://brew.sh/) has been made.

**Pkgsrc package** Thanks to [Jason Bacon](https://github.com/outpaddling), a vsearch [pkgsrc](https://www.pkgsrc.org) package is available for NetBSD and other UNIX-like systems. Install the binary package with `pkgin install vsearch`, or build from source with additional optimizations.

**QIIME 2 plugin** Thanks to the [QIIME 2](https://github.com/qiime2) team, there is now a plugin called [q2-vsearch](https://github.com/qiime2/q2-vsearch) for [QIIME 2](https://qiime2.org).


## Converting output to a biom file for use in QIIME and other software

With the `from-uc`command in [biom](http://biom-format.org/) 2.1.5 or later, it is possible to convert data in a `.uc` file produced by vsearch into a biom file that can be read by QIIME and other software. It is described [here](https://gist.github.com/gregcaporaso/f3c042e5eb806349fa18).

Please note that VSEARCH version 2.2.0 and later are able to directly output OTU tables in biom 1.0 format as well as the classic and mothur formats.


## Implementation details and initial assessment

Please see the paper for details:

Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584
doi: [10.7717/peerj.2584](https://doi.org/10.7717/peerj.2584)


## Dependencies

Compiling VSEARCH requires either GCC (`g++`) or `clang`, `make` and the autotools (`ui-auto` on Debian-based distributions). Optionally, the header files for the following two optional libraries are required if support for gzip and bzip2 compressed FASTA and FASTQ input files is needed:

* libz (zlib library) (`zlib.h` header file, available as `zlib1g-dev` on Debian-based distributions) (optional)
* libbz2 (bzip2lib library) (`bzlib.h` header file, available as `libbz2-dev`on Debian-based distributions) (optional)

VSEARCH will automatically check whether these libraries are available and load them dynamically.

On Windows these libraries are called `zlib1.dll` and `libbz2.dll`. These DLL's are included with the released distribution of vsearch 2.27.0 and later.

To create the PDF file with the manual the ps2pdf tool is required. It is part of the `ghostscript` package.


## VSEARCH license and third party licenses

The VSEARCH code is dual-licensed either under the GNU General Public License version 3 or under the BSD 2-clause license. Please see LICENSE.txt for details.

VSEARCH includes code from several other projects. We thank the authors for making their source code available.

VSEARCH includes code from Google's [CityHash project](https://github.com/google/cityhash) by Geoff Pike and Jyrki Alakuijala, providing some excellent hash functions available under a MIT license.

VSEARCH includes code derived from Tatusov and Lipman's DUST program that is in the public domain.

VSEARCH includes public domain code written by Alexander Peslyak for the MD5 message digest algorithm.

VSEARCH includes public domain code written by Steve Reid and others for the SHA1 message digest algorithm.

The VSEARCH distribution includes code from GNU Autoconf which normally is available under the GNU General Public License, but may be distributed with the special autoconf configure script exception.

VSEARCH may include code from the [zlib](https://www.zlib.net) library copyright Jean-loup Gailly and Mark Adler, distributed under the [zlib license](https://www.zlib.net/zlib_license.html).

VSEARCH may include code from the [bzip2](https://www.sourceware.org/bzip2/) library copyright Julian R. Seward, distributed under a BSD-style license.


## Code

The code is written mostly in C++.

File | Description
---|---
**align_simd.cc** | SIMD parallel global alignment of 1 query with 8 database sequences
**allpairs.cc** | All-vs-all optimal global pairwise alignment (no heuristics)
**arch.cc** | Architecture specific code (Mac/Linux)
**attributes.cc** | Extraction and printing of attributes in FASTA headers
**bitmap.cc** | Implementation of bitmaps
**chimera.cc** | Chimera detection
**city.cc** | CityHash code
**cluster.cc** | Clustering (cluster\_fast and cluster\_smallmem)
**cpu.cc** | Code dependent on specific cpu features (e.g. ssse3)
**cut.cc** | Restriction site cutting
**db.cc** | Handles the database file read, access etc
**dbhash.cc** | Database hashing for exact searches
**dbindex.cc** | Indexes the database by identifying unique kmers in the sequences
**derep.cc** | Dereplication, full-length
**derep_prefix.cc** | Dereplication, prefix
**derep_smallmem.cc** | Dereplication, small memory usage
**dynlibs.cc** | Dynamic loading of compression libraries
**eestats.cc** | Produce statistics for fastq_eestats command
**fasta.cc** | FASTA file parser
**fasta2fastq.cc** | FASTA to FASTQ conversion
**fastq.cc** | FASTQ file parser
**fastq_chars.cc** | FASTQ statistics
**fastq_join.cc** | FASTQ paired-end reads joining
**fastqops.cc** | FASTQ file statistics etc
**fastx.cc** | Detection of FASTA and FASTQ files, wrapper for FASTA and FASTQ parsers
**filter.cc** | Trimming and filtering of sequences in FASTA and FASTQ files
**getseq.cc** | Extraction of sequences based on header labels
**kmerhash.cc** | Hash for kmers used by paired-end read merger
**linmemalign.cc** | Linear memory global sequence aligner
**mask.cc** | Masking (DUST)
**md5.c** | MD5 message digest
**mergepairs.cc** | Paired-end read merging
**minheap.cc** | A minheap implementation for the list of top kmer matches
**msa.cc** | Simple multiple sequence alignment and consensus sequence computation for clusters
**orient.cc** | Orient direction of sequences based on reference database
**otutable.cc** | Generate OTU tables in various formats
**rereplicate.cc** | Rereplication
**results.cc** | Output results in various formats (alnout, userout, blast6, uc)
**search.cc** | Implements search using global alignment
**search_exact.cc** | Exact search functions
**searchcore.cc** | Core search functions for searching, clustering and chimera detection
**sff_convert.cc** | SFF to FASTQ file conversion
**sha1.c** | SHA1 message digest
**showalign.cc** | Output an alignment in a human-readable way given a CIGAR-string and the sequences
**shuffle.cc** | Shuffle sequences
**sintax.cc** | Taxonomic classification using Sintax method
**sortbylength.cc** | Code for sorting by length
**sortbysize.cc** | Code for sorting by size (abundance)
**subsample.cc** | Subsampling reads from a FASTA file
**tax.cc** | Taxonomy information parsing
**udb.cc** | UDB database file handling
**unique.cc** | Find unique kmers in a sequence
**userfields.cc** | Code for parsing the userfields option argument
**util.cc** | Various common utility functions
**vsearch.cc** | Main program file, general initialization, reads arguments and parses options, writes info.
**utils/maps.cc** | Utilities, maps for encoding of nucleotides
**utils/seqcmp.cpp** | Utilities, sequence comparison

VSEARCH may be compiled with zlib or bzip2 integration that allows it to read compressed FASTA files. The [zlib](http://www.zlib.net/) and the [bzip2](https://www.sourceware.org/bzip2/) libraries are needed for this.


## Bugs

All bug reports are highly appreciated.
You may submit a bug report here on GitHub as an [issue](https://github.com/torognes/vsearch/issues) (preferred),
you could post a message on the [VSEARCH Web Forum](https://groups.google.com/forum/#!forum/vsearch-forum)
or you could send an email to [torognes@ifi.uio.no](mailto:torognes@ifi.uio.no?subject=bug_in_vsearch).


## Limitations

VSEARCH is designed for rather short sequences, and will be slow when sequences are longer than about 5,000 bp. This is because it always performs optimal global alignment on selected sequences.


## The VSEARCH team

The main contributors to VSEARCH:

* Torbj&oslash;rn Rognes <torognes@ifi.uio.no> (Coding, testing, documentation, evaluation)
* Fr&eacute;d&eacute;ric Mah&eacute; <mahe@rhrk.uni-kl.de> (Documentation, testing, coding, feature suggestions)
* Tom&aacute;&scaron; Flouri <tomas.flouri@h-its.org> (Coding, testing)
* Christopher Quince <c.quince@warwick.ac.uk> (Initiator, feature suggestions, evaluation)
* Ben Nichols <b.nichols.1@research.gla.ac.uk> (Evaluation)


## Acknowledgements

Special thanks to the following people for patches, suggestions, computer access etc:

* Davide Albanese
* Colin Brislawn
* Michael R. Crusoe
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
the BioMarks project ([Logares et al. 2014](https://doi.org/10.1016/j.cub.2014.02.050)),
the [TARA OCEANS project](https://oceans.taraexpeditions.org/en/) ([Karsenti et al. 2011](https://doi.org/10.1371/journal.pbio.1001177))
and the [Protist Ribosomal Reference Database (PR<sup>2</sup>)](https://github.com/pr2database/pr2database) ([Guillou et al. 2013](https://doi.org/10.1093/nar/gks1160)).


## References

* Edgar RC (2010)
**Search and clustering orders of magnitude faster than BLAST.**
*Bioinformatics*, 26 (19): 2460-2461.
doi:[10.1093/bioinformatics/btq461](https://doi.org/10.1093/bioinformatics/btq461)

* Edgar RC (2016)
**SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS sequences.**
*bioRxiv*.
doi:[10.1101/074161](https://doi.org/10.1101/074161)

* Edgar RC (2016)
**UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing.**
*bioRxiv*.
doi:[10.1101/081257](https://doi.org/10.1101/081257)

* Edgar RC, Flyvbjerg H (2015)
**Error filtering, pair assembly and error correction for next-generation sequencing reads.**
*Bioinformatics*, 31 (21): 3476-3482.
doi:[10.1093/bioinformatics/btv401](https://doi.org/10.1093/bioinformatics/btv401)

* Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R (2011)
**UCHIME improves sensitivity and speed of chimera detection.**
*Bioinformatics*, 27 (16): 2194-2200.
doi:[10.1093/bioinformatics/btr381](https://doi.org/10.1093/bioinformatics/btr381)

* Guillou L, Bachar D, Audic S, Bass D, Berney C, Bittner L, Boutte C, Burgaud G, de Vargas C, Decelle J, del Campo J, Dolan J, Dunthorn M, Edvardsen B, Holzmann M, Kooistra W, Lara E, Lebescot N, Logares R, Mahé F, Massana R, Montresor M, Morard R, Not F, Pawlowski J, Probert I, Sauvadet A-L, Siano R, Stoeck T, Vaulot D, Zimmermann P & Christen R (2013)
**The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote Small Sub-Unit rRNA sequences with curated taxonomy.**
*Nucleic Acids Research*, 41 (D1), D597-D604.
doi:[10.1093/nar/gks1160](https://doi.org/10.1093/nar/gks1160)

* Karsenti E, González Acinas S, Bork P, Bowler C, de Vargas C, Raes J, Sullivan M B, Arendt D, Benzoni F, Claverie J-M, Follows M, Jaillon O, Gorsky G, Hingamp P, Iudicone D, Kandels-Lewis S, Krzic U, Not F, Ogata H, Pesant S, Reynaud E G, Sardet C, Sieracki M E, Speich S, Velayoudon D, Weissenbach J, Wincker P & the Tara Oceans Consortium (2011)
**A holistic approach to marine eco-systems biology.**
*PLoS Biology*, 9(10), e1001177.
doi:[10.1371/journal.pbio.1001177](https://doi.org/10.1371/journal.pbio.1001177)

* Logares R, Audic S, Bass D, Bittner L, Boutte C, Christen R, Claverie J-M, Decelle J, Dolan J R, Dunthorn M, Edvardsen B, Gobet A, Kooistra W H C F, Mahé F, Not F, Ogata H, Pawlowski J, Pernice M C, Romac S, Shalchian-Tabrizi K, Simon N, Stoeck T, Santini S, Siano R, Wincker P, Zingone A, Richards T, de Vargas C & Massana R (2014) **The patterning of rare and abundant community assemblages in coastal marine-planktonic microbial eukaryotes.**
*Current Biology*, 24(8), 813-821.
doi:[10.1016/j.cub.2014.02.050](https://doi.org/10.1016/j.cub.2014.02.050)

* Rognes T (2011)
**Faster Smith-Waterman database searches by inter-sequence SIMD parallelisation.**
*BMC Bioinformatics*, 12: 221.
doi:[10.1186/1471-2105-12-221](https://doi.org/10.1186/1471-2105-12-221)
