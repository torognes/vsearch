bin_PROGRAMS = $(top_builddir)/bin/vsearch

AM_CXXFLAGS=-I${srcdir}/cityhash -O3 -mtune=native -Wall -Wsign-compare -g

libcpu_sse2_a_SOURCES = cpu.cc
libcpu_sse2_a_CXXFLAGS = -msse2

libcpu_ssse3_a_SOURCES = cpu.cc
libcpu_ssse3_a_CXXFLAGS = -mssse3 -DSSSE3

noinst_LIBRARIES = libcpu_sse2.a libcpu_ssse3.a

__top_builddir__bin_vsearch_LDADD = libcpu_ssse3.a libcpu_sse2.a

__top_builddir__bin_vsearch_SOURCES = cityhash/city.h \
cityhash/config.h \
abundance.h \
align.h \
align_simd.h \
allpairs.h \
arch.h \
bitmap.h \
chimera.h \
cluster.h \
cpu.h \
db.h \
dbhash.h \
dbindex.h \
derep.h \
dynlibs.h \
fasta.h \
fastq.h \
fastqops.h \
fastx.h \
linmemalign.h \
maps.h \
mask.h \
md5.h \
minheap.h \
msa.h \
results.h \
search.h \
searchcore.h \
searchexact.h \
showalign.h \
sha1.h \
shuffle.h \
sortbylength.h \
sortbysize.h \
subsample.h \
unique.h \
userfields.h \
util.h \
vsearch.h \
xstring.h \
cityhash/city.cc \
abundance.cc \
align.cc \
align_simd.cc \
allpairs.cc \
arch.cc \
bitmap.cc \
chimera.cc \
cluster.cc \
db.cc \
dbhash.cc \
dbindex.cc \
derep.cc \
dynlibs.cc \
fasta.cc \
fastq.cc \
fastqops.cc \
fastx.cc \
linmemalign.cc \
maps.cc \
mask.cc \
md5.c \
minheap.cc \
msa.cc \
results.cc \
search.cc \
searchcore.cc \
searchexact.cc \
sha1.c \
showalign.cc \
shuffle.cc \
sortbylength.cc \
sortbysize.cc \
subsample.cc \
unique.cc \
userfields.cc \
util.cc \
vsearch.cc

