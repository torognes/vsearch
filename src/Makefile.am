bin_PROGRAMS = $(top_builddir)/bin/vsearch

AM_CXXFLAGS=-Wall -Wsign-compare
AM_CFLAGS=$(AM_CXXFLAGS)

export MACOSX_DEPLOYMENT_TARGET

libcpu_sse2_a_SOURCES = cpu.cc
libcpu_sse2_a_CXXFLAGS = $(AM_CXXFLAGS) -msse2

libcpu_ssse3_a_SOURCES = cpu.cc
libcpu_ssse3_a_CXXFLAGS = $(AM_CXXFLAGS) -mssse3 -DSSSE3

libcityhash_a_SOURCES = cityhash/city.cc
libcityhash_a_CXXFLAGS = $(AM_CXXFLAGS) -Wno-sign-compare -Icityhash

noinst_LIBRARIES = libcpu_sse2.a libcpu_ssse3.a libcityhash.a

__top_builddir__bin_vsearch_LDADD = libcpu_ssse3.a libcpu_sse2.a libcityhash.a

__top_builddir__bin_vsearch_SOURCES = \
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

