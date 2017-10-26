bin_PROGRAMS = $(top_builddir)/bin/vsearch

if TARGET_PPC
AM_CXXFLAGS=-Wall -Wsign-compare -O3 -g -mcpu=power8
else
AM_CXXFLAGS=-Wall -Wsign-compare -O3 -g
endif

AM_CFLAGS=$(AM_CXXFLAGS)

export MACOSX_DEPLOYMENT_TARGET

VSEARCHHEADERS=\
abundance.h \
align.h \
align_simd.h \
allpairs.h \
arch.h \
bitmap.h \
chimera.h \
city.h \
citycrc.h \
cluster.h \
cpu.h \
db.h \
dbhash.h \
dbindex.h \
derep.h \
dynlibs.h \
eestats.h \
fasta.h \
fastq.h \
fastqops.h \
fastx.h \
kmerhash.h \
linmemalign.h \
maps.h \
mask.h \
md5.h \
mergepairs.h \
minheap.h \
msa.h \
otutable.h \
rerep.h \
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
udb.h \
unique.h \
userfields.h \
util.h \
vsearch.h \
xstring.h

if TARGET_PPC
libcpu_a_SOURCES = cpu.cc $(VSEARCHHEADERS)
noinst_LIBRARIES = libcpu.a libcityhash.a
else
libcpu_sse2_a_SOURCES = cpu.cc $(VSEARCHHEADERS)
libcpu_sse2_a_CXXFLAGS = $(AM_CXXFLAGS) -msse2
libcpu_ssse3_a_SOURCES = cpu.cc $(VSEARCHHEADERS)
libcpu_ssse3_a_CXXFLAGS = $(AM_CXXFLAGS) -mssse3 -DSSSE3
noinst_LIBRARIES = libcpu_sse2.a libcpu_ssse3.a libcityhash.a
endif

libcityhash_a_SOURCES = city.cc city.h

if TARGET_WIN

libcityhash_a_CXXFLAGS = -Wall -Wno-sign-compare -O3 -g -D_MSC_VER
__top_builddir__bin_vsearch_LDFLAGS = -static
__top_builddir__bin_vsearch_LDADD = libregex.a libcityhash.a libcpu_ssse3.a libcpu_sse2.a

else

libcityhash_a_CXXFLAGS = -Wall -Wno-sign-compare -O3 -g

if TARGET_PPC
__top_builddir__bin_vsearch_LDADD = libcityhash.a libcpu.a
else
__top_builddir__bin_vsearch_LDADD = libcityhash.a libcpu_ssse3.a libcpu_sse2.a
endif

endif

__top_builddir__bin_vsearch_SOURCES = $(VSEARCHHEADERS) \
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
eestats.cc \
fasta.cc \
fastq.cc \
fastqops.cc \
fastx.cc \
kmerhash.cc \
linmemalign.cc \
maps.cc \
mask.cc \
md5.c \
mergepairs.cc \
minheap.cc \
msa.cc \
otutable.cc \
rerep.cc \
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
udb.cc \
unique.cc \
userfields.cc \
util.cc \
vsearch.cc
