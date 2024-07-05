bin_PROGRAMS = $(top_builddir)/bin/vsearch

AM_CFLAGS = -Wall -Wextra -Wpedantic

if TARGET_PPC
AM_CFLAGS += -mcpu=powerpc64le -maltivec
else
if TARGET_AARCH64
AM_CFLAGS += -march=armv8-a+simd -mtune=generic
else
AM_CFLAGS += -march=x86-64 -mtune=generic
endif
endif

# Conditionally set NDEBUG based on ENABLE_DEBUG
if ENABLE_DEBUG
AM_CFLAGS += -UNDEBUG -Wcast-align -Wdate-time \
             -Wdouble-promotion -Wduplicated-branches -Wduplicated-cond -Wfloat-equal \
             -Wformat=1 -Wformat-overflow -Wlogical-op -Wnon-virtual-dtor -Wnull-dereference \
             -Woverloaded-virtual -Wuninitialized \
             -Wunsafe-loop-optimizations -Wunused -Wvla
else
AM_CFLAGS += -DNDEBUG
endif

AM_CXXFLAGS = $(AM_CFLAGS) -std=c++11

export MACOSX_DEPLOYMENT_TARGET=10.9

VSEARCHHEADERS=\
align.h \
align_simd.h \
allpairs.h \
arch.h \
attributes.h \
bitmap.h \
chimera.h \
city.h \
citycrc.h \
cluster.h \
cpu.h \
cut.h \
db.h \
dbhash.h \
dbindex.h \
derep.h \
derep_prefix.h \
derepsmallmem.h \
dynlibs.h \
eestats.h \
fasta2fastq.h \
fasta.h \
fastq.h \
fastq_chars.h \
fastq_join.h \
fastqops.h \
fastx.h \
filter.h \
getseq.h \
kmerhash.h \
linmemalign.h \
maps.h \
mask.h \
md5.h \
mergepairs.h \
minheap.h \
msa.h \
orient.h \
otutable.h \
rereplicate.h \
results.h \
search.h \
searchcore.h \
search_exact.h \
sffconvert.h \
showalign.h \
sha1.h \
shuffle.h \
sintax.h \
sortbylength.h \
sortbysize.h \
subsample.h \
tax.h \
udb.h \
unique.h \
userfields.h \
util.h \
utils/seqcmp.h \
vsearch.h \
xstring.h

if TARGET_PPC
libcpu_a_SOURCES = cpu.cc $(VSEARCHHEADERS)
noinst_LIBRARIES = libcpu.a libcityhash.a
else
if TARGET_AARCH64
libcpu_a_SOURCES = cpu.cc $(VSEARCHHEADERS)
noinst_LIBRARIES = libcpu.a libcityhash.a
else
libcpu_sse2_a_SOURCES = cpu.cc $(VSEARCHHEADERS)
libcpu_sse2_a_CXXFLAGS = $(AM_CXXFLAGS) -msse2
libcpu_ssse3_a_SOURCES = cpu.cc $(VSEARCHHEADERS)
libcpu_ssse3_a_CXXFLAGS = $(AM_CXXFLAGS) -mssse3 -DSSSE3
noinst_LIBRARIES = libcpu_sse2.a libcpu_ssse3.a libcityhash.a
endif
endif

libcityhash_a_SOURCES = city.cc city.h

if TARGET_WIN

libcityhash_a_CXXFLAGS = $(AM_CXXFLAGS) -Wno-sign-compare -D_MSC_VER
__top_builddir__bin_vsearch_LDFLAGS = -static
__top_builddir__bin_vsearch_LDADD = libcityhash.a libcpu_ssse3.a libcpu_sse2.a

else

libcityhash_a_CXXFLAGS = $(AM_CXXFLAGS) -Wno-sign-compare

if TARGET_PPC
__top_builddir__bin_vsearch_LDADD = libcityhash.a libcpu.a
else
if TARGET_AARCH64
__top_builddir__bin_vsearch_LDADD = libcityhash.a libcpu.a
else
__top_builddir__bin_vsearch_LDADD = libcityhash.a libcpu_ssse3.a libcpu_sse2.a
endif
endif

endif

__top_builddir__bin_vsearch_SOURCES = $(VSEARCHHEADERS) \
align.cc \
align_simd.cc \
allpairs.cc \
arch.cc \
attributes.cc \
bitmap.cc \
chimera.cc \
cluster.cc \
cut.cc \
db.cc \
dbhash.cc \
dbindex.cc \
derep.cc \
derep_prefix.cc \
derepsmallmem.cc \
dynlibs.cc \
eestats.cc \
fasta2fastq.cc \
fasta.cc \
fastq.cc \
fastq_chars.cc \
fastq_join.cc \
fastqops.cc \
fastx.cc \
filter.cc \
getseq.cc \
kmerhash.cc \
linmemalign.cc \
maps.cc \
mask.cc \
md5.c \
mergepairs.cc \
minheap.cc \
msa.cc \
orient.cc \
otutable.cc \
rereplicate.cc \
results.cc \
search.cc \
searchcore.cc \
search_exact.cc \
sffconvert.cc \
sha1.c \
showalign.cc \
shuffle.cc \
sintax.cc \
sortbylength.cc \
sortbysize.cc \
subsample.cc \
tax.cc \
udb.cc \
unique.cc \
userfields.cc \
util.cc \
utils/seqcmp.cc \
vsearch.cc
