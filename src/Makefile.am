bin_PROGRAMS = $(top_builddir)/bin/vsearch

AM_CFLAGS = -Wall -Wextra -Wpedantic -Wno-ignored-attributes

# Conditionally set profiling based on ENABLE_PROFILING
if ENABLE_PROFILING
AM_CFLAGS += -pg -O1
endif

if TARGET_PPC
AM_CFLAGS += -mcpu=powerpc64le -maltivec
else
if TARGET_AARCH64
AM_CFLAGS += -march=armv8-a+simd -mtune=generic
else
if TARGET_X86_64
AM_CFLAGS += -march=x86-64 -mtune=generic
endif
endif
endif

# Conditionally set NDEBUG based on ENABLE_DEBUG
if ENABLE_DEBUG
# next warnings to activate:
# -fsanitize=undefined,address (run-time warnings; at first, try one at a time)
# -Wconversion
# -Wold-style-cast
# -Wshadow
# -Wsign-conversion
# -Wunused-macros
# -Wuseless-cast
# GCC 14: -Wnrvo (Named Return Value Optimization)
AM_CFLAGS += -UNDEBUG -D_GLIBCXX_DEBUG -Wcast-align -Wdate-time \
             -Wdouble-promotion -Wduplicated-branches -Wduplicated-cond -Wfloat-equal \
             -Wformat=1 -Wformat-overflow -Wlogical-op -Wnon-virtual-dtor -Wnull-dereference \
             -Woverloaded-virtual -Wuninitialized \
             -Wunsafe-loop-optimizations -Wunused -Wvla \
             -Wno-ignored-attributes
             # -Wno-ignored-attributes: ignore warnings due to
             # std::array<__m128i> (__m128i is declared with
             # alignment attributes which are ignored by the
             # std::array template
else
AM_CFLAGS += -DNDEBUG
endif

AM_CXXFLAGS = $(AM_CFLAGS) -std=c++11

export MACOSX_DEPLOYMENT_TARGET=10.9

VSEARCHHEADERS=\
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
derep_smallmem.h \
dynlibs.h \
eestats.h \
fasta2fastq.h \
fasta.h \
fastq.h \
fastq_chars.h \
fastq_join.h \
fastq_stats.h \
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
sff_convert.h \
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
utils/maps.hpp \
utils/seqcmp.h \
utils/span.hpp \
utils/taxonomic_fields.h \
vsearch.h \
xstring.h

if TARGET_X86_64
libcpu_sse2_a_SOURCES = cpu.cc $(VSEARCHHEADERS)
libcpu_sse2_a_CXXFLAGS = $(AM_CXXFLAGS) -msse2
libcpu_ssse3_a_SOURCES = cpu.cc $(VSEARCHHEADERS)
libcpu_ssse3_a_CXXFLAGS = $(AM_CXXFLAGS) -mssse3 -DSSSE3
noinst_LIBRARIES = libcpu_sse2.a libcpu_ssse3.a libcityhash.a
else
libcpu_a_SOURCES = cpu.cc $(VSEARCHHEADERS)
noinst_LIBRARIES = libcpu.a libcityhash.a
endif

libcityhash_a_SOURCES = city.cc city.h

if TARGET_WIN

libcityhash_a_CXXFLAGS = $(AM_CXXFLAGS) -Wno-sign-compare -D_MSC_VER
__top_builddir__bin_vsearch_LDFLAGS = -static
__top_builddir__bin_vsearch_LDADD = libcityhash.a libcpu_ssse3.a libcpu_sse2.a

else

libcityhash_a_CXXFLAGS = $(AM_CXXFLAGS) -Wno-sign-compare

if TARGET_X86_64
__top_builddir__bin_vsearch_LDADD = libcityhash.a libcpu_ssse3.a libcpu_sse2.a
else
__top_builddir__bin_vsearch_LDADD = libcityhash.a libcpu.a
endif

endif

__top_builddir__bin_vsearch_SOURCES = $(VSEARCHHEADERS) \
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
derep_smallmem.cc \
dynlibs.cc \
eestats.cc \
fasta2fastq.cc \
fasta.cc \
fastq.cc \
fastq_chars.cc \
fastq_join.cc \
fastq_stats.cc \
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
sff_convert.cc \
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
utils/maps.cpp \
utils/seqcmp.cc \
vsearch.cc
