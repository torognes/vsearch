bin_PROGRAMS = $(top_builddir)/bin/vsearch

AM_CFLAGS = -Wall -Wextra -Wpedantic -Wno-ignored-attributes -I$(top_srcdir)/third_party/simde

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
if TARGET_X86_64_SIMD
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
fastq_mergepairs.h \
fastq_stats.h \
fastqops.h \
fastx.h \
filter.h \
getseq.h \
kmerhash.h \
linmemalign.h \
mask.h \
md5.h \
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
utils/check_output_filehandle.hpp \
utils/cigar.hpp \
utils/cigar_operations.hpp \
utils/compare_strings_nocase.hpp \
utils/fatal.hpp \
utils/kmer_hash_struct.hpp \
utils/maps.hpp \
utils/open_file.hpp \
utils/os_byteswap.hpp \
utils/progress.hpp \
utils/seqcmp.hpp \
utils/span.hpp \
utils/taxonomic_fields.h \
utils/xpthread.hpp \
vsearch.h \
vsearch_api.h \
xstring.h

if TARGET_X86_64
libcpu_sse2_a_SOURCES = cpu.cc $(VSEARCHHEADERS)
libcpu_sse2_a_CXXFLAGS = $(AM_CXXFLAGS) -msse2
libcpu_ssse3_a_SOURCES = cpu.cc $(VSEARCHHEADERS)
libcpu_ssse3_a_CXXFLAGS = $(AM_CXXFLAGS) -mssse3 -DSSSE3
# PIC variants for static library embedding
libcpu_sse2_pic_a_SOURCES = cpu.cc $(VSEARCHHEADERS)
libcpu_sse2_pic_a_CXXFLAGS = $(AM_CXXFLAGS) -msse2 -fPIC
libcpu_ssse3_pic_a_SOURCES = cpu.cc $(VSEARCHHEADERS)
libcpu_ssse3_pic_a_CXXFLAGS = $(AM_CXXFLAGS) -mssse3 -DSSSE3 -fPIC
noinst_LIBRARIES = libcpu_sse2.a libcpu_ssse3.a libcpu_sse2_pic.a libcpu_ssse3_pic.a libcityhash.a libcityhash_pic.a libvsearch_core.a
else
libcpu_a_SOURCES = cpu.cc $(VSEARCHHEADERS)
libcpu_pic_a_SOURCES = cpu.cc $(VSEARCHHEADERS)
libcpu_pic_a_CXXFLAGS = $(AM_CXXFLAGS) -fPIC
noinst_LIBRARIES = libcpu.a libcpu_pic.a libcityhash.a libcityhash_pic.a libvsearch_core.a
endif

libcityhash_a_SOURCES = city.cc city.h
libcityhash_pic_a_SOURCES = city.cc city.h

if TARGET_WIN

libcityhash_a_CXXFLAGS = $(AM_CXXFLAGS) -Wno-sign-compare -D_MSC_VER
libcityhash_pic_a_CXXFLAGS = $(AM_CXXFLAGS) -Wno-sign-compare -D_MSC_VER -fPIC
__top_builddir__bin_vsearch_LDFLAGS = -static
__top_builddir__bin_vsearch_LDADD = libcityhash.a libcpu_ssse3.a libcpu_sse2.a
VSEARCH_PIC_LIBS = libcityhash_pic.a libcpu_ssse3_pic.a libcpu_sse2_pic.a

else

libcityhash_a_CXXFLAGS = $(AM_CXXFLAGS) -Wno-sign-compare
libcityhash_pic_a_CXXFLAGS = $(AM_CXXFLAGS) -Wno-sign-compare -fPIC

if TARGET_X86_64
__top_builddir__bin_vsearch_LDADD = libcityhash.a libcpu_ssse3.a libcpu_sse2.a
VSEARCH_PIC_LIBS = libcityhash_pic.a libcpu_ssse3_pic.a libcpu_sse2_pic.a
else
__top_builddir__bin_vsearch_LDADD = libcityhash.a libcpu.a
VSEARCH_PIC_LIBS = libcityhash_pic.a libcpu_pic.a
endif

endif

# Common source files shared between the executable and the static library.
# cpu.cc and city.cc are excluded — they need per-target SIMD/warning flags
# and are built as separate sub-libraries (libcpu_*, libcityhash_*).
VSEARCH_COMMON_SOURCES = \
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
fastq_mergepairs.cc \
fastq_stats.cc \
fastqops.cc \
fastx.cc \
filter.cc \
getseq.cc \
kmerhash.cc \
linmemalign.cc \
mask.cc \
md5.c \
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
utils/fatal.cpp \
utils/maps.cpp \
utils/check_output_filehandle.cpp \
utils/cigar.cpp \
utils/compare_strings_nocase.cpp \
utils/open_file.cpp \
utils/os_byteswap.cpp \
utils/xpthread.cpp \
utils/seqcmp.cpp \
vsearch.cc

__top_builddir__bin_vsearch_SOURCES = $(VSEARCHHEADERS) $(VSEARCH_COMMON_SOURCES)

# Static library core — same sources as the executable, compiled with -fPIC
# for embedding into shared libraries (e.g., DuckDB extensions).
libvsearch_core_a_CXXFLAGS = $(AM_CXXFLAGS) -fPIC -DVSEARCH_NO_MAIN
libvsearch_core_a_CFLAGS = $(AM_CFLAGS) -fPIC -DVSEARCH_NO_MAIN
libvsearch_core_a_SOURCES = $(VSEARCHHEADERS) $(VSEARCH_COMMON_SOURCES)

# Merged static library — combines core PIC objects with SIMD cpu and cityhash PIC objects.
# Build with: make libvsearch.a
# All objects are -fPIC so the archive can be linked into shared libraries.
VSEARCH_MERGE_TMP = .libvsearch_merge
libvsearch.a: libvsearch_core.a $(VSEARCH_PIC_LIBS)
	rm -rf $(VSEARCH_MERGE_TMP) $@
	mkdir -p $(VSEARCH_MERGE_TMP)
	cd $(VSEARCH_MERGE_TMP) && $(AR) x $(abs_builddir)/libvsearch_core.a
	for lib in $(VSEARCH_PIC_LIBS); do \
	  cd $(abs_builddir)/$(VSEARCH_MERGE_TMP) && $(AR) x $(abs_builddir)/$$lib; \
	done
	$(AR) cr $@ $(VSEARCH_MERGE_TMP)/*.o
	$(RANLIB) $@
	rm -rf $(VSEARCH_MERGE_TMP)

clean-local:
	rm -rf $(VSEARCH_MERGE_TMP) libvsearch.a
