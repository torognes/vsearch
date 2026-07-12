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
if TARGET_X86_64_SIMD
AM_CFLAGS += -march=x86-64 -mtune=generic
endif
endif
endif

# Conditionally set NDEBUG based on ENABLE_DEBUG
if ENABLE_DEBUG
# next warnings to activate:
# -fsanitize=undefined,address (run-time warnings; at first, try one at a time)
# GCC 14: -Wnrvo (Named Return Value Optimization)
# evaluated and left disabled:
# -Wformat=2: its only delta over -Wformat=1 is -Wformat-nonliteral,
#   which fires on the non-variadic printf wrappers fatal()/warn() and
#   cannot be silenced by __attribute__((format)) on a non-variadic
#   function
# -Wold-style-cast is a debug-build flag, and debug builds target
#   x86_64, where it is verified clean. The __PPC__ and __aarch64__
#   SIMD paths in align_simd.cc keep some C-style casts: PPC AltiVec
#   vector-value reinterpretations (in v_mask_gt, v_store and the
#   dprofile merge pipeline) have no portable C++ named-cast
#   equivalent. These paths are only compiled when cross-building for
#   POWER/ARM.
AM_CFLAGS += -UNDEBUG -D_GLIBCXX_DEBUG -Wcast-align -Wdate-time \
             -Wdouble-promotion -Wduplicated-branches -Wduplicated-cond -Wfloat-equal \
             -Wformat=1 -Wformat-overflow -Wlogical-op -Wnull-dereference \
             -Wshadow -Wuninitialized \
             -Wunsafe-loop-optimizations -Wunused -Wunused-macros -Wvla \
             -Wno-ignored-attributes
             # -Wno-ignored-attributes: ignore warnings due to
             # std::array<__m128i> (__m128i is declared with
             # alignment attributes which are ignored by the
             # std::array template
else
AM_CFLAGS += -DNDEBUG
endif

AM_CXXFLAGS = $(AM_CFLAGS) -std=c++11 -fno-exceptions

# C++-only debug warnings. Keeping these out of AM_CFLAGS avoids the
# "valid for C++/ObjC++ but not for C" diagnostics that GCC emits for
# the bundled C sources (md5.c, sha1.c).
# -Wcast-qual is also restricted here because sha1.c intentionally
# casts away const (it byte-swaps its input block in place).
# -Wconversion and -Wsign-conversion are likewise restricted to C++:
# the bundled C sources md5.c and sha1.c rely on many implicit integer
# conversions that are intrinsic to the upstream hashing code.
# -Wuseless-cast: align_simd.cc, the historical blocker (the casts were
# load-bearing in the v_load/v_store/ALIGNCORE macros), is clean now
# that those are functions. results.cc keeps one flagged cast on
# purpose: uint64_t -> size_t is a no-op on the 64-bit targets but is
# retained as an explicit, portable narrowing.
if ENABLE_DEBUG
AM_CXXFLAGS += -Wcast-qual -Wconversion -Wnon-virtual-dtor -Wold-style-cast \
               -Woverloaded-virtual -Wsign-conversion -Wuseless-cast
endif

export MACOSX_DEPLOYMENT_TARGET=10.9

VSEARCHHEADERS=\
align_simd.h \
commands/allpairs_global.hpp \
arch/x86_64/cpu_features.hpp \
attributes.h \
bitmap.h \
core/chimera.hpp \
core/chimera_internal.hpp \
commands/uchime_denovo.hpp \
commands/uchime2_denovo.hpp \
commands/uchime3_denovo.hpp \
commands/uchime_ref.hpp \
commands/chimeras_denovo.hpp \
vendored/city.h \
cli.h \
core/cluster.hpp \
core/cluster_internal.hpp \
commands/cluster_fast.hpp \
commands/cluster_smallmem.hpp \
commands/cluster_size.hpp \
commands/cluster_unoise.hpp \
cpu.h \
commands/cut.hpp \
db.h \
dbhash.h \
dbindex.h \
core/derep.hpp \
core/derep_internal.hpp \
commands/derep_fulllength.hpp \
commands/derep_id.hpp \
commands/fastx_uniques.hpp \
commands/derep_prefix.hpp \
commands/derep_smallmem.hpp \
dynlibs.h \
commands/fastq_eestats.hpp \
commands/fastq_eestats2.hpp \
core/eestats.hpp \
commands/fasta2fastq.hpp \
fasta.h \
fastq.h \
commands/fastq_chars.hpp \
commands/fastq_convert.hpp \
commands/fastq_join.hpp \
commands/fastq_mergepairs.hpp \
core/mergepairs.hpp \
core/mergepairs_internal.hpp \
commands/fastq_stats.hpp \
fastx.h \
commands/fastx_revcomp.hpp \
commands/fastx_subsample.hpp \
commands/fastx_syncpairs.hpp \
commands/fastq_filter.hpp \
commands/fastx_filter.hpp \
core/filter.hpp \
commands/fastx_getseq.hpp \
commands/fastx_getseqs.hpp \
commands/fastx_getsubseq.hpp \
core/getseq.hpp \
commands/help.hpp \
kmerhash.h \
linmemalign.h \
commands/fastx_mask.hpp \
commands/maskfasta.hpp \
core/mask.hpp \
vendored/md5.h \
minheap.h \
msa.h \
commands/orient.hpp \
otutable.h \
commands/rereplicate.hpp \
results.h \
commands/usearch_global.hpp \
core/search.hpp \
core/search_internal.hpp \
searchcore.h \
commands/search_exact.hpp \
commands/sff_convert.hpp \
showalign.h \
vendored/sha1.h \
commands/shuffle.hpp \
commands/sintax.hpp \
commands/sortbylength.hpp \
commands/sortbysize.hpp \
system.h \
tax.h \
commands/makeudb_usearch.hpp \
commands/udb2fasta.hpp \
commands/udbinfo.hpp \
commands/udbstats.hpp \
core/udb.hpp \
unique.h \
util.h \
utils/cigar.hpp \
utils/cigar_operations.hpp \
utils/compare_strings_nocase.hpp \
utils/dynlib_loader.hpp \
utils/fatal.hpp \
utils/kmer_hash_struct.hpp \
utils/logfile.hpp \
utils/maps.hpp \
utils/number_of_strands.hpp \
utils/open_file.hpp \
utils/os_byteswap.hpp \
utils/progress.hpp \
utils/seqcmp.hpp \
utils/span.hpp \
utils/taxonomic_fields.h \
utils/threads.hpp \
utils/timestamp.hpp \
utils/userfields.hpp \
utils/worker_loop.hpp \
commands/version.hpp \
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

libcityhash_a_SOURCES = vendored/city.cc vendored/city.h
libcityhash_pic_a_SOURCES = vendored/city.cc vendored/city.h

if TARGET_WIN

# -Wno-conversion -Wno-sign-conversion: city.cc is vendored upstream
# (Google CityHash) and is kept verbatim for easy resynchronisation;
# its many implicit integer conversions are not ours to fix.
libcityhash_a_CXXFLAGS = $(AM_CXXFLAGS) -Wno-sign-compare -Wno-conversion -Wno-sign-conversion -D_MSC_VER -I$(srcdir)/vendored
libcityhash_pic_a_CXXFLAGS = $(AM_CXXFLAGS) -Wno-sign-compare -Wno-conversion -Wno-sign-conversion -D_MSC_VER -fPIC -I$(srcdir)/vendored
__top_builddir__bin_vsearch_LDFLAGS = -static
__top_builddir__bin_vsearch_LDADD = libcityhash.a libcpu_ssse3.a libcpu_sse2.a
VSEARCH_PIC_LIBS = libcityhash_pic.a libcpu_ssse3_pic.a libcpu_sse2_pic.a

else

# -Wno-conversion -Wno-sign-conversion: see note in the TARGET_WIN branch above.
libcityhash_a_CXXFLAGS = $(AM_CXXFLAGS) -Wno-sign-compare -Wno-conversion -Wno-sign-conversion -I$(srcdir)/vendored
libcityhash_pic_a_CXXFLAGS = $(AM_CXXFLAGS) -Wno-sign-compare -Wno-conversion -Wno-sign-conversion -fPIC -I$(srcdir)/vendored

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
commands/allpairs_global.cpp \
attributes.cc \
bitmap.cc \
core/chimera.cpp \
commands/uchime_denovo.cpp \
commands/uchime2_denovo.cpp \
commands/uchime3_denovo.cpp \
commands/uchime_ref.cpp \
commands/chimeras_denovo.cpp \
cli.cc \
core/cluster.cpp \
commands/cluster_fast.cpp \
commands/cluster_smallmem.cpp \
commands/cluster_size.cpp \
commands/cluster_unoise.cpp \
commands/cut.cpp \
db.cc \
dbhash.cc \
dbindex.cc \
core/derep.cpp \
commands/derep_fulllength.cpp \
commands/derep_id.cpp \
commands/fastx_uniques.cpp \
commands/derep_prefix.cpp \
commands/derep_smallmem.cpp \
dynlibs.cc \
commands/fastq_eestats.cpp \
commands/fastq_eestats2.cpp \
core/eestats.cpp \
commands/fasta2fastq.cpp \
fasta.cc \
fastq.cc \
commands/fastq_chars.cpp \
commands/fastq_convert.cpp \
commands/fastq_join.cpp \
commands/fastq_mergepairs.cpp \
core/mergepairs.cpp \
commands/fastq_stats.cpp \
fastx.cc \
commands/fastx_revcomp.cpp \
commands/fastx_subsample.cpp \
commands/fastx_syncpairs.cpp \
commands/fastq_filter.cpp \
commands/fastx_filter.cpp \
core/filter.cpp \
commands/fastx_getseq.cpp \
commands/fastx_getseqs.cpp \
commands/fastx_getsubseq.cpp \
core/getseq.cpp \
commands/help.cpp \
kmerhash.cc \
linmemalign.cc \
commands/fastx_mask.cpp \
commands/maskfasta.cpp \
core/mask.cpp \
vendored/md5.c \
minheap.cc \
msa.cc \
commands/orient.cpp \
otutable.cc \
commands/rereplicate.cpp \
results.cc \
commands/usearch_global.cpp \
core/search.cpp \
searchcore.cc \
commands/search_exact.cpp \
commands/sff_convert.cpp \
vendored/sha1.c \
showalign.cc \
commands/shuffle.cpp \
commands/sintax.cpp \
commands/sortbylength.cpp \
commands/sortbysize.cpp \
tax.cc \
commands/makeudb_usearch.cpp \
commands/udb2fasta.cpp \
commands/udbinfo.cpp \
commands/udbstats.cpp \
core/udb.cpp \
unique.cc \
util.cc \
utils/fatal.cpp \
utils/maps.cpp \
utils/cigar.cpp \
utils/compare_strings_nocase.cpp \
utils/open_file.cpp \
utils/os_byteswap.cpp \
utils/seqcmp.cpp \
utils/number_of_strands.cpp \
utils/timestamp.cpp \
utils/logfile.cpp \
utils/userfields.cpp \
commands/version.cpp \
vsearch.cc

# Per-architecture sources. Runtime SIMD detection via CPUID only makes
# sense on x86_64 (elsewhere the ISA level is fixed at compile time), so
# cpu_features is compiled only there; its single caller in vsearch.cc is
# guarded by #ifdef __x86_64__ to match.
if TARGET_X86_64
VSEARCH_ARCH_SOURCES = arch/x86_64/cpu_features.cpp
else
VSEARCH_ARCH_SOURCES =
endif

# Per-OS sources. Host-OS-specific code lives in one backend per OS family,
# selected here so the common code stays free of platform #ifdefs: the runtime
# shared-library loader (dynlib::open/symbol/close) and the system_/x* helpers
# (cores, aligned alloc, file descriptors). Windows uses the Win32 API; every
# other target uses the POSIX interfaces. TARGET_WIN is set by configure when
# <windows.h> is available (native or mingw cross-build).
#
# The physical-memory queries (system_get_memused/memtotal) differ per OS and
# get their own os/<os>/system_memory.cc backend, picked below; os/linux is the
# default and also serves as the generic-Unix fallback.
if TARGET_WIN
VSEARCH_OS_SOURCES = os/windows/dynlib_loader.cc os/windows/system.cc os/windows/system_memory.cc
else
VSEARCH_OS_SOURCES = os/posix/dynlib_loader.cc os/posix/system.cc
if TARGET_APPLE
VSEARCH_OS_SOURCES += os/macos/system_memory.cc
else
if TARGET_FREEBSD
VSEARCH_OS_SOURCES += os/freebsd/system_memory.cc
else
VSEARCH_OS_SOURCES += os/linux/system_memory.cc
endif
endif
endif

__top_builddir__bin_vsearch_SOURCES = $(VSEARCHHEADERS) $(VSEARCH_COMMON_SOURCES) $(VSEARCH_ARCH_SOURCES) $(VSEARCH_OS_SOURCES)

.DEFAULT_GOAL = ../bin/vsearch

# Static library core — same sources as the executable, compiled with -fPIC
# for embedding into shared libraries (e.g., DuckDB extensions).
libvsearch_core_a_CXXFLAGS = $(AM_CXXFLAGS) -fPIC -DVSEARCH_NO_MAIN
libvsearch_core_a_CFLAGS = $(AM_CFLAGS) -fPIC -DVSEARCH_NO_MAIN
libvsearch_core_a_SOURCES = $(VSEARCHHEADERS) $(VSEARCH_COMMON_SOURCES) $(VSEARCH_ARCH_SOURCES) $(VSEARCH_OS_SOURCES)

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
