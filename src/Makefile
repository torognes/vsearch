# Makefile for vsearch
#
#   make                       release, native
#   make DEBUG=1               what --enable-debug used to give
#   make CXX=x86_64-w64-mingw32-g++    cross-compile to Windows
#
# Everything configure used to decide is decided here: the target triple
# comes from the compiler, the two compression header probes are compile
# probes, and the warning flags that older compilers reject are probed
# rather than guessed.

# The top-level Makefile passes this down; "make -C src" on its own (which is
# how the library is built) reads the same file itself.
VERSION ?= $(shell cat ../VERSION 2>/dev/null)
ifeq ($(VERSION),)
  $(error cannot read the version number from ../VERSION)
endif

# GNU make defaults CXX to "g++", which does not exist on a stock FreeBSD:
# the base compiler is clang, installed as "c++" (and "cc").  autoconf's
# AC_PROG_CXX searched a list and found it; do the same rather than fail with
# "g++: not found".  Only when the user has not named a compiler.
ifeq ($(origin CXX),default)
  ifeq (,$(shell command -v $(CXX) 2>/dev/null))
    CXX := c++
  endif
endif

##############################################################################
# 1. Target identification.  Ask the compiler, not the host: this is what
#    makes "CXX=<triple>-g++ make" a complete cross-compilation recipe, with
#    no --host= and no toolchain file.
##############################################################################

TRIPLE  := $(shell $(CXX) -dumpmachine 2>/dev/null)
MACHINE := $(firstword $(subst -, ,$(TRIPLE)))

ifneq (,$(filter x86_64 amd64, $(MACHINE)))
  ARCH := x86_64
else ifneq (,$(filter aarch64 arm64, $(MACHINE)))
  ARCH := aarch64
else ifneq (,$(filter powerpc64le ppc64le, $(MACHINE)))
  ARCH := ppc64le
else
  ARCH := generic
endif

ifneq (,$(findstring mingw, $(TRIPLE)))
  OS := windows
else ifneq (,$(findstring darwin, $(TRIPLE)))
  OS := macos
else ifneq (,$(findstring freebsd, $(TRIPLE)))
  OS := freebsd
else
  OS := linux
endif

# Everything else in the toolchain has to follow the target too, and make's
# built-in defaults ("cc", "ar") are the *host* ones.
#
# CC matters more than it looks: vsearch has two C translation units
# (vendored/md5.c and vendored/sha1.c).  With only CXX set to a cross
# compiler, they are built by the host cc, and the mingw linker then accepts
# the resulting ELF objects into a PE link without a word -- it produces a
# vsearch.exe that is quietly corrupt.  So derive CC from CXX unless the user
# named one.  The clang++ substitution has to come before the g++ one,
# because "clang++" also ends in "g++".
ifeq ($(origin CC),default)
  CC := $(patsubst %g++,%gcc,$(patsubst %clang++,%clang,$(patsubst %c++,%cc,$(CXX))))
endif

# The archiver likewise: binutils' ar cannot index COFF objects from a mingw
# cross build.  GCC knows where its own binutils are, and answers plain
# "ar"/"ranlib" for a native build.  AR has a built-in default in make, so
# test its origin rather than using ?=.
ifeq ($(origin AR),default)
  AR := $(shell $(CXX) -print-prog-name=ar)
endif
ifeq ($(origin RANLIB),undefined)
  RANLIB := $(shell $(CXX) -print-prog-name=ranlib)
endif

##############################################################################
# 2. Compiler capability probes.  configure probed these by behaviour rather
#    than by version number so that each workaround disappears by itself once
#    the supported floor moves past the compiler that needed it; keep that.
#    An unknown *positive* -W option is a hard error in GCC, while an unknown
#    -Wno-* is accepted silently and then reported once per translation unit,
#    so a -Wno- flag has to be probed through its positive spelling.
##############################################################################

# $(call cxx-accepts,FLAG) -> FLAG, or empty
cxx-accepts = $(shell printf 'int main(){return 0;}' | \
  $(CXX) -std=c++11 $(1) -fsyntax-only -x c++ - >/dev/null 2>&1 && echo '$(1)')

# GCC before 5 reports -Wmissing-field-initializers (implied by -Wextra) for
# "T value {}", which initializes every member and leaves nothing to report
# (GCC PR 36750) -- about a thousand diagnostics on 4.8/4.9.  Detected by
# behaviour, not by version number.
WARN_SUPPRESS := $(patsubst -Wignored-attributes,-Wno-ignored-attributes, \
                   $(call cxx-accepts,-Wignored-attributes))
WARN_SUPPRESS += $(shell printf 'struct A{int f;int s;};int main(){A const v{};return v.f;}' | \
  $(CXX) -std=c++11 -Wextra -Werror -fsyntax-only -x c++ - >/dev/null 2>&1 || \
  echo -Wno-missing-field-initializers)

# Debug-only warnings introduced by GCC 6 and GCC 7.
WARN_DEBUG_RECENT := $(foreach w, -Wduplicated-branches -Wduplicated-cond \
                                  -Wformat-overflow -Wnull-dereference, \
                       $(call cxx-accepts,$(w)))

##############################################################################
# 3. Compression support.  zlib and bzip2 are dlopen()ed at run time, so this
#    is a header probe and nothing links against them.  Disable with
#    "make ZLIB=0" / "make BZIP2=0" (configure spelled these --disable-zlib
#    and --disable-bzip2).
##############################################################################

ZLIB  ?= 1
BZIP2 ?= 1

# Probed with -include rather than a generated "#include" line, because a
# '#' inside a $(shell ...) call still starts a makefile comment on GNU make
# before 4.3 and fails the parse with "unterminated call to function 'shell'".
has-header = $(shell $(CXX) -fsyntax-only -include $(1) -x c++ /dev/null \
  >/dev/null 2>&1 && echo yes)

ifeq ($(ZLIB),1)
  HAVE_ZLIB := $(call has-header,zlib.h)
endif
ifeq ($(BZIP2),1)
  HAVE_BZIP2 := $(call has-header,bzlib.h)
endif

##############################################################################
# 4. Build flavours.  One of RELEASE (the default), DEBUG, PROFILE or
#    COVERAGE; SANITIZE is orthogonal and combines with any of them, exactly
#    as --enable-sanitize was independent of --enable-debug.
##############################################################################

ifeq ($(or $(RELEASE),$(DEBUG),$(PROFILE),$(COVERAGE)),)
  RELEASE := 1
endif

OPT     :=
LINKOPT :=
LIBS    :=

ifdef RELEASE
  FLAVOUR := release
  OPT += -g -O3 -DNDEBUG
else ifdef PROFILE
  FLAVOUR := profile
  OPT += -g -O1 -DNDEBUG -pg
  LINKOPT += -pg
else ifdef DEBUG
  FLAVOUR := debug
  OPT += -g -O0 -UNDEBUG -D_GLIBCXX_DEBUG
else ifdef COVERAGE
  FLAVOUR := coverage
  OPT += -g -O0 -DNDEBUG -fprofile-arcs -ftest-coverage
  LINKOPT += --coverage
  LIBS += -lgcov
endif

# Extra warnings, debug builds only.  -Wold-style-cast is verified clean on
# x86_64; the __PPC__ and __aarch64__ SIMD paths in core/align_simd.cpp keep
# some C-style casts, because AltiVec vector-value reinterpretations have no
# portable C++ named-cast equivalent.
ifdef DEBUG
  WARN_C   := -Wcast-align -Wdate-time -Wdouble-promotion -Wfloat-equal \
              -Wformat=2 -Wlogical-op -Wshadow -Wuninitialized \
              -Wunsafe-loop-optimizations -Wunused -Wunused-macros -Wvla \
              $(WARN_DEBUG_RECENT)
  # C++-only: the bundled C sources (md5.c, sha1.c) rely on many implicit
  # integer conversions and sha1.c intentionally casts away const, so these
  # stay out of the C flags.
  WARN_CXX := -Wcast-qual -Wconversion -Wnon-virtual-dtor -Wold-style-cast \
              -Woverloaded-virtual -Wsign-conversion -Wuseless-cast
endif

ifdef SANITIZE
  # The runtimes come with the native toolchain, and the Windows target links
  # -static, which AddressSanitizer does not support.
  ifeq ($(OS),windows)
    $(error SANITIZE is not supported for the Windows target)
  endif
  OPT     += -fsanitize=undefined,address -fno-omit-frame-pointer
  LINKOPT += -fsanitize=undefined,address
  FLAVOUR := $(FLAVOUR)-sanitize
endif

##############################################################################
# 5. Per-architecture and per-OS flags and sources.
##############################################################################

ifeq ($(ARCH),x86_64)
  # Not on Apple: configure's TARGET_X86_64_SIMD was "x86_64 and not apple".
  ifneq ($(OS),macos)
    ISA := -march=x86-64 -mtune=generic
  endif
  # Runtime SIMD detection via CPUID only makes sense on x86_64; elsewhere the
  # ISA level is fixed at compile time and cpu_features_detect() is a noop.
  ARCH_SOURCES := arch/x86_64/cpu_features.cpp
  # Two separate single-ISA translation units (no -DSSSE3 toggle); each is
  # compiled with just its own -m flag and cpu_features_detect() picks between
  # the resulting variants at run time.
  CPU_SOURCES  := arch/x86_64/SSE2/increment_counters.cpp \
                  arch/x86_64/SSSE3/increment_counters.cpp
else
  ARCH_SOURCES := arch/generic/cpu_features.cpp
  ifeq ($(ARCH),aarch64)
    ISA := -march=armv8-a+simd -mtune=generic
    CPU_SOURCES := arch/aarch64/increment_counters.cpp
  else ifeq ($(ARCH),ppc64le)
    ISA := -mcpu=powerpc64le -maltivec
    CPU_SOURCES := arch/ppc64le/increment_counters.cpp
  else
    CPU_SOURCES := arch/generic/increment_counters.cpp
  endif
endif

ifeq ($(OS),windows)
  OS_SOURCES := os/windows/dynlib_loader.cc os/windows/system.cc \
                os/windows/system_memory.cc os/windows/system_cpu.cc
  LIBS       += -lpsapi
  LINKOPT    += -static
  BIN        := vsearch.exe
else
  OS_SOURCES := os/posix/dynlib_loader.cc os/posix/system.cc
  BIN        := vsearch
  ifeq ($(OS),macos)
    OS_SOURCES += os/macos/system_memory.cc os/macos/system_cpu.cc
    export MACOSX_DEPLOYMENT_TARGET = 11.0
  else ifeq ($(OS),freebsd)
    OS_SOURCES += os/freebsd/system_memory.cc os/freebsd/system_cpu.cc
  else
    # os/linux doubles as the generic-Unix fallback; cgroup.cc is the
    # control-group reader its two backends share.
    OS_SOURCES += os/linux/system_memory.cc os/linux/system_cpu.cc \
                  os/linux/cgroup.cc
    LIBS       += -ldl
  endif
endif

##############################################################################
# 6. Sources.  Everything under core/, commands/, utils/ and the two root
#    translation units is compiled into both the executable and the library,
#    except the three files that belong to only one of them.  The vendored
#    sources and the per-ISA kernels need their own flags and are named
#    explicitly below.
##############################################################################

COMMON_SOURCES := \
  $(filter-out utils/fatal_exit.cpp utils/fatal_throw.cpp vsearch.cc, \
    $(wildcard core/*.cpp commands/*.cpp utils/*.cpp os/*.cpp *.cpp *.cc)) \
  vendored/city.cc vendored/md5.c vendored/sha1.c \
  $(ARCH_SOURCES) $(OS_SOURCES) $(CPU_SOURCES)

EXE_SOURCES := $(COMMON_SOURCES) vsearch.cc utils/fatal_exit.cpp
LIB_SOURCES := $(COMMON_SOURCES) utils/fatal_throw.cpp

##############################################################################
# 7. Flags and object layout.  Objects go under $(BUILDDIR), mirroring the
#    source tree, so several flavours and cross targets coexist from a single
#    checkout -- what the out-of-tree (VPATH) builds were for.
##############################################################################

# The flavour is part of the path, so a debug and a release build -- or two
# cross targets -- coexist from one checkout without a "make clean" between
# them.  $(BINDIR) is separate because every recipe, CI job and test runner
# expects the binary at bin/vsearch; a cross job that wants to keep several
# passes "make BINDIR=bin/<triple>" at the top level, which hands an
# absolute path down so that it means the same here as it does there.
BUILDDIR ?= build/$(TRIPLE)-$(FLAVOUR)
BINDIR   ?= ../bin
EXE_DIR  := $(BUILDDIR)/exe
LIB_DIR  := $(BUILDDIR)/lib

# $(CXXFLAGS)/$(CFLAGS)/$(CPPFLAGS) come last so a user assignment wins.
vs_CPPFLAGS := -I. -I$(BUILDDIR) -DHAVE_CONFIG_H
vs_COMMON   := -Wall -Wextra -Wpedantic $(WARN_SUPPRESS) $(ISA) $(OPT) -pthread
vs_CFLAGS   := $(vs_COMMON) $(WARN_C)
vs_CXXFLAGS := $(vs_COMMON) $(WARN_C) $(WARN_CXX) -std=c++11 -fno-exceptions

# The library gets -fPIC so the archive can be linked into shared libraries,
# and -fexceptions because fatal() throws a VsearchError in a library session:
# the throw must unwind cleanly through every frame, so *every* library
# translation unit needs cleanup landing pads, not just fatal_throw.cpp.
# -fexceptions comes after -fno-exceptions and wins.
LIB_EXTRA := -fPIC -fexceptions

DEPFLAGS := -MMD -MP

EXE_OBJS := $(addprefix $(EXE_DIR)/, $(addsuffix .o, $(basename $(EXE_SOURCES))))
LIB_OBJS := $(addprefix $(LIB_DIR)/, $(addsuffix .o, $(basename $(LIB_SOURCES))))
DEPFILES := $(EXE_OBJS:.o=.d) $(LIB_OBJS:.o=.d)

# Per-file flags, applied to both object trees at once.  city.cc is vendored
# upstream (Google CityHash) and its many implicit integer conversions are
# intrinsic to that code, not ours to fix.
$(EXE_DIR)/vendored/city.o $(LIB_DIR)/vendored/city.o: \
  PRIVATE := -Wno-sign-compare -Wno-conversion -Wno-sign-conversion -Ivendored
$(EXE_DIR)/arch/x86_64/SSE2/increment_counters.o \
$(LIB_DIR)/arch/x86_64/SSE2/increment_counters.o:  PRIVATE := -msse2
$(EXE_DIR)/arch/x86_64/SSSE3/increment_counters.o \
$(LIB_DIR)/arch/x86_64/SSSE3/increment_counters.o: PRIVATE := -mssse3
$(LIB_OBJS): PRIVATE += $(LIB_EXTRA)

##############################################################################
# 8. Rules.
##############################################################################

PROG   := $(BINDIR)/$(BIN)
CONFIG := $(BUILDDIR)/config.h

.PHONY: all lib clean print-bin
.DEFAULT_GOAL := all

all: $(PROG)
lib: libvsearch.a

$(PROG): $(EXE_OBJS)
	@mkdir -p $(@D)
	$(CXX) $(vs_CXXFLAGS) $(CXXFLAGS) $(LINKOPT) $(LDFLAGS) -o $@ $^ $(LIBS)

# One archive, built straight from the PIC objects.  automake needed six
# noinst_LIBRARIES and an "ar x" merge dance to get here.
libvsearch.a: $(LIB_OBJS)
	$(RM) $@
	$(AR) cr $@ $^
	$(RANLIB) $@

$(EXE_DIR)/%.o: %.cpp $(CONFIG) Makefile
	@mkdir -p $(@D)
	$(CXX) $(vs_CPPFLAGS) $(CPPFLAGS) $(vs_CXXFLAGS) $(CXXFLAGS) $(PRIVATE) $(DEPFLAGS) -c -o $@ $<
$(EXE_DIR)/%.o: %.cc $(CONFIG) Makefile
	@mkdir -p $(@D)
	$(CXX) $(vs_CPPFLAGS) $(CPPFLAGS) $(vs_CXXFLAGS) $(CXXFLAGS) $(PRIVATE) $(DEPFLAGS) -c -o $@ $<
$(EXE_DIR)/%.o: %.c $(CONFIG) Makefile
	@mkdir -p $(@D)
	$(CC) $(vs_CPPFLAGS) $(CPPFLAGS) $(vs_CFLAGS) $(CFLAGS) $(PRIVATE) $(DEPFLAGS) -c -o $@ $<

$(LIB_DIR)/%.o: %.cpp $(CONFIG) Makefile
	@mkdir -p $(@D)
	$(CXX) $(vs_CPPFLAGS) $(CPPFLAGS) $(vs_CXXFLAGS) $(CXXFLAGS) $(PRIVATE) $(DEPFLAGS) -c -o $@ $<
$(LIB_DIR)/%.o: %.cc $(CONFIG) Makefile
	@mkdir -p $(@D)
	$(CXX) $(vs_CPPFLAGS) $(CPPFLAGS) $(vs_CXXFLAGS) $(CXXFLAGS) $(PRIVATE) $(DEPFLAGS) -c -o $@ $<
$(LIB_DIR)/%.o: %.c $(CONFIG) Makefile
	@mkdir -p $(@D)
	$(CC) $(vs_CPPFLAGS) $(CPPFLAGS) $(vs_CFLAGS) $(CFLAGS) $(PRIVATE) $(DEPFLAGS) -c -o $@ $<

# config.h, the six macros that are actually load-bearing.  The feature-test
# macros were relocated here from the top of vsearch.hpp; the __STDC_* ones
# are defined empty (not 1) to match GCC's own <stdint.h> in C++ mode.
$(CONFIG): Makefile
	@mkdir -p $(@D)
	@printf '%s\n' \
	  '/* Generated by src/Makefile -- do not edit. */' \
	  '#define PACKAGE "vsearch"' \
	  '#define PACKAGE_VERSION "$(VERSION)"' \
	  '#define _GNU_SOURCE 1' \
	  '#define __STDC_CONSTANT_MACROS' \
	  '#define __STDC_FORMAT_MACROS' \
	  '#define __STDC_LIMIT_MACROS' \
	  $(if $(HAVE_ZLIB), '#define HAVE_ZLIB_H 1') \
	  $(if $(HAVE_BZIP2),'#define HAVE_BZLIB_H 1') > $@

# The one place that knows the target triple answers what the binary is
# called, so the top-level Makefile does not repeat the mingw detection.
print-bin:
	@echo $(BIN)

# "make print-<VARIABLE>" echoes one variable, so CI can assert on the flags
# without parsing a build log -- a "make -n | grep" says nothing at all once
# the tree is up to date, which would be a vacuous pass.
print-%:
	@echo '$($*)'

clean:
	$(RM) -r build libvsearch.a $(PROG)

-include $(DEPFILES)
