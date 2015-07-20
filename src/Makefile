# Copyright (C) 2014-2015 Torbjorn Rognes
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# Contact: Torbjorn Rognes <torognes@ifi.uio.no>, 
# Department of Informatics, University of Oslo, 
# PO Box 1080 Blindern, NO-0316 Oslo, Norway

# Profiling
#PROFILING=-g -fprofile-arcs -ftest-coverage
#PROFILING=-g -pg
#PROFILING=-g

# Compiler warnings
WARN=-Wall -Wsign-compare
#WARN=-Weverything

CXX=g++
CXXFLAGS=-O3 -msse2 -mtune=core2 -Icityhash $(WARN) $(PROFILING)
LINKFLAGS=$(PROFILING)
LIBS=-lpthread
TEST_LIBS=-lcheck -lm -lrt

-include tests/subdir.mk

PROG=vsearch vsearch_check

OBJS=cityhash/city.o \
	align.o align_simd.o arch.o bitmap.o chimera.o cluster.o cpu_sse2.o \
	cpu_ssse3.o db.o dbindex.o derep.o maps.o mask.o minheap.o msa.o \
	query.o results.o search.o searchcore.o showalign.o shuffle.o \
	sortbylength.o sortbysize.o unique.o userfields.o util.o vsearch.o \
	allpairs.o linmemalign.o subsample.o fastqread.o

DEPS=Makefile cityhash/city.h cityhash/config.h \
	align.h align_simd.h arch.h bitmap.h chimera.h cluster.h cpu.h db.h \
	dbindex.h derep.h maps.h mask.h minheap.h msa.h query.h results.h \
	search.h searchcore.h showalign.h shuffle.h sortbylength.h \
	sortbysize.h unique.h userfields.h util.h vsearch.h \
	allpairs.h linmemalign.h xstring.h subsample.h fastqread.h

.SUFFIXES:.o .cc

%.o : %.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

vsearch : $(OBJS)
	$(CXX) $(LINKFLAGS) -o $@ vsearch_main.c $(OBJS) $(LIBS)

vsearch_check : vsearch $(TESTS)
	$(CXX) $(LINKFLAGS) -o $@ $(OBJS) $(TESTS) $(TEST_LIBS) $(LIBS)

all : $(PROG)

clean :
	rm -f *~ $(OBJS) *.{S,gcda,gcno,gcov} $(PROG) gmon.out cityhash/*.{S,gcda,gcno,gcov}

cpu_ssse3.o : cpu.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -mssse3 -DSSSE3 -c -o $@ $<

cpu_sse2.o : cpu.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
