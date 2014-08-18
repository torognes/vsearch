# Copyright (C) 2014 Torbjorn Rognes
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

# Profiling options
#COMMON=-pg -g
COMMON=-g

COMPILEOPT=-Wall -msse4.1 -mtune=core2 -Icityhash
#COMPILEOPT=-Weverything -msse4.1 -mtune=core2

LIBS=-lpthread
LINKFLAGS=$(COMMON)

# Intel options
#CXX=icpc
#CXXFLAGS=$(COMPILEOPT) $(COMMON) -Wno-missing-declarations -fast

# GNU options
CXX=g++
CXXFLAGS=$(COMPILEOPT) $(COMMON) -O3

PROG=vsearch

OBJS=vsearch.o query.o db.o dbindex.o search.o nw.o showalign.o util.o \
	kmercount.o userfields.o results.o sortbysize.o derep.o maps.o \
	arch.o sortbylength.o cityhash/city.o

DEPS=vsearch.h Makefile cityhash/city.h cityhash/config.h


.SUFFIXES:.o .cc

%.o : %.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

all : $(PROG)

vsearch : $(OBJS)
	$(CXX) $(LINKFLAGS) -o $@ $(OBJS) $(LIBS)

clean :
	rm -f *.o *~ $(PROG) gmon.out cityhash/*.o
