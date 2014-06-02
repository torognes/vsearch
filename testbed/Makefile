CC = g++
CXXFLAGS = -Wall
INCLUDE=../src

all: test-bzquery

test-bzquery: ../src/bzquery.cc ../src/util.cc ../src/query.cc ../src/db.cc test-bzquery.cc
	$(CC) -I $(INCLUDE) -DHAVE_BZLIB_H $(CXXFLAGS) $+ -o $@ -lbz2
