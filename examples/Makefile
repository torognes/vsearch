# Build examples against libvsearch.a
#
# Usage:
#   make                    # build all examples
#   make test               # build and run all examples, compare with vsearch
#   make clean              # remove binaries
#
# Prerequisites:
#   - libvsearch.a built: make -C ../src libvsearch.a
#   - vsearch binary built: make -C ../src (for ground truth comparison)

VSEARCH_SRC = ../src
VSEARCH_LIB = $(VSEARCH_SRC)/libvsearch.a
VSEARCH_BIN = ../bin/vsearch

CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall -I$(VSEARCH_SRC)
LDFLAGS = $(VSEARCH_LIB) -lpthread -ldl

SHELL = /bin/bash

EXAMPLES = example_chimera example_search example_merge example_cluster example_dust example_derep

all: $(EXAMPLES)

example_chimera: example_chimera.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

example_search: example_search.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

example_merge: example_merge.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

example_cluster: example_cluster.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

example_dust: example_dust.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

example_derep: example_derep.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

test: $(EXAMPLES)
	@echo "=== Running example_chimera ==="
	./example_chimera > /tmp/test_chimera.tsv
	@diff <(sort /tmp/test_chimera.tsv) <(sort data/expected_chimera.tsv) \
		&& echo "PASS: chimera output matches vsearch" \
		|| (echo "FAIL: chimera output differs"; diff /tmp/test_chimera.tsv data/expected_chimera.tsv | head -20; exit 1)
	@echo ""
	@echo "=== Running example_search ==="
	./example_search > /tmp/test_search.tsv
	@diff <(sort /tmp/test_search.tsv) <(sort data/expected_search.tsv) \
		&& echo "PASS: search output matches vsearch" \
		|| (echo "FAIL: search output differs"; diff /tmp/test_search.tsv data/expected_search.tsv | head -20; exit 1)
	@echo ""
	@echo "=== Running example_cluster ==="
	./example_cluster > /tmp/test_cluster.uc
	@diff <(grep "^[SH]" /tmp/test_cluster.uc | sort) <(grep "^[SH]" data/expected_cluster.uc | sort) \
		&& echo "PASS: cluster output matches vsearch" \
		|| (echo "FAIL: cluster output differs"; exit 1)
	@echo ""
	@echo "=== Running example_merge ==="
	./example_merge > /tmp/test_merge.fasta
	@diff /tmp/test_merge.fasta data/expected_merge.fasta \
		&& echo "PASS: merge output matches vsearch" \
		|| (echo "FAIL: merge output differs"; diff /tmp/test_merge.fasta data/expected_merge.fasta | head -20; exit 1)
	@echo ""
	@echo "=== Running example_dust ==="
	./example_dust > /tmp/test_dust.fasta
	@diff /tmp/test_dust.fasta data/expected_dust.fasta \
		&& echo "PASS: dust output matches vsearch" \
		|| (echo "FAIL: dust output differs"; diff /tmp/test_dust.fasta data/expected_dust.fasta | head -20; exit 1)
	@echo ""
	@echo "=== Running example_derep ==="
	./example_derep > /tmp/test_derep.fasta
	@diff /tmp/test_derep.fasta data/expected_derep.fasta \
		&& echo "PASS: derep output matches vsearch" \
		|| (echo "FAIL: derep output differs"; diff /tmp/test_derep.fasta data/expected_derep.fasta | head -20; exit 1)
	@echo ""
	@echo "All tests passed."

clean:
	rm -f $(EXAMPLES)

.PHONY: all test clean
