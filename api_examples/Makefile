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
CXXFLAGS = -std=c++11 -O3 -Wall -I$(VSEARCH_SRC)
LDFLAGS = $(VSEARCH_LIB) -lpthread -ldl

# Optional AddressSanitizer + UndefinedBehaviorSanitizer build:
#   make SANITIZE=1 test
# The recipes compile and link in one $(CXX) $(CXXFLAGS) ... $(LDFLAGS) step, so
# the flag on CXXFLAGS covers both. libvsearch.a must be built with the same
# sanitizers (see the library-asan CI job); mixing instrumented and
# uninstrumented objects does not link cleanly.
ifdef SANITIZE
CXXFLAGS += -fsanitize=address,undefined -fno-omit-frame-pointer -g
endif

# The same rule applies to _GLIBCXX_DEBUG, and getting it wrong is quieter than
# a link error. ../src/Makefile.am adds -D_GLIBCXX_DEBUG under --enable-debug;
# that macro replaces std::vector with std::__debug::vector, which is a
# different size, so struct Parameters has a different layout on each side of
# the boundary. An example built without it against a library built with it
# links cleanly and then reads every Parameters member past the first vector
# from the wrong offset: the examples hang, segfault, or abort with
# std::bad_alloc / std::bad_array_new_length out of a Parameters copy.
#
# So it is detected rather than assumed, from the library itself (which is the
# thing being linked, and stays right even for a hand-built libvsearch.a)
# rather than from ../src/Makefile. Debug mode instantiates std::__debug::*
# templates, whose mangled names begin _ZNSt7__debug; a release library has
# none. Override the detection with GLIBCXX_DEBUG=1 or =0 if ever needed.
GLIBCXX_DEBUG := $(shell test -r $(VSEARCH_LIB) && \
                   nm $(VSEARCH_LIB) 2>/dev/null | grep -q _ZNSt7__debug && \
                   echo 1 || echo 0)
ifeq ($(GLIBCXX_DEBUG),1)
CXXFLAGS += -D_GLIBCXX_DEBUG
endif

SHELL = /bin/bash

EXAMPLES = example_chimera example_search example_merge example_merge_error example_cluster example_dust example_derep example_reinit example_lifecycle example_dbinfo example_recover

all: $(EXAMPLES)

example_chimera: example_chimera.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

example_search: example_search.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

example_merge: example_merge.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

example_merge_error: example_merge_error.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

example_cluster: example_cluster.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

example_dust: example_dust.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

example_derep: example_derep.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

example_reinit: example_reinit.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

example_lifecycle: example_lifecycle.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

example_dbinfo: example_dbinfo.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

example_recover: example_recover.cc $(VSEARCH_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

test: $(EXAMPLES)
	@echo "=== Running example_chimera ==="
	@./example_chimera > /tmp/test_chimera.tsv || { echo "FAIL: example_chimera exited non-zero (crash or sanitizer error)"; exit 1; }
	@diff <(sort /tmp/test_chimera.tsv) <(sort data/expected_chimera.tsv) \
		&& echo "PASS: chimera output matches vsearch" \
		|| (echo "FAIL: chimera output differs"; diff /tmp/test_chimera.tsv data/expected_chimera.tsv | head -20; exit 1)
	@echo ""
	@echo "=== Running example_search ==="
	@./example_search > /tmp/test_search.tsv || { echo "FAIL: example_search exited non-zero (crash or sanitizer error)"; exit 1; }
	@diff <(sort /tmp/test_search.tsv) <(sort data/expected_search.tsv) \
		&& echo "PASS: search output matches vsearch" \
		|| (echo "FAIL: search output differs"; diff /tmp/test_search.tsv data/expected_search.tsv | head -20; exit 1)
	@echo ""
	@echo "=== Running example_cluster ==="
	@./example_cluster > /tmp/test_cluster.uc || { echo "FAIL: example_cluster exited non-zero (crash or sanitizer error)"; exit 1; }
	@diff <(grep "^[SH]" /tmp/test_cluster.uc | sort) <(grep "^[SH]" data/expected_cluster.uc | sort) \
		&& echo "PASS: cluster output matches vsearch" \
		|| (echo "FAIL: cluster output differs"; exit 1)
	@echo ""
	@echo "=== Running example_merge ==="
	@./example_merge > /tmp/test_merge.fasta || { echo "FAIL: example_merge exited non-zero (crash or sanitizer error)"; exit 1; }
	@diff /tmp/test_merge.fasta data/expected_merge.fasta \
		&& echo "PASS: merge output matches vsearch" \
		|| (echo "FAIL: merge output differs"; diff /tmp/test_merge.fasta data/expected_merge.fasta | head -20; exit 1)
	@echo ""
	@echo "=== Running example_merge_error ==="
	@./example_merge_error \
		&& echo "PASS: merge error-signal contract holds" \
		|| (echo "FAIL: merge error-signal test failed"; exit 1)
	@echo ""
	@echo "=== Running example_dust ==="
	@./example_dust > /tmp/test_dust.fasta || { echo "FAIL: example_dust exited non-zero (crash or sanitizer error)"; exit 1; }
	@diff /tmp/test_dust.fasta data/expected_dust.fasta \
		&& echo "PASS: dust output matches vsearch" \
		|| (echo "FAIL: dust output differs"; diff /tmp/test_dust.fasta data/expected_dust.fasta | head -20; exit 1)
	@echo ""
	@echo "=== Running example_derep ==="
	@./example_derep > /tmp/test_derep.fasta || { echo "FAIL: example_derep exited non-zero (crash or sanitizer error)"; exit 1; }
	@diff /tmp/test_derep.fasta data/expected_derep.fasta \
		&& echo "PASS: derep output matches vsearch" \
		|| (echo "FAIL: derep output differs"; diff /tmp/test_derep.fasta data/expected_derep.fasta | head -20; exit 1)
	@echo ""
	@echo "=== Running example_reinit ==="
	@./example_reinit \
		&& echo "PASS: re-initialization produces correct results" \
		|| (echo "FAIL: re-initialization test failed"; exit 1)
	@echo ""
	@echo "=== Running example_lifecycle ==="
	@./example_lifecycle \
		&& echo "PASS: library memory/error contracts hold" \
		|| (echo "FAIL: lifecycle contract test failed"; exit 1)
	@echo ""
	@echo "=== Running example_dbinfo ==="
	@./example_dbinfo \
		&& echo "PASS: database query/indexing surface behaves as documented" \
		|| (echo "FAIL: dbinfo test failed"; exit 1)
	@echo ""
	@echo "=== Running example_recover ==="
	@./example_recover \
		&& echo "PASS: library recovers from a fatal error and continues" \
		|| (echo "FAIL: recover-and-continue test failed"; exit 1)
	@echo ""
	@echo "All tests passed."

clean:
	rm -f $(EXAMPLES)

.PHONY: all test clean
