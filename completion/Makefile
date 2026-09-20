# Makefile for the vsearch shell completion scripts
#
# The scripts in completions/ are generated from spec/vsearch.yaml by the
# generators in scripts/, which need Python and PyYAML.  Unlike the manual
# pages they are tracked in git, so an ordinary build and install need
# nothing at all; only a maintainer regenerating them after a spec change
# does, through regenerate-completions below.

# Command-line assignments win over this, so 'make PYTHON=python3.11' works.
PYTHON ?= python3

.PHONY: all check regenerate-completions clean
all:

# Hooked into 'make check' at the top level: the spec against src/cli.cc, the
# tracked scripts against what the generators emit, and each script against
# the shell it targets.  A missing dependency is skipped with a notice;
# STRICT=1 makes every skip a failure, which is what CI sets.
check:
	PYTHON="$(PYTHON)" bash scripts/check_completions.sh ..

regenerate-completions:
	cd scripts && $(PYTHON) gen_bash.py > ../completions/vsearch
	cd scripts && $(PYTHON) gen_zsh.py  > ../completions/_vsearch
	cd scripts && $(PYTHON) gen_fish.py > ../completions/vsearch.fish

clean:
