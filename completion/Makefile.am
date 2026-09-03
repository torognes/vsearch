# Makefile for the vsearch shell completion scripts
#
# The scripts in completions/ are generated from spec/vsearch.yaml by the
# generators in scripts/, which need Python and PyYAML. Unlike the manual
# pages, they are tracked in git: applying the manual's rule here would make
# a shipped feature depend on Python at build time, which is a heavier and
# less expected dependency for a C++ tool than pandoc is for a manual. So an
# ordinary build and install need nothing at all, only a maintainer
# regenerating them after a spec change does, through the
# regenerate-completions target below.
#
# Being tracked, they can fall behind the spec; 'make check' catches that,
# along with any drift between the spec and src/cli.cc.

# Where each shell looks for completions, substituted by configure and
# overridable with --with-bashcompdir=DIR and friends. Each file is already
# named the way it has to be installed (bash looks up the command name, zsh
# an "_"-prefixed function file, fish a "<command>.fish"), so these need no
# rename hook.
bashcompdir = @bashcompdir@
zshcompdir  = @zshcompdir@
fishcompdir = @fishcompdir@

COMPLETIONS = \
  completions/vsearch \
  completions/_vsearch \
  completions/vsearch.fish

# Installed only when configure says so -- but distributed unconditionally,
# through EXTRA_DIST below rather than a dist_ prefix here. A conditional
# dist_ primary is dropped from DISTFILES when its condition is false, which
# would make 'configure --disable-completion && make dist' produce a tarball
# silently missing the scripts.
if INSTALL_COMPLETION
bashcomp_DATA = completions/vsearch
zshcomp_DATA  = completions/_vsearch
fishcomp_DATA = completions/vsearch.fish
endif

# Command-line assignments win over this, so 'make PYTHON=python3.11' works.
PYTHON = python3

EXTRA_DIST = \
  $(COMPLETIONS) \
  README.md \
  spec/vsearch.yaml \
  scripts/spec.py \
  scripts/gen_bash.py \
  scripts/gen_zsh.py \
  scripts/gen_fish.py \
  scripts/extract_cli_options.py \
  scripts/check_option_drift.py \
  scripts/check_completions.sh

# Regenerate the three scripts from the spec, whatever configure decided.
# This is how a maintainer refreshes the tracked scripts after editing
# spec/vsearch.yaml; the generators stop with an error if PyYAML is missing.
# They write to stdout, and into $(srcdir) rather than the build directory,
# because their output is what gets committed.
regenerate-completions:
	cd $(srcdir)/scripts && $(PYTHON) gen_bash.py > ../completions/vsearch
	cd $(srcdir)/scripts && $(PYTHON) gen_zsh.py > ../completions/_vsearch
	cd $(srcdir)/scripts && $(PYTHON) gen_fish.py > ../completions/vsearch.fish

# Hooked into 'make check': the spec against src/cli.cc, the tracked scripts
# against what the generators emit, and each script against the shell it
# targets. A missing dependency is skipped with a notice; STRICT=1 makes
# every skip a failure, which is what CI sets.
check-local:
	PYTHON="$(PYTHON)" bash $(srcdir)/scripts/check_completions.sh $(top_srcdir)

.PHONY: regenerate-completions
