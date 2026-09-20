# Makefile for VSEARCH
#
#   make                  build bin/vsearch (and the manual, if pandoc is here)
#   make DEBUG=1          debug flavour; see src/Makefile for the others
#   make install          install under $(PREFIX), honouring $(DESTDIR)
#   make dist             build the release source tarball
#
# Cross-compiling needs no separate configuration step, only a compiler:
# src/Makefile derives CC, ar and ranlib from it.
#
#   make CXX=x86_64-w64-mingw32-g++

# The version number lives in one tracked file, so that the Makefiles, the
# release-tag check in CI and anyone grepping for it all read the same thing.
# It used to be AC_INIT's second argument in configure.ac.
VERSION := $(shell cat VERSION)
ifeq ($(VERSION),)
  $(error cannot read the version number from ./VERSION)
endif

PREFIX      ?= /usr/local
exec_prefix := $(PREFIX)
datarootdir := $(PREFIX)/share
bindir      := $(exec_prefix)/bin
mandir      := $(datarootdir)/man
docdir      := $(datarootdir)/doc/vsearch
bashcompdir ?= $(datarootdir)/bash-completion/completions
zshcompdir  ?= $(datarootdir)/zsh/site-functions
fishcompdir ?= $(datarootdir)/fish/vendor_completions.d

INSTALL         ?= install
INSTALL_PROGRAM ?= $(INSTALL) -m 0755
INSTALL_DATA    ?= $(INSTALL) -m 0644
MKDIR_P         ?= $(INSTALL) -d
RM              ?= rm -f

# The manual and the completion scripts are optional, in the same two senses
# configure's --disable-manpages and --disable-completion had.
MANPAGES   ?= 1
COMPLETION ?= 1

# Where the finished binary lands, relative to this directory. A cross
# matrix that wants to keep several at once passes BINDIR=bin/<triple>; the
# path handed to src/ is absolute, so that it means the same thing there.
BINDIR ?= bin

# Ask src/ what the binary is called rather than guessing: only src/Makefile
# knows the target triple, and a stale vsearch.exe from an earlier cross build
# must not be installed alongside the native one.
BIN := $(shell $(MAKE) -s -C src print-bin)

# The '+' prefix on every use below marks the line as a recursive make.
# GNU make before 4.4 looks for the literal text "$(MAKE)" in the recipe
# *before* expanding it, so hiding it behind a variable stops the
# jobserver being passed down: "make -j24" then builds src/ with -j1 and
# says "jobserver unavailable" on its way past. The '+' says it outright.
MAKE_SRC = $(MAKE) -C src VERSION=$(VERSION) BINDIR=$(abspath $(BINDIR))

.PHONY: all vsearch lib manual install install-bin install-man \
        install-completion install-doc uninstall check clean distclean dist

all: vsearch $(if $(filter 1,$(MANPAGES)),manual)

vsearch:
	+$(MAKE_SRC)

manual:
	$(MAKE) -C man

# The library archive, for embedding vsearch in another program.
lib:
	+$(MAKE_SRC) lib

install: install-bin install-doc \
         $(if $(filter 1,$(MANPAGES)),install-man) \
         $(if $(filter 1,$(COMPLETION)),install-completion)

install-bin: vsearch
	$(MKDIR_P) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) $(BINDIR)/$(BIN) $(DESTDIR)$(bindir)

# Sectioned, so that "man -M <dir> vsearch" works straight from an unpacked
# binary tarball.  The section is the page's own suffix.
install-man:
	@for section in 1 5 7; do \
	  pages=`ls man/manpages/*.$$section 2>/dev/null` || continue; \
	  test -n "$$pages" || continue; \
	  $(MKDIR_P) $(DESTDIR)$(mandir)/man$$section; \
	  $(INSTALL_DATA) $$pages $(DESTDIR)$(mandir)/man$$section; \
	done

install-completion:
	$(MKDIR_P) $(DESTDIR)$(bashcompdir)
	$(INSTALL_DATA) completion/completions/vsearch $(DESTDIR)$(bashcompdir)/vsearch
	$(MKDIR_P) $(DESTDIR)$(zshcompdir)
	$(INSTALL_DATA) completion/completions/_vsearch $(DESTDIR)$(zshcompdir)/_vsearch
	$(MKDIR_P) $(DESTDIR)$(fishcompdir)
	$(INSTALL_DATA) completion/completions/vsearch.fish $(DESTDIR)$(fishcompdir)/vsearch.fish

install-doc:
	$(MKDIR_P) $(DESTDIR)$(docdir)
	$(INSTALL_DATA) README.md LICENSE.txt LICENSE_GNU_GPL3.txt $(DESTDIR)$(docdir)
	test -f NEWS && $(INSTALL_DATA) NEWS $(DESTDIR)$(docdir) || true

uninstall:
	$(RM) $(DESTDIR)$(bindir)/$(BIN)
	$(RM) $(DESTDIR)$(mandir)/man1/vsearch.1 $(DESTDIR)$(mandir)/man[157]/vsearch-*.[157]
	$(RM) $(DESTDIR)$(bashcompdir)/vsearch $(DESTDIR)$(zshcompdir)/_vsearch
	$(RM) $(DESTDIR)$(fishcompdir)/vsearch.fish
	$(RM) -r $(DESTDIR)$(docdir)

# Checks the completion spec against the option tables in src/cli.cc, the
# tracked scripts against what the generators emit, and each script against
# the shell it targets.
check:
	$(MAKE) -C completion check

clean:
	+$(MAKE_SRC) clean
	$(RM) -r dist

distclean: clean
	$(MAKE) -C man maintainer-clean
	$(RM) -r $(BINDIR)

# The release source tarball: the committed tree minus what .gitattributes
# marks export-ignore, plus the two generated, distributed artefacts (the
# manual pages and NEWS) so that building from the tarball needs no pandoc.
# A maintainer operation -- it needs a git checkout, exactly as 'make dist'
# always did.  'git archive HEAD' rather than 'git ls-files': it honours
# export-ignore, and a release tarball should carry what is committed, not
# whatever happens to be in the working tree.
DISTDIR := dist/vsearch-$(VERSION)

dist: manual
	@test -f man/manpages/vsearch.1 || { \
	  echo "dist: refusing to ship a distribution with no manual pages." >&2; exit 1; }
	@test -f NEWS || { echo "dist: NEWS is missing." >&2; exit 1; }
	@sources=`ls man/index.1.md man/commands/vsearch-*.md man/formats/vsearch-*.md \
	             man/misc/vsearch-*.md 2>/dev/null | wc -l`; \
	 pages=`ls man/manpages/vsearch*.[157] 2>/dev/null | wc -l`; \
	 test "$$pages" -eq "$$sources" || { \
	   echo "dist: $$pages manual pages for $$sources markdown sources." >&2; exit 1; }
	$(RM) -r $(DISTDIR)
	$(MKDIR_P) $(DISTDIR)
	git archive HEAD | tar -xf - -C $(DISTDIR)
	$(MKDIR_P) $(DISTDIR)/man/manpages
	cp man/manpages/vsearch*.[157] $(DISTDIR)/man/manpages
	cp NEWS $(DISTDIR)
	tar czf dist/vsearch-$(VERSION).tar.gz -C dist vsearch-$(VERSION)
	@echo "dist/vsearch-$(VERSION).tar.gz"
