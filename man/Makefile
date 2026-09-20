# Makefile for the vsearch manual pages
#
# The manual is a set of modular pages: a hub page (vsearch.1, generated
# from index.1.md) that lists every command, one page per command in
# section 1, one page per file format in section 5, and one page per
# reference topic in section 7.  They are generated from the markdown
# sources in this folder by scripts/generate_manpages.sh, which needs
# pandoc.
#
# Being generated, the pages are not tracked in git: a checkout needs
# pandoc to build the manual, and builds without it simply get no manual
# -- vsearch itself still builds.  'make dist' does ship them, so a build
# from a release tarball produces and installs the manual with no pandoc
# anywhere in sight.

# The page list is derived from the markdown sources rather than written
# out, so adding a command needs no edit here.  Each source already
# carries its section in the name (vsearch-cut.1.md -> vsearch-cut.1);
# only the hub page is renamed.
MD_SOURCES := $(wildcard commands/vsearch-*.md formats/vsearch-*.md misc/vsearch-*.md)
MANPAGES   := manpages/vsearch.1 \
              $(addprefix manpages/, $(basename $(notdir $(MD_SOURCES))))

FRAGMENTS := $(wildcard commands/fragments/*.md formats/fragments/*.md \
                        misc/fragments/*.md)

HAVE_PANDOC := $(shell command -v pandoc 2>/dev/null)

.PHONY: all regenerate-manpages clean
all: $(MANPAGES) $(if $(HAVE_PANDOC),../NEWS)

ifdef HAVE_PANDOC

manpages/vsearch.1: index.1.md $(FRAGMENTS)
	bash scripts/generate_manpages.sh $< $@
manpages/vsearch-%.1: commands/vsearch-%.1.md $(FRAGMENTS)
	bash scripts/generate_manpages.sh $< $@
manpages/vsearch-%.5: formats/vsearch-%.5.md $(FRAGMENTS)
	bash scripts/generate_manpages.sh $< $@
manpages/vsearch-%.7: misc/vsearch-%.7.md $(FRAGMENTS)
	bash scripts/generate_manpages.sh $< $@

# NEWS, the release history as plain text, comes from the very fragment
# vsearch-history(7) is built from, so the two cannot drift.
../NEWS: NEWS.md misc/fragments/version_history.md
	bash scripts/generate_news.sh $(abspath $@)

else

# No pandoc.  A page that is already here (the release-tarball case) counts
# as up to date and this never runs; a missing one says why rather than
# "No rule to make target".
$(MANPAGES):
	@echo "Error: the manual pages are generated from markdown by pandoc," >&2
	@echo "       which was not found.  Install pandoc, or build from a" >&2
	@echo "       release tarball, which ships them ready-made." >&2
	@exit 1

endif

# Regenerate every page, whatever was detected above.  This is how a
# maintainer refreshes the pages after editing the markdown; the script
# stops with an error if pandoc is missing.
regenerate-manpages:
	bash scripts/generate_manpages.sh

clean:
	@: # the pages are generated but also distributed: only 'make
	@: # maintainer-clean' removes what only pandoc could put back

maintainer-clean:
	$(RM) $(MANPAGES) ../NEWS
