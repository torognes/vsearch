AUTOMAKE_OPTIONS = foreign
SUBDIRS = src man
EXTRA_DIST = autogen.sh

# NEWS is the release history as plain text, rendered by pandoc from the
# very fragment vsearch-history(7) is built from, so the two cannot drift.
# Like the manual pages it is generated -- hence untracked -- but shipped,
# so a build from a release tarball installs it with no pandoc in sight.
if BUILD_NEWS
dist_doc_DATA = NEWS
endif

# generated, but also distributed: 'make clean' in an unpacked tarball must
# not delete what only pandoc could put back
MAINTAINERCLEANFILES = NEWS

if HAVE_PANDOC
NEWS: $(srcdir)/man/NEWS.md $(srcdir)/man/misc/fragments/version_history.md
	cd $(srcdir)/man && bash scripts/generate_news.sh $(abs_builddir)/NEWS
else
# automake distributes a top-level NEWS whether or not we ask it to, so
# 'make dist' demands one even when configure found no pandoc. Without a
# rule here that request fails as "No rule to make target 'NEWS'", which
# says nothing about the cause. No prerequisites: an existing NEWS, the
# one a release tarball carries, counts as up to date and this never runs.
NEWS:
	@echo "Error: NEWS is generated from man/NEWS.md by pandoc, which" >&2
	@echo "       configure did not find. Install pandoc and re-run" >&2
	@echo "       configure to build a distribution." >&2
	@exit 1
endif
