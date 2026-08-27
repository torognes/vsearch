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

# Without pandoc the conditionals above quietly drop the manual and NEWS
# from DISTFILES, so 'make dist' would produce a tarball missing both and
# say nothing. Count what actually reached the distribution: one manual
# page per markdown source, plus NEWS. This also catches a page whose
# markdown source exists but which was never added to GENERATED_MANS.
dist-hook:
	@sources=`ls $(srcdir)/man/index.1.md \
	             $(srcdir)/man/commands/vsearch-*.md \
	             $(srcdir)/man/formats/vsearch-*.md \
	             $(srcdir)/man/misc/vsearch-*.md 2>/dev/null | wc -l`; \
	pages=`ls $(distdir)/man/manpages/vsearch*.1 \
	          $(distdir)/man/manpages/vsearch*.5 \
	          $(distdir)/man/manpages/vsearch*.7 2>/dev/null | wc -l`; \
	status=0; \
	if test "$$pages" -ne "$$sources"; then \
	  echo "dist-hook: $$pages manual pages for $$sources markdown sources." >&2; \
	  status=1; \
	fi; \
	if test ! -f "$(distdir)/NEWS"; then \
	  echo "dist-hook: NEWS is missing from the distribution." >&2; \
	  status=1; \
	fi; \
	test "$$status" -eq 0 || { \
	  echo "dist-hook: refusing to ship an incomplete distribution; pandoc" >&2; \
	  echo "           generates the manual and NEWS -- install it, re-run" >&2; \
	  echo "           configure, and try again." >&2; \
	  exit 1; \
	}
