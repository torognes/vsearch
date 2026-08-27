# Makefile for the vsearch manual pages
#
# The manual is a set of modular pages: a hub page (vsearch.1, generated
# from index.1.md) that lists every command, one page per command in
# section 1, one page per file format in section 5, and one page per
# reference topic in section 7. They are generated from the markdown
# sources in this folder by scripts/generate_manpages.sh, which needs
# pandoc.
#
# The generated pages are tracked in git and shipped in release
# tarballs, exactly as configure and the Makefile.in files are, so
# building, installing and uninstalling the manual requires no pandoc
# and no perl. configure looks for pandoc all the same: when it finds
# one, the rules at the bottom of this file keep each page in step with
# the markdown source it comes from.
#
# The monolithic man/vsearch.1 page these replace is no longer built,
# installed or distributed; it is kept in the source tree, marked as
# superseded, until it is removed.

# $(wildcard) and '%'-style pattern rules below are GNU make extensions.
# vsearch already requires GNU make (3.82 or later on Linux), and the
# FreeBSD builds already use gmake, so silence automake's portability
# warnings about them rather than open-coding sixty rules.
AUTOMAKE_OPTIONS = -Wno-portability

GENERATED_MANS = \
  manpages/vsearch.1 \
  manpages/vsearch-allpairs_global.1 \
  manpages/vsearch-chimeras_denovo.1 \
  manpages/vsearch-cigar.5 \
  manpages/vsearch-cluster_fast.1 \
  manpages/vsearch-cluster_size.1 \
  manpages/vsearch-cluster_smallmem.1 \
  manpages/vsearch-cluster_unoise.1 \
  manpages/vsearch-cut.1 \
  manpages/vsearch-derep_fulllength.1 \
  manpages/vsearch-derep_id.1 \
  manpages/vsearch-derep_prefix.1 \
  manpages/vsearch-derep_smallmem.1 \
  manpages/vsearch-expected_error.7 \
  manpages/vsearch-fasta.5 \
  manpages/vsearch-fasta2fastq.1 \
  manpages/vsearch-fastq.5 \
  manpages/vsearch-fastq_chars.1 \
  manpages/vsearch-fastq_convert.1 \
  manpages/vsearch-fastq_eestats.1 \
  manpages/vsearch-fastq_eestats2.1 \
  manpages/vsearch-fastq_filter.1 \
  manpages/vsearch-fastq_join.1 \
  manpages/vsearch-fastq_mergepairs.1 \
  manpages/vsearch-fastq_stats.1 \
  manpages/vsearch-fastx_filter.1 \
  manpages/vsearch-fastx_getseq.1 \
  manpages/vsearch-fastx_getseqs.1 \
  manpages/vsearch-fastx_getsubseq.1 \
  manpages/vsearch-fastx_mask.1 \
  manpages/vsearch-fastx_revcomp.1 \
  manpages/vsearch-fastx_subsample.1 \
  manpages/vsearch-fastx_syncpairs.1 \
  manpages/vsearch-fastx_uniques.1 \
  manpages/vsearch-help.1 \
  manpages/vsearch-makeudb_usearch.1 \
  manpages/vsearch-maskfasta.1 \
  manpages/vsearch-nucleotides.7 \
  manpages/vsearch-orient.1 \
  manpages/vsearch-pairwise_alignment_parameters.7 \
  manpages/vsearch-rereplicate.1 \
  manpages/vsearch-sam.5 \
  manpages/vsearch-search_exact.1 \
  manpages/vsearch-sff.5 \
  manpages/vsearch-sff_convert.1 \
  manpages/vsearch-shuffle.1 \
  manpages/vsearch-sintax.1 \
  manpages/vsearch-sortbylength.1 \
  manpages/vsearch-sortbysize.1 \
  manpages/vsearch-uchime2_denovo.1 \
  manpages/vsearch-uchime3_denovo.1 \
  manpages/vsearch-uchime_denovo.1 \
  manpages/vsearch-uchime_ref.1 \
  manpages/vsearch-udb.5 \
  manpages/vsearch-udb2fasta.1 \
  manpages/vsearch-udbinfo.1 \
  manpages/vsearch-udbstats.1 \
  manpages/vsearch-usearch_global.1 \
  manpages/vsearch-userfields.7 \
  manpages/vsearch-version.1

# automake reads the section from the suffix, installs into man1/, man5/
# and man7/ accordingly, removes the same files on 'make uninstall', and
# ships them pregenerated in the tarball.
dist_man_MANS = $(GENERATED_MANS)

# The markdown sources are not needed to build the manual, only to
# regenerate it, but they are what a reader should patch, so they travel
# with the tarball too. Listed with $(wildcard) rather than by hand: a
# new command must not need an edit here on top of the one in
# GENERATED_MANS. 'make dist' is a maintainer operation, and vsearch
# already requires GNU make.
MANPAGE_FRAGMENTS = \
  $(wildcard $(srcdir)/commands/fragments/*.md) \
  $(wildcard $(srcdir)/formats/fragments/*.md)

EXTRA_DIST = \
  index.1.md \
  $(wildcard $(srcdir)/commands/vsearch-*.md) \
  $(wildcard $(srcdir)/formats/vsearch-*.md) \
  $(wildcard $(srcdir)/misc/vsearch-*.md) \
  $(MANPAGE_FRAGMENTS) \
  scripts/expand_includes.awk \
  scripts/generate_manpages.sh

# Regenerate every page, whatever configure decided. This is how a
# maintainer refreshes the tracked pages after editing the markdown;
# the script stops with an error if pandoc is missing.
regenerate-manpages:
	cd $(srcdir) && $(SHELL) scripts/generate_manpages.sh

.PHONY: regenerate-manpages

if HAVE_PANDOC

# configure found pandoc, so keep each page in step with its source.
# One rule per section, because the sources of the three sections live
# in three different folders.
manpages/vsearch.1: $(srcdir)/index.1.md $(MANPAGE_FRAGMENTS)
	$(SHELL) $(srcdir)/scripts/generate_manpages.sh $(srcdir)/index.1.md $@

manpages/vsearch-%.1: $(srcdir)/commands/vsearch-%.1.md $(MANPAGE_FRAGMENTS)
	$(SHELL) $(srcdir)/scripts/generate_manpages.sh $< $@

manpages/vsearch-%.5: $(srcdir)/formats/vsearch-%.5.md $(MANPAGE_FRAGMENTS)
	$(SHELL) $(srcdir)/scripts/generate_manpages.sh $< $@

manpages/vsearch-%.7: $(srcdir)/misc/vsearch-%.7.md $(MANPAGE_FRAGMENTS)
	$(SHELL) $(srcdir)/scripts/generate_manpages.sh $< $@

endif
