# Makefile for creating PDF manual from man file

# The modular manual pages are generated from the markdown sources in
# commands/, formats/ and misc/ by scripts/generate_manpages.sh (which
# needs pandoc). They are shipped pregenerated in release tarballs, so
# building from a tarball does not require pandoc; in a git checkout
# the stamp rule below builds them on demand.
GENERATED_MANS = \
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

dist_man_MANS = vsearch.1 $(GENERATED_MANS)

$(GENERATED_MANS): manpages.stamp
	@:
manpages.stamp:
	cd $(srcdir) && bash scripts/generate_manpages.sh
	touch $@

MAINTAINERCLEANFILES = manpages.stamp

doc_DATA =
CLEANFILES =

if HAVE_MAN_HTML

doc_DATA += vsearch_manual.html

vsearch_manual.html : vsearch.1
	sed -e 's/\\-/-/g' $< | \
	iconv -f UTF-8 -t ISO-8859-1 | \
	groff -t -m mandoc -m www -Thtml > $@

CLEANFILES += vsearch_manual.html

endif


if HAVE_PS2PDF

doc_DATA += vsearch_manual.pdf

vsearch_manual.pdf : vsearch.1
	sed -e 's/\\-/-/g' $< | \
	iconv -f UTF-8 -t ISO-8859-1 | \
	groff -W space -t -m mandoc -T ps -P -pa4 | ps2pdf - $@

CLEANFILES += vsearch_manual.pdf

endif
