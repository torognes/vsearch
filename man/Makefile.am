# Makefile for creating PDF manual from man file

dist_man_MANS = vsearch.1

doc_DATA =
CLEANFILES =

if HAVE_MAN_HTML

doc_DATA += vsearch_manual.html

vsearch_manual.html : vsearch.1
	sed -e 's/\\-/-/g' $< | \
	if [ $$(uname) == "Darwin" ] ; then \
		iconv -f UTF-8 -t ISO-8859-1 ; \
	else \
		cat ; \
	fi | \
	groff -t -m mandoc -m www -Thtml > $@

CLEANFILES += vsearch_manual.html

endif


if HAVE_PS2PDF

doc_DATA += vsearch_manual.pdf

vsearch_manual.pdf : vsearch.1
	sed -e 's/\\-/-/g' $< | \
	if [ $$(uname) == "Darwin" ] ; then \
		iconv -f UTF-8 -t ISO-8859-1 ; \
	else \
		cat ; \
	fi | \
	groff -W space -t -m mandoc -T ps -P -pa4 | ps2pdf - $@

CLEANFILES += vsearch_manual.pdf

endif
