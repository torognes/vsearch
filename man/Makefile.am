# Makefile for creating PDF manual from man file

dist_man_MANS = vsearch.1

if HAVE_PS2PDF

doc_DATA = vsearch_manual.pdf

vsearch_manual.pdf : vsearch.1
	TEMP=$$(mktemp temp.XXXXXXXX) ; \
	if [ $$(uname) == "Darwin" ] ; then \
		sed -e 's/\\-/-/g' $< | \
		iconv -f UTF-8 -t ISO-8859-1 > $$TEMP ; \
	else \
		sed -e 's/\\-/-/g' $< > $$TEMP ; \
	fi ; \
	man -t ./$$TEMP | ps2pdf -sPAPERSIZE=a4 - $@ ; \
	rm $$TEMP

CLEANFILES=vsearch_manual.pdf

endif
