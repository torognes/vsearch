#!/bin/bash

## assume script is launched from vsearch/man/
## assume any internal link is relative to the md file itself (important)

## usage:
##   bash scripts/generate_manpages.sh
##       (re)generate every manual page in ./manpages/
##   bash scripts/generate_manpages.sh <markdown_source> <manual_page>
##       convert one markdown source into one manual page; this is how
##       man/Makefile.am rebuilds a single page

## the path is computed, so only "shellcheck -x" can follow it
# shellcheck source-path=SCRIPTDIR
# shellcheck source=manpage_tools.sh
# shellcheck disable=SC1091
source "$(dirname "${0}")/manpage_tools.sh" || exit 1

## pandoc is the only tool here that is not already required to build
## vsearch; the generated manual pages are tracked and shipped, so an
## ordinary build never reaches this script (see man/Makefile.am)
require_commands pandoc awk || exit 1

## called by name through write_output(), which shellcheck cannot see
# shellcheck disable=SC2317
convert_markdown_to_groff() {
    expand_markdown_includes "${1}" | pandoc - --standalone --to man
}

## ${1}: markdown source, ${2}: manual page to write
generate_manpage() {
    write_output "${2}" convert_markdown_to_groff "${1}"
}


case "${#}" in
    0) ;;
    2) generate_manpage "${1}" "${2}" ; exit "${?}" ;;
    *) >&2 echo "Usage: ${0} [markdown_source manual_page]" ; exit 1 ;;
esac


## index.1.md documents the vsearch command itself and points to all
## the other pages: it is the hub page, generated as vsearch.1 and
## installed as vsearch(1). It replaces the monolithic vsearch.1 page
## kept in the parent folder.
STATUS=0
while read -r raw_md ; do
    generate_manpage "${raw_md}" "./manpages/$(manpage_name "${raw_md}")" || STATUS=1
done < <(manpage_sources)

exit "${STATUS}"
