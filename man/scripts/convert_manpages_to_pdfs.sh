#!/bin/bash

## assume script is launched from vsearch/man/

## usage: bash scripts/convert_manpages_to_pdfs.sh
##
## Render every manual page in ./manpages/ as a PDF next to it. This is
## a maintainer convenience, not part of the build: the PDF renderings
## are not tracked, not installed and not distributed.

## every function below is reached by name through write_output(),
## which shellcheck cannot see
# shellcheck disable=SC2317

## the path is computed, so only "shellcheck -x" can follow it
# shellcheck source-path=SCRIPTDIR
# shellcheck source=manpage_tools.sh
# shellcheck disable=SC1091
source "$(dirname "${0}")/manpage_tools.sh" || exit 1

## ps2pdf is part of ghostscript
require_commands man ps2pdf sed || exit 1

remove_dash_protections() {
    sed 's/\\-/-/g' "${1}"
}

roff_to_postscript() {
    man --troff --local-file -
}

postscript_to_pdf() {
    ps2pdf -sPAPERSIZE=a4 -
}

generate_pdf() {
    remove_dash_protections "${1}" | \
        roff_to_postscript | \
        postscript_to_pdf
}


cd ./manpages/ || exit 1

STATUS=0
FOUND=0

for MANPAGE in ./vsearch*.[157] ; do
    ## an unmatched glob comes back as the pattern itself
    [ -f "${MANPAGE}" ] || continue
    FOUND=1
    write_output "${MANPAGE}.pdf" generate_pdf "${MANPAGE}" || STATUS=1
done

if [ "${FOUND}" -eq 0 ] ; then
    >&2 echo "Error: no manual page found in ./manpages/"
    exit 1
fi

exit "${STATUS}"
