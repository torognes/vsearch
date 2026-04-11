#!/bin/bash

## assume script is launched from vsearch/man/

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

for MANPAGE in ./vsearch*.[157] ; do
    generate_pdf "${MANPAGE}" > "${MANPAGE}.pdf"
done

exit 0
