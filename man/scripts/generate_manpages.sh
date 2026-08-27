#!/bin/bash

## assume script is launched from vsearch/man/
## assume any internal link is relative to the md file itself (important)

## pandoc is the only tool that is not already required to build
## vsearch; the generated manual pages are tracked and shipped, so an
## ordinary build never reaches this script (see man/Makefile.am)
for dependency in pandoc awk ; do
    command -v "${dependency}" > /dev/null || \
        { >&2 echo "Error: missing ${dependency}" ; exit 1 ; }
done

## the expander runs from inside commands/, formats/ and misc/, so
## resolve its location before changing directory
EXPAND_INCLUDES="$(cd "$(dirname "${0}")" && pwd)/expand_includes.awk"

build_markdown_file() {
    awk -f "${EXPAND_INCLUDES}" "${1}"
}

convert_markdown_to_groff() {
    pandoc - --standalone --to man
}

## ${1}: markdown source, ${2}: manual page to write. Write through a
## temporary file: the manual pages are tracked in git and installed
## as-is when pandoc is missing, so a failed conversion must not leave
## a truncated or empty page behind.
generate_manpage() {
    if build_markdown_file "${1}" | convert_markdown_to_groff > "${2}.tmp" &&
            mv -f "${2}.tmp" "${2}" ; then
        return 0
    fi
    rm -f "${2}.tmp"
    >&2 echo "Error: cannot generate ${2}"
    return 1
}


# create folder
mkdir -p manpages

STATUS=0

## index.1.md documents the vsearch command itself and points to all
## the other pages: it is the hub page, installed as vsearch(1). It
## replaces the monolithic vsearch.1 page kept in the parent folder.
generate_manpage ./index.1.md ./manpages/vsearch.1 || STATUS=1

## one page per command (section 1), file format (section 5) and
## reference topic (section 7)
for raw_md in ./{commands,formats,misc}/vsearch*.md ; do
    FOLDER="$(dirname "${raw_md}")"
    FILENAME="$(basename "${raw_md}")"
    (cd "${FOLDER}" || exit 1
     generate_manpage "${FILENAME}" "../manpages/${FILENAME/\.md/}"
    ) || STATUS=1
done

exit "${STATUS}"
