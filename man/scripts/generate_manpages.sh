#!/bin/bash

## assume script is launched from vsearch/man/
## assume any internal link is relative to the md file itself (important)

## with two arguments, convert a single markdown source into a single
## manual page (this is how the Makefile rebuilds one page); with no
## argument, (re)generate every page

## a failure in the awk stage must not be masked by a successful
## pandoc stage: without this, a missing markdown source still produced
## a (near-empty) manual page and an exit status of 0
set -o pipefail

## pandoc is the only tool that is not already required to build
## vsearch; the generated manual pages are tracked and shipped, so an
## ordinary build never reaches this script (see man/Makefile.am)
for dependency in pandoc awk ; do
    command -v "${dependency}" > /dev/null || \
        { >&2 echo "Error: missing ${dependency}" ; exit 1 ; }
done

## the expander runs from the folder holding the markdown source, so
## resolve its location before any change of directory
EXPAND_INCLUDES="$(cd "$(dirname "${0}")" && pwd)/expand_includes.awk"

build_markdown_file() {
    awk -f "${EXPAND_INCLUDES}" "${1}"
}

convert_markdown_to_groff() {
    pandoc - --standalone --to man
}

## ${1}: markdown source, ${2}: manual page to write. The includes are
## relative to the source, hence the change of directory; the target is
## made absolute first so that it survives it. Writing through a
## temporary file and renaming on success matters because the manual
## pages are tracked in git and installed as-is when pandoc is missing:
## a failed conversion must not leave a truncated or empty page behind.
generate_manpage() {
    local folder target
    if [ ! -r "${1}" ] ; then
        >&2 echo "Error: cannot read ${1}"
        return 1
    fi
    folder="$(dirname "${2}")"
    mkdir -p "${folder}" || return 1
    target="$(cd "${folder}" && pwd)/$(basename "${2}")"

    if (cd "$(dirname "${1}")" && build_markdown_file "$(basename "${1}")") | \
            convert_markdown_to_groff > "${target}.tmp" &&
            [ -s "${target}.tmp" ] &&
            mv -f "${target}.tmp" "${target}" ; then
        return 0
    fi
    rm -f "${target}.tmp"
    >&2 echo "Error: cannot generate ${2}"
    return 1
}


if [ "${#}" -eq 2 ] ; then
    generate_manpage "${1}" "${2}"
    exit "${?}"
fi

if [ "${#}" -ne 0 ] ; then
    >&2 echo "Usage: ${0} [markdown_source manual_page]"
    exit 1
fi


STATUS=0

## index.1.md documents the vsearch command itself and points to all
## the other pages: it is the hub page, installed as vsearch(1). It
## replaces the monolithic vsearch.1 page kept in the parent folder.
generate_manpage ./index.1.md ./manpages/vsearch.1 || STATUS=1

## one page per command (section 1), file format (section 5) and
## reference topic (section 7)
for raw_md in ./{commands,formats,misc}/vsearch*.md ; do
    generate_manpage "${raw_md}" "./manpages/$(basename "${raw_md/\.md/}")" || STATUS=1
done

exit "${STATUS}"
