#!/bin/bash

## assume script is launched from vsearch/man/

## usage: bash scripts/visualize_manpage.sh commands/vsearch-cut.1.md
##
## Render one markdown source and page through it, without writing
## anything: the quickest way to see the effect of an edit.

## the path is computed, so only "shellcheck -x" can follow it
# shellcheck source-path=SCRIPTDIR
# shellcheck source=manpage_tools.sh
# shellcheck disable=SC1091
source "$(dirname "${0}")/manpage_tools.sh" || exit 1

require_commands pandoc awk man || exit 1

if [ "${#}" -ne 1 ] ; then
    >&2 echo "Usage: ${0} markdown_source"
    exit 1
fi

expand_markdown_includes "${1}" | \
    pandoc - --standalone --to man | \
    man --local-file -
