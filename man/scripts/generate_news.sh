#!/bin/bash

## assume script is launched from vsearch/man/

## usage: bash scripts/generate_news.sh [news_file]
##
## Render NEWS.md -- a short preamble plus the release-history fragment
## that vsearch-history(7) is built from -- as the plain-text NEWS file
## the GNU coding standards ask for. The two therefore cannot drift.
## Defaults to ../NEWS; the Makefile passes an absolute path so that a
## VPATH build writes into the build tree and not into the sources.

## the path is computed, so only "shellcheck -x" can follow it
# shellcheck source-path=SCRIPTDIR
# shellcheck source=manpage_tools.sh
# shellcheck disable=SC1091
source "$(dirname "${0}")/manpage_tools.sh" || exit 1

require_commands pandoc awk || exit 1

if [ "${#}" -gt 1 ] ; then
    >&2 echo "Usage: ${0} [news_file]"
    exit 1
fi

## 'markdown-smart' keeps the output pure ASCII: without it pandoc turns
## the quotes and dashes of the release notes into their typographic
## variants, which is not what a plain-text NEWS file wants
## called by name through write_output(), which shellcheck cannot see
# shellcheck disable=SC2317
convert_markdown_to_plain_text() {
    expand_markdown_includes "${1}" | \
        pandoc - --from markdown-smart --to plain --columns=79
}

write_output "${1:-../NEWS}" convert_markdown_to_plain_text ./NEWS.md
