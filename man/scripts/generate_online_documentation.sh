#!/bin/bash

## assume script is launched from vsearch/man/
## assume any internal link is relative to the md file itself (important)

## usage: bash scripts/generate_online_documentation.sh
##
## Render every markdown source of the manual as GitHub-flavoured
## markdown under ../docs/, the folder the GitHub Pages workflow feeds
## to Jekyll.

## the path is computed, so only "shellcheck -x" can follow it
# shellcheck source-path=SCRIPTDIR
# shellcheck source=manpage_tools.sh
# shellcheck disable=SC1091
source "$(dirname "${0}")/manpage_tools.sh" || exit 1

## check dependencies
require_commands pandoc awk || exit 1

## pandoc resolves the escaped option hyphens (\-\-cut) to literal
## double-hyphens. The literal hyphens are kept as-is by kramdown,
## the GitHub Pages renderer, thanks to the typographic_symbols
## setting in _config.yml (without it, kramdown would turn them
## into en-dashes).
## called by name through write_output(), which shellcheck cannot see
# shellcheck disable=SC2317
convert_markdown_to_github_markdown() {
    expand_markdown_includes "${1}" | pandoc - --to gfm
}


# create folder
mkdir -p ../docs/{commands,formats,misc} || exit 1

# test: maybe the config file needs to be placed at the root of the documentation?
cp -f ../_config.yml ../docs/ || exit 1

STATUS=0

# future: use vsearch.1.md as the starting page (index.html)
write_output ../docs/index.md \
             convert_markdown_to_github_markdown ./index.1.md || STATUS=1

# mirror the organization of manpages
while read -r raw_md ; do
    write_output "../docs/${raw_md#./}" \
                 convert_markdown_to_github_markdown "${raw_md}" || STATUS=1
done < <(manpage_sources | grep -v '^\./index\.1\.md$')

exit "${STATUS}"
