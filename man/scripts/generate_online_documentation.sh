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
##
## GFM proper has no definition-list syntax, so a plain "--to gfm" would
## degrade each OPTIONS entry (term, then ": description") into a term, a
## hard line break and the description, all inside a single paragraph.
## The +definition_lists extension keeps them as definition lists, which
## kramdown understands and renders as <dl>/<dt>/<dd>.
## just-the-docs builds its sidebar from YAML front matter, not from the
## folder layout. A page that declares no title is listed under its first
## heading, which for a manual page is always "NAME", so every page has to
## announce the manual page it stands for and the section it belongs to.
##
## ${1}: title, ${2}: parent section (empty for the hub page), ${3}: rank
## in the sidebar (empty to let just-the-docs sort the section
## alphabetically, which is what the command pages want), ${4}: the word
## "children" when the page groups other pages under it.
emit_front_matter() {
    echo "---"
    echo "title: \"${1}\""
    [ -n "${2}" ] && echo "parent: \"${2}\""
    [ -n "${3}" ] && echo "nav_order: ${3}"
    [ "${4:-}" = "children" ] && echo "has_children: true"
    echo "---"
    echo
}

## ${1}: markdown source (for instance ./commands/vsearch-cut.1.md).
## Write the title the online manual gives it: the manual page name in
## the usual "name(section)" form, which is also how the pages refer to
## one another.
page_title() {
    local name
    name="$(manpage_name "${1}")"
    printf '%s(%s)\n' "${name%.*}" "${name##*.}"
}

## ${1}: markdown source. Write the sidebar section it belongs to, or
## nothing at all for the hub page, which sits at the top level.
page_section() {
    case "${1}" in
        ./commands/*) echo "Commands" ;;
        ./formats/*)  echo "File formats" ;;
        ./misc/*)     echo "Reference topics" ;;
        *)            echo "" ;;
    esac
}

## ${1}: title, ${2}: rank in the sidebar, ${3}: one-line summary.
##
## just-the-docs needs a page for every section its children name. These
## three carry no manual content of their own: they only group the pages
## below them.
## called by name through write_output(), which shellcheck cannot see
# shellcheck disable=SC2317
section_page() {
    emit_front_matter "${1}" "" "${2}" children
    echo "# ${1}"
    echo
    echo "${3}"
}

## ${1}: markdown source, ${2}: title, ${3}: parent section, ${4}: rank
## called by name through write_output(), which shellcheck cannot see
# shellcheck disable=SC2317
convert_markdown_to_github_markdown() {
    emit_front_matter "${2}" "${3}" "${4}"
    expand_markdown_includes "${1}" | pandoc - --to gfm+definition_lists
}


# create folder
mkdir -p ../docs/{commands,formats,misc} || exit 1

# test: maybe the config file needs to be placed at the root of the documentation?
cp -f ../_config.yml ../docs/ || exit 1

STATUS=0

# future: use vsearch.1.md as the starting page (index.html)
# the hub page opens the sidebar, above the three sections
write_output ../docs/index.md \
             convert_markdown_to_github_markdown ./index.1.md \
             "$(page_title ./index.1.md)" "" 1 || STATUS=1

# one grouping page per section, in manual-section order
write_output ../docs/commands/index.md section_page \
             "Commands" 2 \
             "One page per vsearch command (section 1 of the manual)." || STATUS=1
write_output ../docs/formats/index.md section_page \
             "File formats" 3 \
             "The file formats vsearch reads and writes (section 5)." || STATUS=1
write_output ../docs/misc/index.md section_page \
             "Reference topics" 4 \
             "Topics shared by several commands (section 7)." || STATUS=1

# mirror the organization of manpages
while read -r raw_md ; do
    write_output "../docs/${raw_md#./}" \
                 convert_markdown_to_github_markdown "${raw_md}" \
                 "$(page_title "${raw_md}")" "$(page_section "${raw_md}")" "" || STATUS=1
done < <(manpage_sources | grep -v '^\./index\.1\.md$')

exit "${STATUS}"
