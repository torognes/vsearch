#!/bin/bash

## assume script is launched from vsearch/man/
## assume any internal link is relative to the md file itself (important)

## usage: bash scripts/generate_online_documentation.sh [CHANNEL]
##
## Render every markdown source of the manual as GitHub-flavoured
## markdown under ../docs/, the folder the GitHub Pages workflow feeds
## to Jekyll.
##
## CHANNEL is "released" (the default) or "development", and says which
## half of the published site this run builds: the released manual, from
## master, at /vsearch/, or the development manual, from dev, at
## /vsearch/dev/. GitHub Pages serves one site per repository, so the two
## are built as two folders of a single site rather than as two
## deployments, and the channel decides the baseurl the theme builds its
## asset URLs from, the site title, and the banner every page carries.
##
## The banner is not decoration. A published page names no version
## anywhere -- the front matter written below replaces the pandoc title
## block that carries it -- so without the banner a reader cannot tell
## the two manuals apart, and neither can a search engine.

## the path is computed, so only "shellcheck -x" can follow it
# shellcheck source-path=SCRIPTDIR
# shellcheck source=manpage_tools.sh
# shellcheck disable=SC1091
source "$(dirname "${0}")/manpage_tools.sh" || exit 1

## check dependencies
require_commands pandoc awk || exit 1

CHANNEL="${1:-released}"
readonly CHANNEL
case "${CHANNEL}" in
    released|development) ;;
    *) >&2 echo "Error: unknown channel '${CHANNEL}' (released or development)"
       exit 1 ;;
esac

## where each half of the site is published, spelled absolutely: the
## banner of one channel links to the other, which a baseurl-relative
## link cannot reach
readonly RELEASED_URL="https://torognes.github.io/vsearch/"
readonly DEVELOPMENT_URL="https://torognes.github.io/vsearch/dev/"

## The version the banner names is read from the title line of the hub
## page ("% vsearch(1) version 2.32.0 | vsearch manual") rather than
## passed in, so that it cannot drift from the manual it labels: both
## come from the same sources in the same checkout.
manual_version() {
    sed -n '1s/^%.* version \([^ |]*\).*/\1/p' ./index.1.md
}
VERSION="$(manual_version)"
readonly VERSION
[ -n "${VERSION}" ] || {
    >&2 echo "Error: no version in the title line of ./index.1.md"
    exit 1
}

## Every page opens by saying which of the two manuals it belongs to and
## pointing at the other. It is written into the page body rather than
## into a theme template because the site uses a remote theme: a body
## line renders wherever the theme puts the content, with nothing to
## override and nothing to keep in step with the theme's own layouts.
## called by name through write_output(), which shellcheck cannot see
# shellcheck disable=SC2317
emit_banner() {
    if [ "${CHANNEL}" = "development" ] ; then
        echo "> Development manual for **vsearch ${VERSION}**, built from the \`dev\`"
        echo "> branch: it describes changes that are not released yet. The"
        echo "> [manual for the current release](${RELEASED_URL}) is"
        echo "> published separately."
    else
        echo "> Manual for **vsearch ${VERSION}**, the current release. Changes that"
        echo "> are not released yet are described in the"
        echo "> [development manual](${DEVELOPMENT_URL})."
    fi
    echo
}

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
    ## named per page rather than as a site-wide default, which would
    ## also wrap the stylesheets the theme generates (see _config.yml)
    echo "layout: default"
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
    emit_banner
    echo "# ${1}"
    echo
    echo "${3}"
}

## The hub page is the only page whose published name differs from its
## source name: index.1.md is written as docs/index.md, so that GitHub
## Pages serves it as the site root. jekyll-relative-links rewrites a
## link only when its target exists, so a link spelled with the source
## name is left untouched and becomes a 404. Rewrite those to the
## published name. Every other page of the manual sits exactly one
## folder below the hub, which is why one '../' form is enough.
rename_hub_page_links() {
    sed 's|](\.\./index\.1\.md)|](../index.md)|g'
}

## ${1}: markdown source, ${2}: title, ${3}: parent section, ${4}: rank
## called by name through write_output(), which shellcheck cannot see
# shellcheck disable=SC2317
convert_markdown_to_github_markdown() {
    emit_front_matter "${2}" "${3}" "${4}"
    emit_banner
    expand_markdown_includes "${1}" \
        | pandoc - --to gfm+definition_lists \
        | rename_hub_page_links
}


# create folder
mkdir -p ../docs/{commands,formats,misc} || exit 1

# jekyll reads its configuration from the folder it builds, and the
# workflow builds ../docs, so the configuration has to be copied there.
# Without it the build reports "Configuration file: none", warns that the
# layout the pages ask for does not exist, and emits an unthemed site.
#
# The development manual is served one folder deeper, and just-the-docs
# builds every asset URL from baseurl, so a configuration left at
# /vsearch would send it to the released manual's stylesheets and render
# it unstyled -- the same failure the site-wide layout once caused, and
# just as silent. The title moves with it, so that the sidebar header and
# the browser tab also say which manual is open.
write_site_config() {
    if [ "${CHANNEL}" != "development" ] ; then
        cp -f ../_config.yml ../docs/_config.yml
        return
    fi
    sed -e 's|^baseurl: /vsearch$|baseurl: /vsearch/dev|' \
        -e 's|^title: vsearch manual$|title: vsearch manual (dev)|' \
        ../_config.yml > ../docs/_config.yml || return 1
    ## a substitution that quietly matched nothing would publish a
    ## development site pointing at the released one's assets
    grep -q '^baseurl: /vsearch/dev$' ../docs/_config.yml || {
        >&2 echo "Error: could not set baseurl in ../docs/_config.yml"
        return 1
    }
}
write_site_config || exit 1

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
