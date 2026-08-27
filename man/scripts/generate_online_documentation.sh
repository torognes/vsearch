#!/bin/bash

## assume script is launched from vsearch/man/
## assume any internal link is relative to the md file itself (important)

## check dependencies
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

convert_markdown_to_github_markdown() {
    pandoc - --to gfm
}

generate_github_markdown() {
    ## pandoc resolves the escaped option hyphens (\-\-cut) to literal
    ## double-hyphens. The literal hyphens are kept as-is by kramdown,
    ## the GitHub Pages renderer, thanks to the typographic_symbols
    ## setting in _config.yml (without it, kramdown would turn them
    ## into en-dashes).
    build_markdown_file "${1}" | \
        convert_markdown_to_github_markdown
}


# create folder
mkdir -p ../docs/{commands,formats,misc}

# test: maybe the config file needs to be placed at the root of the documentation?
ln ../_config.yml ../docs/  2> /dev/null

# future: use vsearch.1.md as the starting page (index.html)
generate_github_markdown ./index.1.md > "../docs/index.md"

# mirror the organization of manpages
for raw_md in ./{commands,formats,misc}/vsearch*.md ; do
    FOLDER="$(dirname "${raw_md}")"
    FILENAME="$(basename "${raw_md}")"
    (cd "${FOLDER}" || exit 1
     generate_github_markdown "${FILENAME}" > "../../docs/${FOLDER}/${FILENAME}"
    )
done

exit 0
