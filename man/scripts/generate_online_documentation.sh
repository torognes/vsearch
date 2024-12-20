#!/bin/bash

## assume script is launched from vsearch/man/
## assume any internal link is relative to the md file itself (important)

## check dependencies
for dependency in pandoc perl ; do
    which "${dependency}" > /dev/null || \
        { >&2 echo "Error: missing ${dependency}" ; exit 1 ; }
done

build_markdown_file() {
    perl -ne \
         's/^#\((.+)\).*/`cat "$1"`/e;print' "${1}"
}

convert_markdown_to_github_markdown() {
    pandoc - --to gfm
}

generate_github_markdown() {
    ## Failed tests:
    # sed 's/\\\-\\\-/\-\-/g'
    # sed 's/\\\-\\\-/\\-\\-/g'
    # sed 's/\\\-\\\-/\\\-\\\-/g'
    # sed 's/\\\-\\\-/\\\\-\\\\-/g'
    # sed 's/\\\-\\\-/\\\\\-\\\\\-/g'
    build_markdown_file "${1}" | \
        sed 's/\\\-\\\-/\\\\\\-\\\\\\-/g' | \
        convert_markdown_to_github_markdown
}


# create folder
mkdir -p ../docs/{commands,formats,misc}

# test: maybe the config file needs to be placed at the root of the documentation? does not work
# retry with:
# ln ../_config.yml ../docs/  2> /dev/null

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
