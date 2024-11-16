#!/bin/bash

## assume script is launched from vsearch/man/
## assume any internal link is relative to the md file itself (important)

build_markdown_file() {
    perl -ne \
         's/^#\((.+)\).*/`cat "$1"`/e;print' "${1}"
}

convert_markdown_to_github_markdown() {
    pandoc - --to gfm
}

generate_github_markdown() {
    build_markdown_file "${1}" | \
        convert_markdown_to_github_markdown
}


# create folder
mkdir -p ../docs/{commands,formats}

generate_github_markdown ./index.1.md > "../docs/index.md"

for raw_md in ./{commands,formats}/vsearch*.md ; do
    FOLDER="$(dirname "${raw_md}")"
    FILENAME="$(basename "${raw_md}")"
    (cd "${FOLDER}" || exit 1
     generate_github_markdown "${FILENAME}" > "../../docs/${FOLDER}/${FILENAME}"
    )
done

exit 0
