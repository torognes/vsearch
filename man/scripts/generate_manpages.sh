#!/bin/bash

## assume script is launched from vsearch/man/
## assume any internal link is relative to the md file itself (important)

build_markdown_file() {
    perl -ne \
         's/^#\((.+)\).*/`cat "$1"`/e;print' "${1}"
}

convert_markdown_to_groff() {
    pandoc - --standalone --to man
}

generate_manpage() {
    build_markdown_file "${1}" | \
        convert_markdown_to_groff \
            > "../manpages/${1/\.md/}"
}


# create folder
mkdir -p manpages

for raw_md in ./{commands,formats,misc}/vsearch*.md ; do
    FOLDER="$(dirname "${raw_md}")"
    FILENAME="$(basename "${raw_md}")"
    (cd "${FOLDER}" || exit 1
     generate_manpage "${FILENAME}" > "../manpages/${FILENAME/\.md/}"
    )
done

exit 0
