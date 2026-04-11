#!/bin/bash

# exit if ${1} is not set

## assume it is launched from the folder containing $1

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

exit 0
