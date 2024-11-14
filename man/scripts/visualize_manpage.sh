#!/bin/bash

bash ./include_markdown_files.sh "${1}" | \
    pandoc - --standalone --to man | \
    /usr/bin/man -l -

exit 0
