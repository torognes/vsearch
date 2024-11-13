#!/usr/bin/env bash

# modified from https://stackoverflow.com/a/18517316

perl -ne 's/^#\((.+)\).*/`cat "$1"`/e;print' "${1}"
# perl -ne 's/^#\\((.+)\\).*/cat \"\\$1\"/e;print' "$@"
# perl -ne 's#^!\[\[(.+?)\]\].*#`'$0' "$1"`#e;print' "$@"

exit 0
