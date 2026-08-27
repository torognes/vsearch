#!/bin/bash

## Helpers shared by the scripts in this folder. Source it, do not run
## it:
##
##     source "$(dirname "${0}")/manpage_tools.sh"
##
## All of these scripts are launched from vsearch/man/ and convert the
## markdown sources of the manual into something else (manual pages,
## GitHub-flavoured markdown, PDF renderings). What they have in common
## is gathered here so that they check their dependencies, resolve the
## '#(path)' includes and write their output the same way.

## a failure in an early stage of a conversion pipeline must not be
## masked by a later stage succeeding: without this, feeding a missing
## markdown source to the pipeline still produced an output file and an
## exit status of 0
set -o pipefail

## resolved once, so that the helpers keep working after the change of
## directory that resolving relative includes requires
MANPAGE_TOOLS_FOLDER="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly MANPAGE_TOOLS_FOLDER

## report every missing command, not just the first one
require_commands() {
    local dependency
    local status=0
    for dependency in "${@}" ; do
        command -v "${dependency}" > /dev/null || {
            >&2 echo "Error: missing ${dependency}"
            status=1
        }
    done
    return "${status}"
}

## ${1}: markdown source. Expand its '#(path)' include directives on
## the standard output. The includes are relative to the file itself,
## hence the change of directory.
expand_markdown_includes() {
    [ -r "${1}" ] || { >&2 echo "Error: cannot read ${1}" ; return 1 ; }
    ( cd "$(dirname "${1}")" &&
          awk -f "${MANPAGE_TOOLS_FOLDER}/expand_includes.awk" "$(basename "${1}")" )
}

## ${1}: file to write, ${2}...: a command (usually one of the
## conversion functions defined by the caller) whose standard output
## becomes that file.
##
## The output goes to a temporary file and is renamed only once the
## command succeeded and produced something. A conversion that fails
## halfway -- a missing tool, a missing source -- therefore leaves the
## previous version of the file untouched, instead of replacing it with
## a truncated or empty one. That matters most for the manual pages:
## they are tracked in git and installed as-is when pandoc is missing.
write_output() {
    local target="${1}"
    local folder
    shift
    folder="$(dirname "${target}")"
    mkdir -p "${folder}" || return 1
    if "${@}" > "${target}.tmp" && [ -s "${target}.tmp" ] &&
            mv -f "${target}.tmp" "${target}" ; then
        return 0
    fi
    rm -f "${target}.tmp"
    >&2 echo "Error: cannot write ${target}"
    return 1
}

## ${1}: markdown source (for instance ./commands/vsearch-cut.1.md).
## Write the name of the manual page it yields (vsearch-cut.1). The hub
## page is the exception: index.1.md becomes vsearch.1, since it is the
## page documenting the vsearch command itself.
manpage_name() {
    local name
    name="$(basename "${1}" .md)"
    if [ "${name}" = "index.1" ] ; then
        echo "vsearch.1"
    else
        echo "${name}"
    fi
}

## every markdown source of the manual, hub page first, then the
## commands (section 1), the file formats (section 5) and the reference
## topics (section 7)
manpage_sources() {
    echo ./index.1.md
    ls ./commands/vsearch-*.md ./formats/vsearch-*.md ./misc/vsearch-*.md
}
