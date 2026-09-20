#!/bin/sh
#
# check-object-formats.sh -- fail if a build tree mixes object formats.
#
#   sh check-object-formats.sh <build-dir> [CXX]
#
# Why this exists
# ---------------
# vsearch has two C translation units (vendored/md5.c, vendored/sha1.c) where
# swarm has none. make's built-in default for CC is the *host* "cc", so a
# cross build driven by CXX alone compiles them for the host -- and the mingw
# linker then accepts the resulting ELF objects into a PE link without a
# single diagnostic. The build exits 0, "file bin/vsearch.exe" says
# "PE32+ executable", and the binary is quietly corrupt.
#
# src/Makefile derives CC (and ar, ranlib) from CXX so this cannot happen by
# accident, but the failure is silent enough that it deserves a guard rather
# than a convention. The same check also catches a stale object left behind
# from a different target.
#
# How it checks
# -------------
# Every object in the tree must report the same "file format" as every other.
# That is target-agnostic -- no table of expected formats per triple -- and it
# is the property that actually has to hold for the link to be meaningful.
#
# Note that objdump *succeeds* on a foreign object: GNU BFD reads many formats,
# so "objdump ran without error" proves nothing. The format string is what
# matters.

# objdump's output is translated: on a French locale it prints "format de
# fichier elf64-x86-64", and an English-only pattern then matches nothing and
# reports every object as unreadable. Pin the locale before parsing anything.
LC_ALL=C
export LC_ALL

set -eu

dir=${1:-}
if [ -z "${dir}" ]; then
    printf 'usage: %s <build-dir> [CXX]\n' "$0" >&2
    exit 2
fi
if [ ! -d "${dir}" ]; then
    printf 'check-object-formats: no such directory: %s\n' "${dir}" >&2
    exit 2
fi

cxx=${2:-${CXX:-g++}}
# Ask the compiler for its own objdump, so a cross build is inspected by its
# own binutils. Falls back to a bare "objdump" for a native build.
objdump=$("${cxx}" -print-prog-name=objdump 2>/dev/null || echo objdump)
command -v "${objdump}" >/dev/null 2>&1 || objdump=objdump

objects=$(find "${dir}" -name '*.o' -print)
count=$(printf '%s\n' "${objects}" | grep -c . || true)

# A check that inspected nothing must fail, not pass: an empty build tree is
# the one way this could report success while proving nothing.
if [ "${count}" -eq 0 ]; then
    printf 'check-object-formats: no objects under %s -- nothing was checked\n' "${dir}" >&2
    exit 1
fi

report=$(
    printf '%s\n' "${objects}" | while read -r obj; do
        [ -n "${obj}" ] || continue
        fmt=$("${objdump}" -f "${obj}" 2>/dev/null |
                  sed -n 's/.*file format \(.*\)/\1/p' | head -n 1)
        printf '%s\t%s\n' "${fmt:-UNREADABLE}" "${obj}"
    done
)

formats=$(printf '%s\n' "${report}" | cut -f1 | sort -u)
distinct=$(printf '%s\n' "${formats}" | grep -c . || true)

if [ "${distinct}" -eq 1 ] && [ "${formats}" != "UNREADABLE" ]; then
    printf 'check-object-formats: OK -- %s objects, all %s\n' "${count}" "${formats}"
    exit 0
fi

printf 'check-object-formats: FAIL -- %s objects in %s do not share one format\n' \
       "${count}" "${dir}" >&2
printf '\n' >&2
printf '%s\n' "${formats}" | while read -r fmt; do
    [ -n "${fmt}" ] || continue
    n=$(printf '%s\n' "${report}" | awk -F'\t' -v f="${fmt}" '$1==f' | grep -c . || true)
    printf '  %-20s %s object(s), e.g. %s\n' "${fmt}" "${n}" \
        "$(printf '%s\n' "${report}" | awk -F'\t' -v f="${fmt}" '$1==f {print $2; exit}')" >&2
done
printf '\n' >&2
printf 'Most likely cause: CC was left as the host compiler while CXX was a\n' >&2
printf 'cross compiler. Build with CXX alone and let src/Makefile derive CC.\n' >&2
exit 1
