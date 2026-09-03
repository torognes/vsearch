#!/bin/bash
#
# Check the shell completion scripts. Run by 'make -C completion check'.
#
#   1. drift:    does spec/vsearch.yaml still agree with src/cli.cc?
#   2. staleness: do the checked-in scripts match what the generators emit?
#   3. syntax:   does each script parse in the shell it targets?
#
# Usage: check_completions.sh [top_srcdir]
#
# A missing dependency (python3, PyYAML, zsh, fish) is reported and
# skipped, so a developer without fish installed is not blocked. Set
# STRICT=1 to turn every skip into a failure -- CI does, because a skip
# there would quietly turn the whole check into a pass.

set -u

completion_dir=$(cd -- "$(dirname -- "$0")/.." && pwd)
top_srcdir=${1:-${completion_dir}/..}
cli_cc="${top_srcdir}/src/cli.cc"

: "${PYTHON:=python3}"
: "${STRICT:=0}"

status=0
skipped=0

fail() { printf '%s\n' "FAIL: $*" >&2 ; status=1 ; }

skip() {
    if [ "${STRICT}" -ne 0 ]; then
        fail "$* (STRICT is set)"
    else
        printf '%s\n' "SKIP: $*"
        skipped=$((skipped + 1))
    fi
}

work_dir=$(mktemp -d) || { fail "cannot create a temporary directory" ; exit 1 ; }
trap 'rm -rf "${work_dir}"' EXIT

# ---------------------------------------------------------------- drift
have_python=0
if ! command -v "${PYTHON}" > /dev/null 2>&1 ; then
    skip "${PYTHON} not found: cannot check the spec against cli.cc"
elif ! "${PYTHON}" -c 'import yaml' > /dev/null 2>&1 ; then
    skip "PyYAML not installed: cannot check the spec against cli.cc"
elif [ ! -f "${cli_cc}" ]; then
    skip "${cli_cc} not found: cannot check the spec against cli.cc"
else
    have_python=1
    printf '%s\n' "== spec versus src/cli.cc =="
    if "${PYTHON}" "${completion_dir}/scripts/extract_cli_options.py" \
           "${cli_cc}" > "${work_dir}/cli.json" ; then
        "${PYTHON}" "${completion_dir}/scripts/check_option_drift.py" \
            "${work_dir}/cli.json" "${completion_dir}/spec/vsearch.yaml" \
            || fail "the completion spec has drifted from src/cli.cc"
    else
        fail "could not read the option tables out of ${cli_cc}"
    fi
    printf '\n'
fi

# ------------------------------------------------------------ staleness
if [ "${have_python}" -eq 1 ]; then
    printf '%s\n' "== checked-in scripts versus the generators =="
    for pair in "gen_bash.py:vsearch" "gen_zsh.py:_vsearch" "gen_fish.py:vsearch.fish" ; do
        generator=${pair%%:*}
        product=${pair##*:}
        if ( cd "${completion_dir}/scripts" \
                 && "${PYTHON}" "${generator}" > "${work_dir}/${product}" ) ; then
            if diff -u "${completion_dir}/completions/${product}" \
                    "${work_dir}/${product}" ; then
                printf '%s\n' "  ${product}: up to date"
            else
                fail "completions/${product} is not what ${generator} emits; run 'make -C completion regenerate-completions'"
            fi
        else
            fail "${generator} failed"
        fi
    done
    printf '\n'
fi

# --------------------------------------------------------------- syntax
printf '%s\n' "== syntax =="
for pair in "bash:vsearch" "zsh:_vsearch" "fish:vsearch.fish" ; do
    shell=${pair%%:*}
    product=${pair##*:}
    if command -v "${shell}" > /dev/null 2>&1 ; then
        if "${shell}" -n "${completion_dir}/completions/${product}" ; then
            printf '%s\n' "  ${product}: parses in ${shell}"
        else
            fail "completions/${product} does not parse in ${shell}"
        fi
    else
        skip "${shell} not installed: cannot check completions/${product}"
    fi
done

printf '\n'
if [ "${status}" -eq 0 ]; then
    if [ "${skipped}" -gt 0 ]; then
        printf '%s\n' "completion checks passed, ${skipped} skipped"
    else
        printf '%s\n' "completion checks passed"
    fi
fi
exit "${status}"
