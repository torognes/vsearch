#!/bin/bash
#
# run_clang_tidy.sh -- quick clang-tidy check of one or more source files, to
# catch warnings a refactoring may have introduced (used like cppcheck, but with
# clang-tidy's much larger check set).
#
# It needs a compile_commands.json describing how the project is built. Generate
# it once with `bear` (regenerate after adding/removing files or changing flags):
#
#     make clean
#     bear -- make -j ARFLAGS="cr"
#
# Then check the file(s) you touched, from anywhere in the tree:
#
#     bash run_clang_tidy.sh src/core/derep.cpp
#     bash run_clang_tidy.sh src/core/cluster.cpp src/core/cluster.hpp
#
# Exit status is 1 if any warning is reported (0 if clean), so it can also gate
# a script or a pre-commit hook.
#
# Environment overrides:
#   CLANG_TIDY     clang-tidy binary          (default: clang-tidy-23, else clang-tidy)
#   CHECKS         value for --checks         (default: the list below)
#   HEADER_FILTER  value for --header-filter  (default: the checked files' own
#                                              headers; set to '.*' for all headers)
#   CXX            compiler used to locate libstdc++ (default: g++)

set -u

CLANG_TIDY="${CLANG_TIDY:-clang-tidy-23}"
command -v "${CLANG_TIDY}" >/dev/null 2>&1 || CLANG_TIDY=clang-tidy

# Keep this list in sync with the one documented in CLAUDE.md.
CHECKS="${CHECKS:--*,bugprone-*,cert-*,clang-analyzer-*,concurrency-*,cppcoreguidelines-*,hicpp-*,misc-*,modernize-*,performance-*,portability-*,readability-*}"

if [ "$#" -eq 0 ]; then
  echo "usage: $(basename "$0") <file.cpp> [more files ...]" >&2
  exit 2
fi

# Locate compile_commands.json by walking up from the current directory.
build_dir=""
dir="${PWD}"
while [ "${dir}" != "/" ]; do
  if [ -f "${dir}/compile_commands.json" ]; then build_dir="${dir}"; break; fi
  dir="$(dirname "${dir}")"
done
if [ -z "${build_dir}" ]; then
  echo "error: no compile_commands.json found at or above ${PWD}" >&2
  echo "       run 'make clean && bear -- make -j ARFLAGS=\"cr\"' first" >&2
  exit 2
fi

# Point clang at GCC's libstdc++ explicitly. Some clang-tidy/GCC combinations do
# not auto-detect it (e.g. clang-tidy 23 against GCC 13); harmless where they do.
extra_args=()
gcc_dir="$(dirname "$("${CXX:-g++}" -print-file-name=crtbegin.o 2>/dev/null)")"
case "${gcc_dir}" in
  /*) extra_args+=("--extra-arg=--gcc-install-dir=${gcc_dir}") ;;
esac

# By default, report diagnostics only in the given files and headers that share
# their basename stem (checking derep.cpp also covers derep.hpp and
# derep_internal.hpp), rather than every included project header -- this keeps
# the output focused like cppcheck. Override with HEADER_FILTER='.*' for all.
if [ -z "${HEADER_FILTER:-}" ]; then
  stems=""
  for file in "$@"; do
    stem="$(basename "${file}")"
    stem="${stem%.*}"
    stems="${stems:+${stems}|}${stem}"
  done
  HEADER_FILTER="(${stems})[^/]*\.(hpp|hh|h)$"
fi

output="$("${CLANG_TIDY}" \
  -p "${build_dir}" \
  --quiet \
  --checks="${CHECKS}" \
  --header-filter="${HEADER_FILTER}" \
  "${extra_args[@]+${extra_args[@]}}" \
  "$@" 2>&1)"
tidy_status=$?

printf '%s\n' "${output}"

# Fail if clang-tidy itself errored (e.g. a diagnostic error, unknown file) or
# if any warning was reported.
if [ "${tidy_status}" -ne 0 ] || printf '%s\n' "${output}" | grep -qE ': warning:'; then
  exit 1
fi
exit 0
