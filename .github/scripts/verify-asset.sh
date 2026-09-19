#!/bin/sh
# Check one packaged release asset before it is published.
#
#   usage: verify-asset.sh <assets-dir> <version>
#   ASSET      asset base name, e.g. linux-aarch64-static   (required)
#   LD         the LDFLAGS the asset was linked with, so that
#              -static-libstdc++ can be checked to have taken effect
#   MAX_GLIBC  highest glibc symbol version the binary may need (default 2.31,
#              the debian:11 floor)
#
# readelf rather than objdump or ldd: it parses any ELF regardless of the host
# architecture, which is what makes the cross-built assets checkable here.
set -eu

assets_dir="${1:?usage: verify-asset.sh <assets-dir> <version>}"
version="${2:?usage: verify-asset.sh <assets-dir> <version>}"
: "${ASSET:?ASSET is required}"
LD="${LD:-}"
MAX_GLIBC="${MAX_GLIBC:-2.31}"

dirname="vsearch-${version}-${ASSET}"
workdir=$(mktemp -d)
trap 'rm -rf "${workdir}"' EXIT

archive="${assets_dir}/${dirname}.tar.gz"
if [ -f "${assets_dir}/${dirname}.zip" ] ; then
  archive="${assets_dir}/${dirname}.zip"
  unzip -q "${archive}" -d "${workdir}"
else
  tar xzf "${archive}" -C "${workdir}"
fi
root="${workdir}/${dirname}"

status=0
fail() { echo "::error::${ASSET}: $*" >&2 ; status=1 ; }

# --- contents -------------------------------------------------------------
pages=$(find "${root}/man" -type f 2>/dev/null | wc -l)
test "${pages}" -gt 0 || fail "no manual pages"
for f in README.md LICENSE.txt LICENSE_GNU_GPL3.txt ; do
  test -f "${root}/${f}" || fail "${f} is missing"
done
# Checked by name, not by count: each script is installed under the name its
# shell looks up, so a rename would leave the folder looking complete while
# completion silently never triggers.
for f in completion/vsearch completion/_vsearch completion/vsearch.fish ; do
  test -s "${root}/${f}" || fail "${f} is missing or empty"
done
# the hand-built assets carried these; a Linux CI job must not reintroduce them
test "$(find "${root}" -name '._*' | wc -l)" -eq 0 || fail "AppleDouble files in the archive"

binary="${root}/bin/vsearch"
test -f "${binary}" || binary="${root}/bin/vsearch.exe"
test -f "${binary}" || { fail "no binary in bin/" ; exit 1 ; }

echo "${ASSET}: ${pages} manual pages, $(find "${root}/completion" -type f 2>/dev/null | wc -l) completion scripts, $(wc -c < "${binary}") byte binary"

# --- ELF checks (skipped for the Windows PE asset) ------------------------
case "${ASSET}" in
  win-*)
    echo "  PE binary, skipping the ELF checks"
    ;;
  *)
    needed=$(readelf -d "${binary}" 2>/dev/null | sed -n 's/.*(NEEDED).*\[\(.*\)\]/\1/p' | tr '\n' ' ')
    echo "  NEEDED: ${needed:-<none>}"

    floor=$(readelf -V "${binary}" 2>/dev/null | grep -o 'GLIBC_[0-9.]*' | sed 's/GLIBC_//' | sort -V | tail -1)
    echo "  needs glibc >= ${floor:-<none>}"
    if [ -n "${floor}" ] ; then
      highest=$(printf '%s\n%s\n' "${floor}" "${MAX_GLIBC}" | sort -V | tail -1)
      test "${highest}" = "${MAX_GLIBC}" || \
        fail "needs glibc ${floor}, above the ${MAX_GLIBC} floor -- was it built in the pinned container?"
    fi

    # The whole point of option (c): the -static assets must carry no
    # libstdc++ dependency, and the others must still be ordinary dynamic
    # builds, so a silent flag regression in either direction is caught.
    case "${LD}" in
      *-static-libstdc++*)
        case "${needed}" in
          *libstdc++*) fail "-static-libstdc++ did not take effect" ;;
          *) echo "  libstdc++: statically linked, as intended" ;;
        esac
        ;;
      *)
        case "${needed}" in
          *libstdc++*) echo "  libstdc++: dynamic, as intended" ;;
          *) fail "expected a dynamic libstdc++ dependency and found none" ;;
        esac
        ;;
    esac
    ;;
esac

exit "${status}"
