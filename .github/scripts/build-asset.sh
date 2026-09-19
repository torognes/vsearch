#!/bin/sh
# Build one release asset from an unpacked source distribution.
#
# Driven entirely by the environment so the workflow matrix stays data and
# this stays the only copy of the logic. Run it inside the build container
# (Linux, Windows) or directly on the runner (macOS); it assumes only a C++
# compiler, make, and the compression headers.
#
#   SRCDIR    unpacked vsearch-<version>/ directory            (required)
#   ASSET     asset base name, e.g. linux-aarch64-static       (required)
#   HOST      cross triple for --host=, empty for a native build
#   LDFLAGS_EXTRA  appended to LDFLAGS, e.g. -static-libstdc++ -static-libgcc
#   BINARY    name of the built binary            (default: vsearch)
#   FORMAT    tar | zip                           (default: tar)
#   OUTDIR    where to write the archive          (default: current directory)
#
# The asset layout deliberately matches what releases have always shipped,
# except that the monolithic man/vsearch.1 and doc/vsearch_manual.pdf are
# replaced by the modular manual, sectioned so that "man -M ./man vsearch"
# works straight from the unpacked archive.
#
# The completion/ folder is the second departure, and it closes a gap rather
# than opening one: README.md has told users to copy completion/vsearch and
# friends out of the binary distribution since ed71d158 (2026-09-03), two
# weeks before this script existed, and nothing was ever packing them.
set -eu

: "${SRCDIR:?SRCDIR is required}"
: "${ASSET:?ASSET is required}"
HOST="${HOST:-}"
LDFLAGS_EXTRA="${LDFLAGS_EXTRA:-}"
BINARY="${BINARY:-vsearch}"
FORMAT="${FORMAT:-tar}"
OUTDIR="${OUTDIR:-$(pwd)}"

version=$(sed -n 's/^AC_INIT(\[vsearch\], *\[\([^]]*\)\].*/\1/p' "${SRCDIR}/configure.ac")
test -n "${version}" || { echo "cannot read the version from configure.ac" >&2; exit 1; }
dirname="vsearch-${version}-${ASSET}"

cd "${SRCDIR}"

# --host= is what selects the architecture and OS backends in src/Makefile.am,
# so it must be passed even when the compiler is already a cross compiler.
set -- CXXFLAGS="-O3" CFLAGS="-O3"
test -z "${HOST}" || set -- "$@" --host="${HOST}"
test -z "${LDFLAGS_EXTRA}" || set -- "$@" LDFLAGS="${LDFLAGS_EXTRA}"

./configure "$@"
make -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)"

# Staging through 'make install' rather than copying by hand: it is the same
# code path distributions use, and it is what places the manual under
# share/man/man{1,5,7} for us.
rm -rf "${PWD}/.stage" "${PWD}/${dirname:?}"
make install DESTDIR="${PWD}/.stage" prefix=/

mkdir -p "${dirname}/bin"
cp ".stage/bin/${BINARY}" "${dirname}/bin/"
cp -R .stage/share/man "${dirname}/man"
cp README.md LICENSE.txt LICENSE_GNU_GPL3.txt "${dirname}/"

# The completion scripts, flattened out of the three per-shell directories
# 'make install' spreads them over: README.md tells the reader of a binary
# distribution to copy completion/vsearch, completion/_vsearch or
# completion/vsearch.fish, not to go walking share/bash-completion/. Each file
# already carries the name its shell looks up, so flattening is only a move.
# Missing means --disable-completion leaked into the build, or the install
# directories moved: either way the asset would silently lose a documented
# folder, so stop here rather than ship it.
mkdir -p "${dirname}/completion"
for script in share/bash-completion/completions/vsearch \
              share/zsh/site-functions/_vsearch \
              share/fish/vendor_completions.d/vsearch.fish ; do
  test -f ".stage/${script}" || { echo "${script} was not installed" >&2; exit 1; }
  cp ".stage/${script}" "${dirname}/completion/"
done

pages=$(find "${dirname}/man" -type f | wc -l)
scripts=$(find "${dirname}/completion" -type f | wc -l)
echo "packaged ${dirname}: $(wc -c < "${dirname}/bin/${BINARY}") byte binary, ${pages} manual pages, ${scripts} completion scripts"
test "${pages}" -gt 0 || { echo "no manual pages in the asset" >&2; exit 1; }

mkdir -p "${OUTDIR}"
if [ "${FORMAT}" = "zip" ] ; then
  zip -qr "${OUTDIR}/${dirname}.zip" "${dirname}"
else
  # --format=ustar and no xattrs: the hand-built assets carried macOS
  # AppleDouble files and com.apple.metadata xattrs, which GNU tar then
  # warns about on every extraction.
  tar --format=ustar -czf "${OUTDIR}/${dirname}.tar.gz" "${dirname}"
fi
ls -l "${OUTDIR}/${dirname}".* 
