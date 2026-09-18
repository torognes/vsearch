#!/bin/bash
#
# run_legacy_gcc.sh -- build vsearch with the oldest supported compilers,
# GCC 4.9 and GCC 4.8, inside a container.
#
# vsearch targets C++11 and still supports these compilers, but nothing in an
# ordinary day of work exercises them: a construct only the modern toolchain
# accepts goes in unnoticed and is found much later, by a user on an old
# distribution. This script is the local counterpart of the linux-gcc-legacy
# job in .github/workflows/build-and-test.yml, so the drift can be caught
# before pushing.
#
# Usage, from anywhere in the tree:
#
#     bash run_legacy_gcc.sh           # both compilers
#     bash run_legacy_gcc.sh 4.9       # just one
#     bash run_legacy_gcc.sh 4.8 4.9   # both, in that order
#
# The official gcc images are old Debian releases (jessie for 4.9, wheezy for
# 4.8) and already carry everything the build needs -- autoconf, automake,
# zlib.h and bzlib.h -- so nothing is installed and the container never needs
# the network. pandoc is absent, so configure reports that it is building
# without the manual pages; that is expected and not an error.
#
# The build runs on a copy of the tree made inside the container, so it leaves
# your own configure output, object files and bin/vsearch untouched. Nothing is
# written to the source directory, which is mounted read-only.
#
# Exit status is 0 when every requested compiler builds, 1 otherwise, so this
# can gate a script or a pre-commit hook.
#
# Environment overrides:
#   ENGINE        container engine             (default: podman, else docker)
#   IMAGE_PREFIX  image name, version appended (default: docker.io/library/gcc:)
#   JOBS          make -j value                (default: the container's nproc)
#   CONFIGURE_ARGS  extra ./configure arguments      (default: none)
#   STRICT        set to 1 to fail on compiler warnings as well as on errors
#                 (default: 0 -- GCC 4.8 emits one known -Wpedantic warning,
#                 for the pointer-to-function cast dlsym forces on us in
#                 src/os/posix/dynlib_loader.cc; GCC 4.9 and later do not)

set -u

DEFAULT_VERSIONS=(4.9 4.8)

ENGINE="${ENGINE:-}"
if [ -z "${ENGINE}" ]; then
  if command -v podman >/dev/null 2>&1; then
    ENGINE=podman
  elif command -v docker >/dev/null 2>&1; then
    ENGINE=docker
  else
    echo "error: neither podman nor docker was found" >&2
    echo "       install one, or set ENGINE to the container engine to use" >&2
    exit 2
  fi
fi

IMAGE_PREFIX="${IMAGE_PREFIX:-docker.io/library/gcc:}"
STRICT="${STRICT:-0}"
CONFIGURE_ARGS="${CONFIGURE_ARGS:-}"

# Locate the source root: the directory holding this script.
source_dir="$(cd "$(dirname "$0")" && pwd)"
if [ ! -f "${source_dir}/configure.ac" ]; then
  echo "error: ${source_dir} does not look like the vsearch source tree" >&2
  exit 2
fi

if [ "$#" -gt 0 ]; then
  versions=("$@")
else
  versions=("${DEFAULT_VERSIONS[@]}")
fi

# The container copies the tree to /build and works there, so the mount can
# stay read-only. The archive pipe skips the repository metadata, which is
# large, unreadable from inside an isolated build and of no use to make.
# ':z' relabels the mount for SELinux; it is ignored where SELinux is not in
# use, and docker accepts it too.
# shellcheck disable=SC2016  # JOBS and CONFIGURE_ARGS are passed in with -e
# and must expand inside the container, not here.
build_script='
set -e
mkdir /build
tar -C /src --exclude=./.git -cf - . | tar -C /build -xf -
cd /build
g++ --version | head -1
./autogen.sh
./configure CFLAGS="-O2" CXXFLAGS="-O2" ${CONFIGURE_ARGS}
make ARFLAGS="cr" -j"${JOBS:-$(nproc)}"
bin/vsearch --version
'

status=0
for version in "${versions[@]}"; do
  image="${IMAGE_PREFIX}${version}"
  echo "=============================================================="
  echo "Building vsearch with GCC ${version} (${image})"
  echo "=============================================================="

  # Stream the build and keep a copy: a legacy build takes a couple of
  # minutes, and watching it is the point of running it locally.
  log="$(mktemp)"
  "${ENGINE}" run --rm \
    -v "${source_dir}":/src:ro,z \
    -e JOBS="${JOBS:-}" \
    -e CONFIGURE_ARGS="${CONFIGURE_ARGS}" \
    "${image}" \
    bash -c "${build_script}" 2>&1 | tee "${log}"
  build_status="${PIPESTATUS[0]}"

  warnings="$(grep -c ': warning:' "${log}")"
  rm -f "${log}"
  if [ "${build_status}" -ne 0 ]; then
    echo "--> GCC ${version}: FAILED (${warnings} warning(s))"
    status=1
  elif [ "${STRICT}" = "1" ] && [ "${warnings}" -gt 0 ]; then
    echo "--> GCC ${version}: built, but ${warnings} warning(s) and STRICT=1"
    status=1
  else
    echo "--> GCC ${version}: OK (${warnings} warning(s))"
  fi
done

exit "${status}"
