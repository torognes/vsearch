#!/bin/sh
# Stage the two compression DLLs the Windows asset must carry (issue #658).
#
#   BUNDLE_DIR     where to stage them                           (required)
#   MINGW_SYSROOT  the cross toolchain's sysroot  (default: /usr/x86_64-w64-mingw32)
#
# vsearch does not link against zlib or bzip2: it loads zlib1.dll and
# libbz2.dll at run time (src/os/dynlibs.cpp), and only from the directory
# vsearch.exe lives in or from System32 (LOAD_LIBRARY_SEARCH_DEFAULT_DIRS, see
# src/os/windows/dynlib_loader.cc). Neither DLL is part of Windows, so an asset
# without them reports "the library was not found" for gzip and, when bzlib.h
# was missing at build time, "Compiled without support for bzip2". That is
# what the 2.32.0 asset shipped; 2.31.0, built by hand, had both.
#
# zlib comes from Debian (libz-mingw-w64-dev, installed by the workflow).
# Debian has no mingw build of bzip2, so the header and the DLL come from the
# MSYS2 package the hand-built assets already used, pinned by checksum: a
# release must not pick up whatever the mirror serves that day. MSYS2 names
# the DLL libbz2-1.dll; vsearch asks for libbz2.dll. Both DLLs import only
# KERNEL32.dll and msvcrt.dll, so nothing else needs shipping.
#
# The bzip2 licence requires binary redistributions to reproduce its notice,
# hence LICENSE_bzip2.txt next to the other licences.
#
# This must run *before* the build: the Makefile enables bzip2 support only
# when it finds bzlib.h (the has-header probe in src/Makefile).
set -eu

: "${BUNDLE_DIR:?BUNDLE_DIR is required}"
MINGW_SYSROOT="${MINGW_SYSROOT:-/usr/x86_64-w64-mingw32}"

BZIP2_PKG=mingw-w64-x86_64-bzip2-1.0.8-4-any.pkg.tar.zst
BZIP2_URL="https://repo.msys2.org/mingw/mingw64/${BZIP2_PKG}"
BZIP2_SHA256=123768f30ae14ba654a6feb70f8526146a331bc85f831a91a549ffb3f6cbffc7

workdir=$(mktemp -d)
trap 'rm -rf "${workdir}"' EXIT

curl -sSfL --retry 3 -o "${workdir}/${BZIP2_PKG}" "${BZIP2_URL}"
echo "${BZIP2_SHA256}  ${workdir}/${BZIP2_PKG}" | sha256sum -c -
zstd -dc "${workdir}/${BZIP2_PKG}" | tar -x -C "${workdir}" \
  mingw64/include/bzlib.h \
  mingw64/bin/libbz2-1.dll \
  mingw64/share/licenses/bzip2/LICENSE

cp "${workdir}/mingw64/include/bzlib.h" "${MINGW_SYSROOT}/include/"

mkdir -p "${BUNDLE_DIR}/bin"
cp "${MINGW_SYSROOT}/lib/zlib1.dll" "${BUNDLE_DIR}/bin/zlib1.dll"
cp "${workdir}/mingw64/bin/libbz2-1.dll" "${BUNDLE_DIR}/bin/libbz2.dll"
cp "${workdir}/mingw64/share/licenses/bzip2/LICENSE" "${BUNDLE_DIR}/LICENSE_bzip2.txt"

ls -l "${BUNDLE_DIR}" "${BUNDLE_DIR}/bin"
