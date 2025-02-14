#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([vsearch], [2.29.4], [torognes@ifi.uio.no], [vsearch], [https://github.com/torognes/vsearch])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([subdir-objects])
AC_LANG([C++])
AC_CONFIG_SRCDIR([src/vsearch.cc])
AC_CONFIG_HEADERS([config.h])

AC_SUBST(MACOSX_DEPLOYMENT_TARGET)
MACOSX_DEPLOYMENT_TARGET="10.9"

# Checks for programs.
AC_PROG_CXX
AC_PROG_RANLIB
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_LIB([dl], [dlopen])
AC_CHECK_LIB([psapi], [GetProcessMemoryInfo])

# Checks for header files.
AC_CHECK_HEADERS([getopt.h fcntl.h float.h regex.h ctype.h locale.h limits.h string.h sys/time.h dlfcn.h pthread.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_CHECK_FUNCS([memmove memcpy posix_memalign gettimeofday localtime memchr memset pow regcomp strcasecmp strchr strcspn sysinfo])

have_bzip2=no
AC_ARG_ENABLE(bzip2, AS_HELP_STRING([--disable-bzip2], [Disable bzip2 support]))
AS_IF([test "x$enable_bzip2" != "xno"], [
  have_bzip2=yes
])
if test "x${have_bzip2}" = "xyes"; then
  AC_CHECK_HEADERS([bzlib.h], [], [have_bzip2=no])
fi

have_zlib=no
AC_ARG_ENABLE(zlib, AS_HELP_STRING([--disable-zlib], [Disable zlib support]))
AS_IF([test "x$enable_zlib" != "xno"], [
  have_zlib=yes
])
if test "x${have_zlib}" = "xyes"; then
  AC_CHECK_HEADERS([zlib.h], [], [have_zlib=no])
fi

have_ps2pdf=no
AC_ARG_ENABLE(pdfman, AS_HELP_STRING([--disable-pdfman], [Disable PDF manual creation]))
AS_IF([test "x$enable_pdfman" != "xno"], [
  have_ps2pdf=yes
  AC_CHECK_PROG(HAVE_PS2PDF, ps2pdf, yes, no)
  if test "x$HAVE_PS2PDF" = "xno"; then
    AC_MSG_WARN([*** ps2pdf is required to build a PDF version of the manual])
    have_ps2pdf=no
  fi
])


# Check for --enable-profiling option
AC_ARG_ENABLE([profiling],
  [AS_HELP_STRING([--enable-profiling], [Enable profiling build])],
  [enable_profiling=$enableval],
  [enable_profiling=no])

# Define AM_CONDITIONAL for profiling
AM_CONDITIONAL([ENABLE_PROFILING], [test "x$enable_profiling" = "xyes"])


# Check for --enable-debug option
AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([--enable-debug], [Enable debug build])],
  [enable_debug=$enableval],
  [enable_debug=no])

# Define AM_CONDITIONAL for debug
AM_CONDITIONAL([ENABLE_DEBUG], [test "x$enable_debug" = "xyes"])

have_man_html=no

case $target in
     aarch64*) target_aarch64="yes" ;;
     powerpc64*) target_ppc="yes" ;;
     x86_64*) target_x86_64="yes" ;;
esac

AC_CHECK_HEADERS([windows.h], [AM_CONDITIONAL(TARGET_WIN, true)], [AM_CONDITIONAL(TARGET_WIN, false)])

AM_CONDITIONAL(HAVE_PS2PDF, test "x${have_ps2pdf}" = "xyes")
AM_CONDITIONAL(HAVE_MAN_HTML, test "x${have_man_html}" = "xyes")
AM_CONDITIONAL(TARGET_PPC, test "x${target_ppc}" = "xyes")
AM_CONDITIONAL(TARGET_AARCH64, test "x${target_aarch64}" = "xyes")
AM_CONDITIONAL(TARGET_X86_64, test "x${target_x86_64}" = "xyes")
AM_PROG_CC_C_O

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 man/Makefile])

AC_OUTPUT
