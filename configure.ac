#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([vsearch], [2.7.1], [torognes@ifi.uio.no])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([subdir-objects])
AC_LANG([C++])
AC_CONFIG_SRCDIR([src/vsearch.cc])
AC_CONFIG_HEADERS([config.h])

AC_SUBST(MACOSX_DEPLOYMENT_TARGET)
MACOSX_DEPLOYMENT_TARGET="10.7"

# Set default gcc and g++ options

CFLAGS='-g'
CXXFLAGS='-g'

# Checks for programs.
AC_PROG_CXX
AC_PROG_RANLIB
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_LIB([dl], [dlopen])
#AC_CHECK_LIB([regex], [regexec])
AC_CHECK_LIB([psapi], [GetProcessMemoryInfo])

# Checks for header files.
AC_CHECK_HEADERS([getopt.h fcntl.h float.h regex.h ctype.h locale.h limits.h string.h sys/time.h dlfcn.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_FUNC_ALLOCA
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memmove memcpy posix_memalign gettimeofday localtime_r memchr memset pow regcomp strcasecmp strchr strcspn sysinfo])

have_pthreads=no
have_bzip2=no
have_zlib=no

AC_SEARCH_LIBS([pthread_create], [pthread], [have_pthreads=yes])

AC_ARG_ENABLE(bzip2, AS_HELP_STRING([--disable-bzip2], [Disable bzip2 support]))
AS_IF([test "x$enable_bzip2" != "xno"], [
  have_bzip2=yes
])

AC_ARG_ENABLE(zlib, AS_HELP_STRING([--disable-zlib], [Disable zlib support]))
AS_IF([test "x$enable_zlib" != "xno"], [
  have_zlib=yes
])

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

if test "x${have_pthreads}" = "xyes"; then
  AC_CHECK_HEADERS([pthread.h], [], [have_pthreads=no])
fi

if test "x${have_bzip2}" = "xyes"; then
  AC_CHECK_HEADERS([bzlib.h], [], [have_bzip2=no])
fi

if test "x${have_zlib}" = "xyes"; then
  AC_CHECK_HEADERS([zlib.h], [], [have_zlib=no])
fi

case $target in
     powerpc64*) target_ppc="yes" ;;
esac

AC_CHECK_HEADERS([windows.h], [AM_CONDITIONAL(TARGET_WIN, true)], [AM_CONDITIONAL(TARGET_WIN, false)])

AM_CONDITIONAL(HAVE_BZLIB, test "x${have_bzip2}" = "xyes")
AM_CONDITIONAL(HAVE_ZLIB, test "x${have_zlib}" = "xyes")
AM_CONDITIONAL(HAVE_PTHREADS, test "x${have_pthreads}" = "xyes")
AM_CONDITIONAL(HAVE_PS2PDF, test "x${have_ps2pdf}" = "xyes")
AM_CONDITIONAL(TARGET_PPC, test "x${target_ppc}" = "xyes")
AM_PROG_CC_C_O

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 man/Makefile])

AC_OUTPUT
