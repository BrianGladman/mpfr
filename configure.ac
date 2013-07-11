dnl Process this file with autoconf to produce a configure script.

AC_COPYRIGHT([
Copyright 1999-2013 Free Software Foundation, Inc.
Contributed by the AriC and Caramel projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation; either version 3 of the License, or (at
your option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
])

dnl Add check-news when it checks for more than 15 lines
AC_INIT([MPFR],[3.2.0-dev])

dnl Automake 1.11.2 is required for AM_PROG_AR (but note that
dnl unless a security patch has been applied to Automake, e.g.
dnl via an OS distribution, Automake versions up to 1.11.6 are
dnl vulnerable for "make distcheck").
AM_INIT_AUTOMAKE([1.11.2 no-define dist-bzip2 dist-xz dist-zip])
AM_MAINTAINER_MODE(enable)

AC_CONFIG_MACRO_DIR([m4])

dnl FIXME: The AC_ARG_ENABLE(decimal-float...) part does things too
dnl early, even when this option is not used. In particular, it must
dnl be put after AC_PROG_CC; another problem is that the GMP CFLAGS
dnl and CC check may modify the compiler.

test_CFLAGS=${CFLAGS+set}

dnl Check if user request its CC and CFLAGS
if test -n "$CFLAGS" || test -n "$CC" ; then
 user_redefine_cc=yes
fi

dnl Basic Autoconf macros. At this point, they must not make Autoconf
dnl choose a compiler because of the CC and CFLAGS setup from gmp.h!

AC_CANONICAL_HOST

AC_PROG_EGREP
AC_PROG_SED

dnl To use a separate config header.
dnl There is still some problem with GMP's HAVE_CONFIG
dnl AC_CONFIG_HEADERS([mpfrconf.h:mpfrconf.in])

dnl Extra arguments to configure
unset gmp_lib_path GMP_CFLAGS GMP_CC
AC_ARG_WITH(gmp_include,
   [  --with-gmp-include=DIR  GMP include directory ],
   CPPFLAGS="$CPPFLAGS -I$withval")
AC_ARG_WITH(gmp_lib,
   [  --with-gmp-lib=DIR      GMP lib directory ], [
   LDFLAGS="$LDFLAGS -L$withval"
   gmp_lib_path="$withval"
  ])
AC_ARG_WITH(gmp,
   [  --with-gmp=DIR          GMP install directory ], [
   if test -z "$with_gmp_lib" && test -z "$with_gmp_include" ; then
      CPPFLAGS="$CPPFLAGS -I$withval/include"
      LDFLAGS="$LDFLAGS -L$withval/lib"
      gmp_lib_path="$withval/lib"
   else
      AC_MSG_FAILURE([Do not use --with-gmp and --with-gmp-include/--with-gmp-lib options simultaneously.])
   fi
  ])

AC_ARG_WITH(gmp_build,
   [  --with-gmp-build=DIR    GMP build directory (please read INSTALL file)],
   [
   if test -z "$gmp_lib_path" && test -z "$with_gmp_include" ; then
      CPPFLAGS="$CPPFLAGS -I$withval -I$withval/tune"
      LDFLAGS="$LDFLAGS -L$withval -L$withval/.libs -L$withval/tune"
      gmp_lib_path="$withval$PATH_SEPARATOR$withval/.libs$PATH_SEPARATOR$withval/tune"
      if test -r $withval/Makefile ; then
         GMP_CFLAGS=`$SED -n 's/^CFLAGS = //p' $withval/Makefile`
         GMP_CC=`$SED -n 's/^CC = //p' $withval/Makefile`
         GMP_SOURCE=`$SED -n 's/^srcdir = *//p' $withval/Makefile`
         case "$GMP_SOURCE" in
           .)  GMP_SOURCE="" ;;
           /*) ;;
           ?*) GMP_SOURCE="$withval/$GMP_SOURCE" ;;
         esac
         if test -d "$GMP_SOURCE" ; then
            CPPFLAGS="$CPPFLAGS -I$GMP_SOURCE -I$GMP_SOURCE/tune"
         fi
      fi
      use_gmp_build=yes
   else
      AC_MSG_FAILURE([Do not use --with-gmp-build and other --with-gmp options simultaneously.])
   fi
   ])

AC_ARG_WITH(mulhigh_size,
   [  --with-mulhigh-size=NUM internal threshold table for mulhigh],
   AC_DEFINE_UNQUOTED([MPFR_MULHIGH_SIZE],$withval, [Mulhigh size]))

AC_ARG_ENABLE(gmp-internals,
   [  --enable-gmp-internals  enable use of GMP undocumented functions [[default=no]]],
   [ case $enableval in
      yes) AC_DEFINE([WANT_GMP_INTERNALS],1,[Want GMP internals]) ;;
      no)  ;;
      *) AC_MSG_ERROR([bad value for --enable-gmp-internals: yes or no]) ;;
     esac])

AC_ARG_ENABLE(assert,
   [  --enable-assert         enable ASSERT checking [[default=no]]],
   [ case $enableval in
      yes) AC_DEFINE([MPFR_WANT_ASSERT],1,[Want assertion]) ;;
      none) AC_DEFINE([MPFR_WANT_ASSERT],-1,[Do no want any assertion]) ;;
      no)  ;;
      full) AC_DEFINE([MPFR_WANT_ASSERT],2,[Want assertion]) ;;
      *) AC_MSG_ERROR([bad value for --enable-assert: yes, no, none or full]) ;;
     esac])
AC_ARG_ENABLE(logging,
   [  --enable-logging        enable MPFR logging (the system must support it)
                          [[default=no]]],
   [ case $enableval in
      yes) AC_DEFINE([MPFR_USE_LOGGING],1,[Log what MPFR does]) ;;
      no)  ;;
      *)   AC_MSG_ERROR([bad value for --enable-logging: yes or no]) ;;
     esac])
AC_ARG_ENABLE(thread-safe,
   [  --disable-thread-safe   explicitly disable TLS support
  --enable-thread-safe    build MPFR as thread safe, i.e. with TLS support
                          (the system must support it) [[default=autodetect]]],
   [ case $enableval in
      yes) ;;
      no)  ;;
      *)   AC_MSG_ERROR([bad value for --enable-thread-safe: yes or no]) ;;
     esac])
AC_ARG_ENABLE(warnings,
   [  --enable-warnings       allow MPFR to output warnings to stderr [[default=no]]],
   [ case $enableval in
      yes) AC_DEFINE([MPFR_USE_WARNINGS],1,[Allow MPFR to output warnings to stderr]) ;;
      no)  ;;
      *)   AC_MSG_ERROR([bad value for --enable-warnings: yes or no]) ;;
     esac])

AC_ARG_ENABLE(tests-timeout,
   [  --enable-tests-timeout=NUM    enable timeout (NUM seconds) for test programs
                          (NUM <= 9999) [[default=no]]; if enabled, env variable
                          $MPFR_TESTS_TIMEOUT overrides NUM (0: no timeout).],
   [ case $enableval in
      no)   ;;
      yes)  AC_DEFINE([MPFR_TESTS_TIMEOUT], 0, [timeout limit]) ;;
      [[0-9]]|[[0-9]][[0-9]]|[[0-9]][[0-9]][[0-9]]|[[0-9]][[0-9]][[0-9]][[0-9]])
       AC_DEFINE_UNQUOTED([MPFR_TESTS_TIMEOUT], $enableval, [timeout limit]) ;;
      *)    AC_MSG_ERROR([bad value for --enable-tests-timeout]) ;;
     esac])

AC_ARG_ENABLE(tune-for-coverage,
   [  --tune-for-coverage     tune MPFR for coverage tests.],
   [ case $enableval in
      no)   ;;
      yes)  AC_DEFINE([MPFR_TUNE_COVERAGE], 1, [tune for coverage]) ;;
      *)    AC_MSG_ERROR([bad value for --enable-tune-for-coverage]) ;;
     esac])


dnl
dnl Setup CC and CFLAGS
dnl

dnl ********************************************************************
dnl Check for CC and CFLAGS in gmp.h

dnl Warning! The following tests must be done before Autoconf selects
dnl a compiler. This means that some macros such as AC_PROG_CC and
dnl AM_PROG_AR must be put after the following code.

if test -z "$user_redefine_cc" && test "$cross_compiling" != yes ; then

dnl We need to guess the C preprocessor instead of using AC_PROG_CPP,
dnl since AC_PROG_CPP implies AC_PROG_CC, which chooses a compiler
dnl (before we have the chance to get it from gmp.h) and does some
dnl checking related to this compiler (such as dependency tracking
dnl options); if the compiler changes due to __GMP_CC in gmp.h, one
dnl would have incorrect settings.
dnl FIXME: Move this in acinclude.m4 ?
if test -z "$GMP_CC$GMP_CFLAGS" ; then
   AC_MSG_CHECKING(for CC and CFLAGS in gmp.h)
   GMP_CC=__GMP_CC
   GMP_CFLAGS=__GMP_CFLAGS
   for cpp in /lib/cpp gcc cc c99
   do
     test $cpp = /lib/cpp || cpp="$cpp -E"
     echo foo > conftest.c
     if $cpp $CPPFLAGS conftest.c > /dev/null 2> /dev/null ; then
       # Get CC
       echo "#include \"gmp.h\"" >  conftest.c
       echo "MPFR_OPTION __GMP_CC"           >> conftest.c
       GMP_CC=`$cpp $CPPFLAGS conftest.c 2> /dev/null | $EGREP MPFR_OPTION | $SED -e 's/MPFR_OPTION //g;s/ *" *//g'`
       # Get CFLAGS
       echo "#include \"gmp.h\"" >  conftest.c
       echo "MPFR_OPTION __GMP_CFLAGS"           >> conftest.c
       unset rmpedantic
       [test "$enable_logging" = yes && rmpedantic='s/[ "]-pedantic[ "]/ /g;']
       GMP_CFLAGS=`$cpp $CPPFLAGS conftest.c 2> /dev/null | $EGREP MPFR_OPTION | $SED -e "$rmpedantic"'s/MPFR_OPTION //g;s/ *" *//g'`
       break
     fi
   done
   rm -f conftest*
   if test "x$GMP_CC" = "x__GMP_CC" || test "x$GMP_CFLAGS" = "x__GMP_CFLAGS" ; then
      AC_MSG_RESULT(no)
      GMP_CFLAGS=
      GMP_CC=
   else
      AC_MSG_RESULT(yes [CC=$GMP_CC CFLAGS=$GMP_CFLAGS])
   fi
fi

dnl But these variables may be invalid, so we must check them first.
dnl Note: we do not use AC_RUN_IFELSE, as it implies AC_PROG_CC.
if test -n "$GMP_CC$GMP_CFLAGS" ; then
   AC_MSG_CHECKING(for CC=$GMP_CC and CFLAGS=$GMP_CFLAGS)
   echo "int main (void) { return 0; }" > conftest.c
   if $GMP_CC $GMP_CFLAGS -o conftest conftest.c 2> /dev/null ; then
     AC_MSG_RESULT(yes)
     CFLAGS=$GMP_CFLAGS
     CC=$GMP_CC
   else
     AC_MSG_RESULT(no)
   fi
   rm -f conftest*
fi

fi

dnl ********************************************************************

AC_PROG_CC
AC_PROG_CPP
AC_LANG(C)

dnl Support unusual archivers.
AM_PROG_AR

dnl For GCC, _Decimal64 was introduced in GCC 4.3 for some targets
dnl (note that it is not guaranteed to be available because it may
dnl be disabled in the GCC build). See:
dnl   http://gcc.gnu.org/gcc-4.3/changes.html
dnl _Decimal64 is not yet defined in GCC for C++:
dnl   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=51364
AC_ARG_ENABLE(decimal-float,
   [  --enable-decimal-float  build conversion functions from/to decimal floats
                          [[default=no]]],
   [ case $enableval in
      yes) AC_DEFINE([MPFR_WANT_DECIMAL_FLOATS],1,
              [Build decimal float functions])
           AC_MSG_CHECKING(if compiler knows _Decimal64)
           AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[_Decimal64 x;]])],
              [AC_MSG_RESULT(yes)],
              [AC_MSG_RESULT(no)
AC_MSG_ERROR([Compiler doesn't know _Decimal64 (ISO/IEC TR 24732).
Please use another compiler or build MPFR without --enable-decimal-float.])]
           )
           AC_MSG_CHECKING(decimal float format)
           AC_RUN_IFELSE([AC_LANG_PROGRAM([[
#include <stdlib.h>
]], [[
union { double d; _Decimal64 d64; } y;
y.d64 = 1234567890123456.0dd;
return y.d == 0.14894469406741037E-123 ? 0 :
       y.d == 0.59075095508629822E-68  ? 1 : 2;
]])], [AC_MSG_RESULT(DPD)
       AC_DEFINE([DPD_FORMAT],1,[])],
      [if test "$?" != 1 ; then
          AC_MSG_FAILURE(neither DPD nor BID)
       fi
       AC_MSG_RESULT(BID)],
      [AC_MSG_RESULT(assuming DPD)
       AC_DEFINE([DPD_FORMAT],1,[])])
           ;;
      no)  ;;
      *)   AC_MSG_ERROR([bad value for --enable-decimal-float: yes or no]) ;;
     esac])

AC_ARG_ENABLE(float128,
   [  --enable-float128       build conversion functions from/to __float128
                          [[default=no]]],
   [ case $enableval in
      yes) AC_DEFINE([MPFR_WANT_FLOAT128],1,
              [Build float128 conversion functions])
           AC_MSG_CHECKING(if compiler knows __float128)
           AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[__float128 x;]])],
              [AC_MSG_RESULT(yes)],
              [AC_MSG_ERROR([Compiler doesn't know __float128; try GCC >= 4.3]
           )]);;
      no)  ;;
      *)   AC_MSG_ERROR([bad value for --enable-float128: yes or no]) ;;
     esac])

dnl Check if compiler is ICC, and if such a case, disable GCC
dnl And add some specific flags.
dnl Don't add Warnings Flags (Otherwise you'll get more than 20000 warnings).
dnl Add -long_double flags? Don't use -pc64 !
dnl Notes (VL):
dnl   * With icc 10.1 20080212 on itanium, the __ICC macro is not defined,
dnl     even when the -icc option is used (contrary to what is documented
dnl     on the icc man page).
dnl   * When ICC is correctly detected (__ICC macro defined), unsetting
dnl     the GCC variable confuses libtool. See:
dnl       http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=485421
dnl   * If need be, the gcc predefined macros __GNUC_* can be disabled
dnl     thanks to the -no-gcc option.
dnl   * Now use -mieee-fp instead of -mp (ICC 13 says: option '-mp' is
dnl     deprecated and will be removed in a future release.). According
dnl     to "icc -help", both options are equivalent; it also suggests
dnl     to use -fp-model <arg> instead of -mp, but what is <arg>? Anyway
dnl     -mieee-fp solves the tset_d failure due to a negative zero. This
dnl     option has been introduced in ICC 9.0 for Linux (2005-06):
dnl     http://www.ncsa.illinois.edu/UserInfo/Resources/Software/Intel/Compilers/9.0/C_ReleaseNotes.htm
AC_MSG_CHECKING(for ICC)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#if !defined(__ICC)
# error "ICC Not Found"
error
#endif
]], [[]])],[
 AC_MSG_RESULT(yes)
 CFLAGS="-fp_port -mieee-fp -wd1572 -wd265 -wd186 -wd239 $CFLAGS"
],[AC_MSG_RESULT(no)])

dnl If compiler is gcc, then use some specific flags.
dnl But don't touch user other flags.
if test "$test_CFLAGS" != set && test -n "$GCC"; then
  CFLAGS="-Wall -Wmissing-prototypes -Wpointer-arith $CFLAGS"
fi

AM_PROG_CC_C_O

case $host in
  *-apple-darwin*)
dnl This allows to take the first GMP library in the library paths,
dnl whether it is dynamic or static. This behavior is more sensible,
dnl in particular because it is the only way to link with a version
dnl only available in static form when another version is available
dnl in dynamic, and also for consistency, because the compiler will
dnl take the first gmp.h found in the include paths (so, we need to
dnl take a library that corresponds to this header file). This is a
dnl common problem with darwin.
    MPFR_LD_SEARCH_PATHS_FIRST ;;
esac

AC_C_CONST
AC_C_VOLATILE
dnl Determine the endianness of integer. But since MPFR doesn't use
dnl AC_CONFIG_HEADERS, configure cannot provide such information for
dnl universal binaries containing variants with different endianness
dnl (i.e. generic code must be used).
AC_C_BIGENDIAN([AC_DEFINE(HAVE_BIG_ENDIAN)],[AC_DEFINE(HAVE_LITTLE_ENDIAN)],
[true],[true])

MPFR_CONFIGS

# (Based on GMP 5.1)
# clock_gettime is in librt on *-*-osf5.1 and on glibc < 2.17, so add -lrt to
# TUNE_LIBS if needed (e.g. if clock_gettime is not already in the C library).
# On linux (tested on x86_32, 2.6.26), clock_getres reports ns accuracy,
# while in a quick test on osf, clock_getres said only 1 millisecond.
old_LIBS="$LIBS"
AC_SEARCH_LIBS(clock_gettime, rt, [
  AC_DEFINE([HAVE_CLOCK_GETTIME],1,[Define to 1 if you have the `clock_gettime' function])])
TUNE_LIBS="$LIBS"
LIBS="$old_LIBS"
AC_SUBST(TUNE_LIBS)

dnl
dnl Setup GMP detection
dnl

dnl Check GMP Header
AC_MSG_CHECKING(for gmp.h)
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include "gmp.h"
]])],[AC_MSG_RESULT(yes)],[
 AC_MSG_RESULT(no)
 AC_MSG_ERROR([gmp.h can't be found, or is unusable.])
])

dnl Configs for Windows DLLs.
dnl libtool requires "-no-undefined" for win32 dll
dnl
dnl "-Wl,output-def" is used to get a .def file for use by MS lib to make
dnl a .lib import library, described in the manual.
dnl
dnl Incidentally, libtool does generate an import library libmpfr.dll.a,
dnl but it's "ar" format and cannot be used by the MS linker.  There
dnl doesn't seem to be any GNU tool for generating or converting to .lib.
AC_SUBST(MPFR_LDFLAGS)
AC_SUBST(LIBMPFR_LDFLAGS)
LT_INIT(win32-dll)
case $host in
  *-*-cygwin* | *-*-mingw* | *-*-pw32* | *-*-os2*)
   AC_MSG_CHECKING(for DLL/static GMP)
   if test "$enable_shared" = yes; then
     MPFR_LDFLAGS="$MPFR_LDFLAGS -no-undefined"
     LIBMPFR_LDFLAGS="$LIBMPFR_LDFLAGS -Wl,--output-def,.libs/libmpfr-4.dll.def"
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include "gmp.h"
#if !__GMP_LIBGMP_DLL
# error "Dead man"
error
#endif
     ]], [[]])],[AC_MSG_RESULT(DLL)],[
  AC_MSG_RESULT(static)
  AC_MSG_ERROR([gmp.h isn't a DLL: use --enable-static --disable-shared]) ])
   else
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include "gmp.h"
#if __GMP_LIBGMP_DLL
# error "Dead man"
error
#endif
     ]], [[]])],[AC_MSG_RESULT(static)],[
  AC_MSG_RESULT(DLL)
  AC_MSG_ERROR([gmp.h is a DLL: use --disable-static --enable-shared]) ])
  fi
  ;;
esac


dnl
dnl For mpfr-longlong.h - TODO: should be replaced (see acinclude.m4).
dnl

GMP_C_ATTRIBUTE_MODE


dnl
dnl Setup GMP detection (continued)
dnl

dnl Check minimal GMP version
dnl We only guarantee that with a *functional* and recent enough GMP version,
dnl MPFR will compile; we do not guarantee that GMP will compile.
dnl In particular fat builds are broken in GMP 4.3.2 and GMP 5.0.0
dnl (at least on 64-bit Core 2 under Linux),
dnl see http://gmplib.org/list-archives/gmp-bugs/2011-August/002345.html.
AC_MSG_CHECKING(for recent GMP)
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include "gmp.h"
#if (__GNU_MP_VERSION*100+__GNU_MP_VERSION_MINOR*10 < 410)
# error "min GMP version is 4.1.0"
error
#endif
]])],[AC_MSG_RESULT(yes)],[
 AC_MSG_RESULT(no)
 AC_MSG_ERROR([GMP 4.1.0 min required])
])

dnl Check if we can use internal header files of GMP (only --with-gmp-build)
if test "$use_gmp_build" = yes ; then
   AC_MSG_CHECKING(for gmp internal files)
   AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
   #include "gmp.h"
   #include "gmp-impl.h"
   #include "longlong.h"
   ]])],[
     AC_MSG_RESULT(yes)
     AC_DEFINE([MPFR_HAVE_GMP_IMPL],1,[Use GMP Internal Files])
   ],[
     AC_MSG_ERROR([header files gmp-impl.h and longlong.h not found])
   ])
fi

dnl Check for valid GMP_NUMB_BITS and BYTES_PER_MP_LIMB
dnl This test doesn't need to link with libgmp (at least it shouldn't).
if test "$use_gmp_build" = yes ; then
  AC_MSG_CHECKING(for valid GMP_NUMB_BITS)
  AC_RUN_IFELSE([AC_LANG_PROGRAM([[
#include <limits.h>
#include "gmp.h"
#include "gmp-impl.h"
]], [[
  return GMP_NUMB_BITS == BYTES_PER_MP_LIMB * CHAR_BIT
         && sizeof(mp_limb_t) == BYTES_PER_MP_LIMB ? 0 : 1;
]])], [AC_MSG_RESULT(yes)], [
       AC_MSG_RESULT(no)
       AC_MSG_ERROR([GMP_NUMB_BITS is incorrect.
You probably need to change some of the GMP or MPFR compile options.])],
       [AC_MSG_RESULT([can't test])])
fi


dnl Check if we can link with GMP
AC_CHECK_LIB(gmp, __gmpz_init, [LIBS="-lgmp $LIBS"],
 [AC_MSG_ERROR(libgmp not found or uses a different ABI (including static vs shared).
Please read the INSTALL file -- see "In case of problem".)])

dnl Check for corresponding 'gmp.h' and libgmp.a
AC_MSG_CHECKING(if gmp.h version and libgmp version are the same)
dnl We do not set LD_LIBRARY_PATH, as it is not possible to set it just
dnl before the test program is run, and we do not want to affect other
dnl programs (such as the compiler), because the behavior could be
dnl incorrect and even have security implications.
dnl WARNING! LD_RUN_PATH is not taken into account by the GNU gold ld,
dnl e.g. from binutils-gold 2.22-5 under Debian; see
dnl   http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=660813
saved_LD_RUN_PATH="$LD_RUN_PATH"
LD_RUN_PATH="${LD_RUN_PATH:+$LD_RUN_PATH$PATH_SEPARATOR}$gmp_lib_path"
export LD_RUN_PATH
AC_RUN_IFELSE([AC_LANG_PROGRAM([[
#include <stdio.h>
#include <string.h>
#include "gmp.h"
]], [[
  char buffer[100];
  sprintf (buffer, "%d.%d.%d", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR,
           __GNU_MP_VERSION_PATCHLEVEL);
  printf ("(%s/%s) ", buffer, gmp_version);
  fflush (stdout);
  /* Also put information in config.log (stderr) */
  fprintf (stderr, "gmp.h: %s / libgmp: %s\n", buffer, gmp_version);
  if (strcmp (buffer, gmp_version) == 0)
    return 0;
  if (__GNU_MP_VERSION_PATCHLEVEL != 0)
    return 1;
  sprintf (buffer, "%d.%d", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR);
  return (strcmp (buffer, gmp_version) != 0) ? 1 : 0;
]])],
   [AC_MSG_RESULT(yes)
    MPFR_CHECK_DBL2INT_BUG
    MPFR_CHECK_PRINTF_SPEC],
   [AC_MSG_RESULT(no)
    AC_MSG_WARN([==========================================================])
    AC_MSG_WARN(['gmp.h' and 'libgmp' seem to have different versions or])
    AC_MSG_WARN([we cannot run a program linked with GMP (if you cannot])
    AC_MSG_WARN([see the version numbers above). A cause may be different])
    AC_MSG_WARN([GMP versions with different ABI's or the use of --with-gmp])
    AC_MSG_WARN([or --with-gmp-include with a system include directory])
    AC_MSG_WARN([(such as /usr/include or /usr/local/include).])
    AC_MSG_WARN([However since we can't use 'libtool' inside the configure,])
    AC_MSG_WARN([we can't be sure. See 'config.log' for details.])
    AC_MSG_WARN([Temporary LD_RUN_PATH was "$LD_RUN_PATH".])
    AC_MSG_WARN([==========================================================])
    ],AC_MSG_RESULT([can not test])
)
LD_RUN_PATH="$saved_LD_RUN_PATH"

dnl Warning! __gmpn_rootrem is an internal GMP symbol; thus its behavior
dnl may change or this symbol may be removed in the future (without being
dnl handled by the library versioning system, which is even worse, as this
dnl can mean undetected incompatibilities in case of GMP library upgrade,
dnl without rebuilding MPFR). So, this symbol must not be used, unless
dnl WANT_GMP_INTERNALS is defined. Only the GMP public API should be used
dnl by default, in particular by binary distributions. Moreover the check
dnl below may yield an incorrect result as libtool isn't used in configure
dnl (see above).
dnl Same for __gmpn_sbpi1_divappr_q.
AC_CHECK_FUNCS([__gmpn_rootrem __gmpn_sbpi1_divappr_q])

dnl Remove also many MACROS (AC_DEFINE) which are unused by MPFR
dnl and polluate (and slow down because libtool has to parse them) the build.
if test -f confdefs.h; then
  $SED '/#define PACKAGE_/d' <confdefs.h >confdefs.tmp
  $SED '/#define HAVE_STRING/d' <confdefs.tmp >confdefs.h
  $SED '/#define HAVE_ALLOCA /d' <confdefs.h >confdefs.tmp
  $SED '/#define HAVE_DLFCN_H/d' <confdefs.tmp >confdefs.h
  $SED '/#define HAVE_MEM/d' <confdefs.h >confdefs.tmp
  $SED '/#define STDC_HEADERS/d' <confdefs.tmp >confdefs.h
  $SED '/#define HAVE_STRTOL/d' <confdefs.h >confdefs.tmp
  $SED '/#define HAVE_STDLIB_H/d' <confdefs.tmp >confdefs.h
  $SED '/#define HAVE_UNISTD_H/d' <confdefs.h >confdefs.tmp
  $SED '/#define HAVE_STDC_HEADERS/d' <confdefs.tmp >confdefs.h
  $SED '/#define HAVE_LONG_DOUBLE/d' <confdefs.h >confdefs.tmp
  $SED '/#define HAVE_SYS_STAT_H/d' <confdefs.tmp >confdefs.h
  $SED '/#define HAVE_SYS_TYPES_H/d' <confdefs.h >confdefs.tmp
  $SED '/#define PROTOTYPES/d' <confdefs.tmp >confdefs.h
  $SED '/#define __PROTOTYPES/d' <confdefs.h >confdefs.tmp

  mv confdefs.tmp confdefs.h
fi

if $EGREP -q -e '-dev$' $srcdir/VERSION; then
  AC_SUBST([DATAFILES])dnl
  DATAFILES=`echo \`$SED -n \
    's/^ *data_check *("\(data\/[[^"]]*\)".*/tests\/\1/p' \
    $srcdir/tests/*.c\``
fi

dnl Output
AC_CONFIG_FILES([Makefile doc/Makefile src/Makefile tests/Makefile tune/Makefile src/mparam.h:src/mparam_h.in])
AC_OUTPUT

dnl NEWS README AUTHORS Changelog
