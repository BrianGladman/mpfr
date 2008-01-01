dnl  MPFR specific autoconf macros

dnl  Copyright 2000, 2002, 2003, 2004, 2005, 2006, 2007, 2008 Free Software Foundation, Inc.
dnl  Contributed by the Arenaire and Cacao projects, INRIA.
dnl
dnl  This file is part of the MPFR Library.
dnl
dnl  The MPFR Library is free software; you can redistribute it and/or modify
dnl  it under the terms of the GNU Lesser General Public License as published
dnl  by the Free Software Foundation; either version 2.1 of the License, or (at
dnl  your option) any later version.
dnl
dnl  The MPFR Library is distributed in the hope that it will be useful, but
dnl  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
dnl  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
dnl  License for more details.
dnl
dnl  You should have received a copy of the GNU Lesser General Public License
dnl  along with the MPFR Library; see the file COPYING.LIB.  If not, write to
dnl  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
dnl  MA 02110-1301, USA.

dnl  autoconf 2.50 is necessary because of the use of AH_VERBATIM,
dnl  but it would be better to make the config file compatible with
dnl  both autoconf 2.13 and autoconf 2.50.
dnl  The following line allows the autoconf wrapper (when installed)
dnl  to work as expected.
AC_PREREQ(2.50)

dnl ------------------------------------------------------------
dnl You must put in MPFR_CONFIGS everything which configure MPFR
dnl except:
dnl   -everything dealing with CC and CFLAGS in particular the ABI
dnl   but the IEEE-754 specific flags must be set here.
dnl   -GMP's linkage.
dnl   -Libtool stuff.
dnl   -Handling of special arguments of MPFR's configure.
AC_DEFUN([MPFR_CONFIGS],
[
AC_REQUIRE([AC_OBJEXT])
AC_REQUIRE([MPFR_CHECK_LIBM])
AC_REQUIRE([AC_HEADER_TIME])
AC_REQUIRE([AC_CANONICAL_HOST])

AC_CHECK_HEADER([limits.h],, AC_MSG_ERROR([limits.h not found]))
AC_CHECK_HEADER([float.h],,  AC_MSG_ERROR([float.h not found]))
AC_CHECK_HEADER([string.h],, AC_MSG_ERROR([string.h not found]))

dnl Check for locales
AC_CHECK_HEADERS([locale.h])

dnl Check for wide characters (wchar_t and wint_t)
AC_CHECK_HEADERS([wchar.h])

dnl Check for stdargs
AC_CHECK_HEADER([stdarg.h],[AC_DEFINE([HAVE_STDARG],1,[Define if stdarg])],
  [AC_CHECK_HEADER([varargs.h],,
    AC_MSG_ERROR([stdarg.h or varargs.h not found]))])

dnl sys/fpu.h - MIPS specific
AC_CHECK_HEADERS([sys/time.h sys/fpu.h])

dnl FIXME: The functions memmove, memset and strtol are really needed by
dnl MPFR, but if they are implemented as macros, this is also OK (in our
dnl case).  So, we do not return an error, but their tests are currently
dnl useless.
dnl gettimeofday is not defined for MinGW
AC_CHECK_FUNCS([memmove memset setlocale strtol gettimeofday])

dnl Check for IEEE-754 switches on Alpha
case $host in
alpha*-*-*)
  saved_CFLAGS="$CFLAGS"
  AC_CACHE_CHECK([for IEEE-754 switches], mpfr_cv_ieee_switches, [
  if test -n "$GCC"; then
    mpfr_cv_ieee_switches="-mfp-rounding-mode=d -mieee-with-inexact"
  else
    mpfr_cv_ieee_switches="-fprm d -ieee_with_inexact"
  fi
  CFLAGS="$CFLAGS $mpfr_cv_ieee_switches"
  AC_TRY_COMPILE(,,, mpfr_cv_ieee_switches="none")
  ])
  if test "$mpfr_cv_ieee_switches" = "none"; then
    CFLAGS="$saved_CFLAGS"
  else
    CFLAGS="$saved_CFLAGS $mpfr_cv_ieee_switches"
  fi
esac

dnl Check for Core2 processor
case $host in
x86_64-*linux*)
  case `sed -n '/^vendor_id/s/^.*: \(.*\)/\1/p' < /proc/cpuinfo` in
    *Intel*) AC_DEFINE(HAVE_HOST_CORE2,1,[Define if processor is Core 2]) ;;
  esac
esac

AC_CHECK_TYPE( [union fpc_csr],
   AC_DEFINE(HAVE_FPC_CSR,1,[Define if union fpc_csr is available]), ,
[
#if HAVE_SYS_FPU_H
#  include <sys/fpu.h>
#endif
])

dnl Check for fesetround
AC_CACHE_CHECK([for fesetround], mpfr_cv_have_fesetround, [
saved_LIBS="$LIBS"
LIBS="$LIBS $MPFR_LIBM"
AC_TRY_LINK([#include <fenv.h>], [fesetround(FE_TONEAREST);],
  mpfr_cv_have_fesetround=yes, mpfr_cv_have_fesetround=no)
LIBS="$saved_LIBS"
])
if test "$mpfr_cv_have_fesetround" = "yes"; then
  AC_DEFINE(MPFR_HAVE_FESETROUND,1,[Define if you have the `fesetround' function via the <fenv.h> header file.])
fi

dnl Check for gcc float-conversion bug; if need be, -ffloat-store is used to
dnl force the conversion to the destination type when a value is stored to
dnl a variable (see ISO C99 standard 5.1.2.3#13, 6.3.1.5#2, 6.3.1.8#2). This
dnl is important concerning the exponent range. Note that this doesn't solve
dnl the double-rounding problem.
if test -n "$GCC"; then
  AC_CACHE_CHECK([for gcc float-conversion bug], mpfr_cv_gcc_floatconv_bug, [
  saved_LIBS="$LIBS"
  LIBS="$LIBS $MPFR_LIBM"
  AC_TRY_RUN([
#include <float.h>
#ifdef MPFR_HAVE_FESETROUND
#include <fenv.h>
#endif
static double get_max (void);
int main() {
  double x = 0.5;
  int i;
  for (i = 1; i <= 11; i++)
    x *= x;
  if (x != 0)
    return 1;
#ifdef MPFR_HAVE_FESETROUND
  /* Useful test for the G4 PowerPC */
  fesetround(FE_TOWARDZERO);
  x = get_max ();
  x *= 2.0;
  if (x != get_max ())
    return 1;
#endif
  return 0;
}
static double get_max (void) { return DBL_MAX; }
  ], [mpfr_cv_gcc_floatconv_bug="no"],
     [mpfr_cv_gcc_floatconv_bug="yes, use -ffloat-store"],
     [mpfr_cv_gcc_floatconv_bug="cannot test, use -ffloat-store"])
  LIBS="$saved_LIBS"
  ])
  if test "$mpfr_cv_gcc_floatconv_bug" != "no"; then
    CFLAGS="$CFLAGS -ffloat-store"
  fi
fi

dnl Check if denormalized numbers are supported
AC_CACHE_CHECK([for denormalized numbers], mpfr_cv_have_denorms, [
AC_TRY_RUN([
#include <math.h>
#include <stdio.h>
int main() {
  double x = 2.22507385850720138309e-308;
  fprintf (stderr, "%e\n", x / 2.0);
  return 2.0 * (x / 2.0) != x;
}
], mpfr_cv_have_denorms=yes, mpfr_cv_have_denorms=no, mpfr_cv_have_denorms=no)
])
if test "$mpfr_cv_have_denorms" = "yes"; then
  AC_DEFINE(HAVE_DENORMS,1,[Define if denormalized floats work.])
fi

dnl Check whether NAN != NAN (as required by the IEEE-754 standard,
dnl but not by the ISO C standard). For instance, this is false with
dnl MIPSpro 7.3.1.3m under IRIX64. By default, assume this is true.
AC_CACHE_CHECK([if NAN == NAN], mpfr_cv_nanisnan, [
AC_TRY_RUN([
#include <stdio.h>
#include <math.h>
#ifndef NAN
# define NAN (0.0/0.0)
#endif
int main() {
  double d;
  d = NAN;
  return d != d;
}
], [mpfr_cv_nanisnan="yes"],
   [mpfr_cv_nanisnan="no"],
   [mpfr_cv_nanisnan="cannot test, assume no"])
])
if test "$mpfr_cv_nanisnan" = "yes"; then
  AC_DEFINE(MPFR_NANISNAN,1,[Define if NAN == NAN.])
fi

dnl Check if the chars '0' to '9' are consecutive values
AC_MSG_CHECKING([if charset has consecutive values])
AC_RUN_IFELSE(AC_LANG_PROGRAM([[
char *number = "0123456789";
char *lower  = "abcdefghijklmnopqrstuvwxyz";
char *upper  = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
]],[[
 int i;
 unsigned char *p;
 for (p = (unsigned char*) number, i = 0; i < 9; i++)
   if ( (*p)+1 != *(p+1) ) return 1;
 for (p = (unsigned char*) lower, i = 0; i < 25; i++)
   if ( (*p)+1 != *(p+1) ) return 1;
 for (p = (unsigned char*) upper, i = 0; i < 25; i++)
   if ( (*p)+1 != *(p+1) ) return 1;
]]), [AC_MSG_RESULT(yes)],[
 AC_MSG_RESULT(no)
 AC_DEFINE(MPFR_NO_CONSECUTIVE_CHARSET,1,[Charset is not consecutive])
], [AC_MSG_RESULT(can not test)])

dnl Must be checked with the LIBM
dnl but we don't want to add the LIBM to MPFR dependency.
dnl Can't use AC_CHECK_FUNCS since the function may be in LIBM but
dnl not exported in math.h
saved_LIBS="$LIBS"
LIBS="$LIBS $MPFR_LIBM"
dnl AC_CHECK_FUNCS([round trunc floor ceil nearbyint])
AC_MSG_CHECKING(for math/round)
AC_LINK_IFELSE([AC_LANG_PROGRAM([[
#include <math.h>
int f (double (*func)(double)) { return 0;}
]], [[
 double a = 17.42;
 a = f (round);
 return 0;
]])], [
   AC_MSG_RESULT(yes)
   AC_DEFINE(HAVE_ROUND, 1,[Have ISO-C99 round function])
],[AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for math/trunc)
AC_LINK_IFELSE([AC_LANG_PROGRAM([[
#include <math.h>
int f (double (*func)(double)) { return 0;}
]], [[
 double a = 17.42;
 a = f(trunc);
 return 0;
]])], [
   AC_MSG_RESULT(yes)
   AC_DEFINE(HAVE_TRUNC, 1,[Have ISO-C99 trunc function])
],[AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for math/floor)
AC_LINK_IFELSE([AC_LANG_PROGRAM([[
#include <math.h>
int f (double (*func)(double)) { return 0;}
]], [[
 double a = 17.42;
 a = f(floor);
 return 0;
]])], [
   AC_MSG_RESULT(yes)
   AC_DEFINE(HAVE_FLOOR, 1,[Have ISO-C99 floor function])
],[AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for math/ceil)
AC_LINK_IFELSE([AC_LANG_PROGRAM([[
#include <math.h>
int f (double (*func)(double)) { return 0;}
]], [[
 double a = 17.42;
 a = f(ceil);
 return 0;
]])], [
   AC_MSG_RESULT(yes)
   AC_DEFINE(HAVE_CEIL, 1,[Have ISO-C99 ceil function])
],[AC_MSG_RESULT(no)])

AC_MSG_CHECKING(for math/rint)
AC_LINK_IFELSE([AC_LANG_PROGRAM([[
#include <math.h>
int f (double (*func)(double)) { return 0;}
]], [[
 double a = 17.42;
 a = f(nearbyint);
 return 0;
]])], [
   AC_MSG_RESULT(yes)
   AC_DEFINE(HAVE_NEARBYINT, 1,[Have ISO-C99 rint function])
],[AC_MSG_RESULT(no)])

LIBS="$saved_LIBS"

dnl Now try to check the long double format
MPFR_C_LONG_DOUBLE_FORMAT
])


dnl  MPFR_C_LONG_DOUBLE_FORMAT
dnl  -------------------------
dnl  Determine the format of a long double.
dnl
dnl  The object file is grepped, so as to work when cross compiling.  A
dnl  start and end sequence is included to avoid false matches, and
dnl  allowance is made for the desired data crossing an "od -b" line
dnl  boundary.  The test number is a small integer so it should appear
dnl  exactly, no rounding or truncation etc.
dnl
dnl  "od -b" is supported even by Unix V7, and the awk script used doesn't
dnl  have functions or anything, so even an "old" awk should suffice.
dnl
dnl  The 10-byte IEEE extended format is generally padded to either 12 or 16
dnl  bytes for alignment purposes.  The SVR4 i386 ABI is 12 bytes, or i386
dnl  gcc -m128bit-long-double selects 16 bytes.  IA-64 is 16 bytes in LP64
dnl  mode, or 12 bytes in ILP32 mode.  The first 10 bytes is the relevant
dnl  part in all cases (big and little endian).
dnl
dnl  Enhancements:
dnl
dnl  Could match more formats, but no need to worry until there's code
dnl  wanting to use them.
dnl
dnl  Don't want to duplicate the double matching from GMP_C_DOUBLE_FORMAT,
dnl  perhaps we should merge with that macro, to match data formats
dnl  irrespective of the C type in question.  Or perhaps just let the code
dnl  use DOUBLE macros when sizeof(double)==sizeof(long double).

AC_DEFUN([MPFR_C_LONG_DOUBLE_FORMAT],
[AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_PROG_AWK])
AC_REQUIRE([AC_OBJEXT])
AC_CHECK_TYPES([long double])
AC_CACHE_CHECK([format of `long double' floating point],
                mpfr_cv_c_long_double_format,
[mpfr_cv_c_long_double_format=unknown
if test "$ac_cv_type_long_double" != yes; then
  mpfr_cv_c_long_double_format="not available"
else
  cat >conftest.c <<\EOF
[
/* "before" is 16 bytes to ensure there's no padding between it and "x".
   We're not expecting any "long double" bigger than 16 bytes or with
   alignment requirements stricter than 16 bytes.  */
struct {
  char         before[16];
  long double  x;
  char         after[8];
} foo = {
  { '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\001', '\043', '\105', '\147', '\211', '\253', '\315', '\357' },
  -123456789.0,
  { '\376', '\334', '\272', '\230', '\166', '\124', '\062', '\020' }
};
]
EOF
  mpfr_compile="$CC $CFLAGS $CPPFLAGS -c conftest.c >&AC_FD_CC 2>&1"
  if AC_TRY_EVAL(mpfr_compile); then
    cat >conftest.awk <<\EOF
[
BEGIN {
  found = 0
}

# got[] holds a sliding window of bytes read the input.  got[0] is the most
# recent byte read, and got[31] the oldest byte read, so when looking to
# match some data the indices are "reversed".
#
{
  for (f = 2; f <= NF; f++)
    {
      # new byte, shift others up
      for (i = 31; i >= 0; i--)
        got[i+1] = got[i];
      got[0] = $f;

      # end sequence
      if (got[7] != "376") continue
      if (got[6] != "334") continue
      if (got[5] != "272") continue
      if (got[4] != "230") continue
      if (got[3] != "166") continue
      if (got[2] != "124") continue
      if (got[1] != "062") continue
      if (got[0] != "020") continue

      # start sequence, with 12-byte body
      if (got[27] == "001" && \
          got[26] == "043" && \
          got[25] == "105" && \
          got[24] == "147" && \
          got[23] == "211" && \
          got[22] == "253" && \
          got[21] == "315" && \
          got[20] == "357")
        {
          saw = " (" got[19] \
                 " " got[18] \
                 " " got[17] \
                 " " got[16] \
                 " " got[15] \
                 " " got[14] \
                 " " got[13] \
                 " " got[12] \
                 " " got[11] \
                 " " got[10] \
                 " " got[9]  \
                 " " got[8] ")"

          if (got[19] == "000" && \
              got[18] == "000" && \
              got[17] == "000" && \
              got[16] == "000" && \
              got[15] == "240" && \
              got[14] == "242" && \
              got[13] == "171" && \
              got[12] == "353" && \
              got[11] == "031" && \
              got[10] == "300")
            {
              print "IEEE extended, little endian"
              found = 1
              exit
            }
        }

      # start sequence, with 16-byte body
      if (got[31] == "001" && \
          got[30] == "043" && \
          got[29] == "105" && \
          got[28] == "147" && \
          got[27] == "211" && \
          got[26] == "253" && \
          got[25] == "315" && \
          got[24] == "357")
        {
          saw = " (" got[23] \
                 " " got[22] \
                 " " got[21] \
                 " " got[20] \
                 " " got[19] \
                 " " got[18] \
                 " " got[17] \
                 " " got[16] \
                 " " got[15] \
                 " " got[14] \
                 " " got[13] \
                 " " got[12] \
                 " " got[11] \
                 " " got[10] \
                 " " got[9]  \
                 " " got[8] ")"

          if (got[23] == "000" && \
              got[22] == "000" && \
              got[21] == "000" && \
              got[20] == "000" && \
              got[19] == "240" && \
              got[18] == "242" && \
              got[17] == "171" && \
              got[16] == "353" && \
              got[15] == "031" && \
              got[14] == "300")
            {
              print "IEEE extended, little endian"
              found = 1
              exit
            }

          if (got[23] == "300" && \
              got[22] == "031" && \
              got[21] == "326" && \
              got[20] == "363" && \
              got[19] == "105" && \
              got[18] == "100" && \
              got[17] == "000" && \
              got[16] == "000" && \
              got[15] == "000" && \
              got[14] == "000" && \
              got[13] == "000" && \
              got[12] == "000" && \
              got[11] == "000" && \
              got[10] == "000" && \
              got[9]  == "000" && \
              got[8]  == "000")
            {
              print "IEEE quad, big endian"
              found = 1
              exit
            }
        }
    }
}

END {
  if (! found)
    print "unknown", saw
}
]
EOF
    mpfr_cv_c_long_double_format=`od -b conftest.$OBJEXT | $AWK -f conftest.awk`
    case $mpfr_cv_c_long_double_format in
    unknown*)
      echo "cannot match anything, conftest.$OBJEXT contains" >&AC_FD_CC
      od -b conftest.$OBJEXT >&AC_FD_CC
      ;;
    esac
  else
    AC_MSG_WARN([oops, cannot compile test program])
  fi
fi
rm -f conftest*
])

AH_VERBATIM([HAVE_LDOUBLE],
[/* Define one of the following to 1 for the format of a `long double'.
   If your format is not among these choices, or you don't know what it is,
   then leave all undefined.
   IEEE_EXT is the 10-byte IEEE extended precision format.
   IEEE_QUAD is the 16-byte IEEE quadruple precision format.
   LITTLE or BIG is the endianness.  */
#undef HAVE_LDOUBLE_IEEE_EXT_LITTLE
#undef HAVE_LDOUBLE_IEEE_QUAD_BIG])

case $mpfr_cv_c_long_double_format in
  "IEEE extended, little endian")
    AC_DEFINE(HAVE_LDOUBLE_IEEE_EXT_LITTLE, 1)
    ;;
  "IEEE quad, big endian")
    AC_DEFINE(HAVE_LDOUBLE_IEEE_QUAD_BIG, 1)
    ;;
  unknown* | "not available")
    ;;
  *)
    AC_MSG_WARN([oops, unrecognised float format: $mpfr_cv_c_long_double_format])
    ;;
esac
])


dnl  MPFR_CHECK_LIBM
dnl  ---------------
dnl  Determine a math library -lm to use.

AC_DEFUN([MPFR_CHECK_LIBM],
[AC_REQUIRE([AC_CANONICAL_HOST])
AC_SUBST(MPFR_LIBM,'')
case $host in
  *-*-beos* | *-*-cygwin* | *-*-pw32*)
    # According to libtool AC CHECK LIBM, these systems don't have libm
    ;;
  *-*-hpux*)
    # -lM means something subtly different to -lm, SVID style error handling
    # or something.  FIXME: Why exactly do we want this?
    AC_CHECK_LIB(M, main, MPFR_LIBM="-lM")
    ;;
  *-*-solaris*)
    # On Solaris the math functions new in C99 are in -lm9x.
    # FIXME: Do we need -lm9x as well as -lm, or just instead of?
    AC_CHECK_LIB(m9x, main, MPFR_LIBM="-lm9x")
    AC_CHECK_LIB(m,   main, MPFR_LIBM="$MPFR_LIBM -lm")
    ;;
  *-ncr-sysv4.3*)
    # FIXME: What does -lmw mean?  Libtool AC CHECK LIBM does it this way.
    AC_CHECK_LIB(mw, _mwvalidcheckl, MPFR_LIBM="-lmw")
    AC_CHECK_LIB(m, main, MPFR_LIBM="$MPFR_LIBM -lm")
    ;;
  *)
    AC_CHECK_LIB(m, main, MPFR_LIBM="-lm")
    ;;
esac
])


dnl  MPFR_LD_SEARCH_PATHS_FIRST
dnl  --------------------------

AC_DEFUN([MPFR_LD_SEARCH_PATHS_FIRST],
[case "$LD $LDFLAGS" in
  *-Wl,-search_paths_first*) ;;
  *) AC_MSG_CHECKING([if the compiler understands -Wl,-search_paths_first])
     saved_LDFLAGS="$LDFLAGS"
     LDFLAGS="-Wl,-search_paths_first $LDFLAGS"
     AC_LINK_IFELSE([AC_LANG_PROGRAM([[]], [[]])],
       [AC_MSG_RESULT(yes)],
       [AC_MSG_RESULT(no)]
        LDFLAGS="$saved_LDFLAGS")
     ;;
 esac
])


dnl  GMP_C_ATTRIBUTE_MODE
dnl  --------------------
dnl  Introduced in gcc 2.2, but perhaps not in all Apple derived versions.
dnl  Needed for mpfr-longlong.h; this is currently necessary for s390.
dnl
dnl  TODO: Replace this with a cleaner type size detection, as this
dnl  solution only works with gcc and assumes CHAR_BIT == 8. Probably use
dnl  <stdint.h>, and <http://gcc.gnu.org/viewcvs/trunk/config/stdint.m4>
dnl  as a fallback.

AC_DEFUN([GMP_C_ATTRIBUTE_MODE],
[AC_CACHE_CHECK([whether gcc __attribute__ ((mode (XX))) works],
               gmp_cv_c_attribute_mode,
[AC_TRY_COMPILE([typedef int SItype __attribute__ ((mode (SI)));], ,
 gmp_cv_c_attribute_mode=yes, gmp_cv_c_attribute_mode=no)
])
if test $gmp_cv_c_attribute_mode = yes; then
 AC_DEFINE(HAVE_ATTRIBUTE_MODE, 1,
 [Define to 1 if the compiler accepts gcc style __attribute__ ((mode (XX)))])
fi
])
