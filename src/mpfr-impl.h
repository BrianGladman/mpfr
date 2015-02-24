/* Utilities for MPFR developers, not exported.

Copyright 1999-2015 Free Software Foundation, Inc.
Contributed by the AriC and Caramel projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#ifndef __MPFR_IMPL_H__
#define __MPFR_IMPL_H__

/* Let's include some standard headers unconditionally as they are
   already needed by several source files or when some options are
   enabled/disabled, and it is easy to forget them (some configure
   options may hide the error).
   Note: If some source file must not have such a header included
   (which is very unlikely and probably means something broken in
   this source file), we should do that with some macro (that would
   also force to disable incompatible features). */
#if defined (__cplusplus)
#include <cstdio>
#include <cstring>
#else
#include <stdio.h>
#include <string.h>
#endif
#include <stdlib.h>
#include <limits.h>

/* The macros defined in mpfr-cvers.h do not depend on anything,
   so that it is better to include this header file early: then
   it can be used by any other header. */
#include "mpfr-cvers.h"

#if _MPFR_EXP_FORMAT == 4
/* mpfr_exp_t will be defined as intmax_t */
# include "mpfr-intmax.h"
#endif

/* Check if we are inside a build of MPFR or inside the test suite.
   This is needed in mpfr.h to export or import the functions.
   It matters only for Windows DLL */
#ifndef __MPFR_TEST_H__
# define __MPFR_WITHIN_MPFR 1
#endif



/******************************************************
 ****************** Include files *********************
 ******************************************************/

/* Include 'config.h' before using ANY configure macros if needed
   NOTE: It isn't MPFR 'config.h', but GMP's one! */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#ifdef  MPFR_HAVE_GMP_IMPL /* Build with gmp internals */

# include "gmp.h"
# include "gmp-impl.h"
# ifdef MPFR_NEED_LONGLONG_H
#  include "longlong.h"
# endif
# include "mpfr.h"
# include "mpfr-gmp.h"

#else /* Build without gmp internals */

# include "gmp.h"
/* if using mini-gmp, include missing definitions in mini-gmp */
# ifdef MPFR_USE_MINI_GMP
#  include "mpfr-mini-gmp.h"
# endif
# include "mpfr.h"
# include "mpfr-gmp.h"
# ifdef MPFR_NEED_LONGLONG_H
#  define LONGLONG_STANDALONE
#  include "mpfr-longlong.h"
# endif

#endif

#undef MPFR_NEED_LONGLONG_H



/******************************************************
 ****************** (U)INTMAX_MAX *********************
 ******************************************************/

/* Let's try to fix UINTMAX_MAX and INTMAX_MAX if these macros don't work
   (e.g. with gcc -ansi -pedantic-errors in 32-bit mode under GNU/Linux),
   see <https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=582698>. */
#ifdef _MPFR_H_HAVE_INTMAX_T
# ifdef MPFR_HAVE_INTMAX_MAX
#  define MPFR_UINTMAX_MAX UINTMAX_MAX
#  define MPFR_INTMAX_MAX INTMAX_MAX
#  define MPFR_INTMAX_MIN INTMAX_MIN
# else
#  define MPFR_UINTMAX_MAX ((uintmax_t) -1)
#  define MPFR_INTMAX_MAX ((intmax_t) (MPFR_UINTMAX_MAX >> 1))
#  define MPFR_INTMAX_MIN (INT_MIN + INT_MAX - MPFR_INTMAX_MAX)
# endif
#endif

#define MPFR_BYTES_PER_MP_LIMB (GMP_NUMB_BITS/CHAR_BIT)

/******************************************************
 ************** Attributes definition *****************
 ******************************************************/

/* For the definition of MPFR_THREAD_ATTR. GCC/ICC detection macros are
   no longer used, as they sometimes gave incorrect information about
   the support of thread-local variables. A configure check is now done. */
#include "mpfr-thread.h"

#if defined(MPFR_HAVE_NORETURN)
/* _Noreturn is specified by ISO C11 (Section 6.7.4);
   in GCC, it is supported as of version 4.7. */
# define MPFR_NORETURN _Noreturn
#elif !defined(noreturn)
/* A noreturn macro could be defined if <stdnoreturn.h> has been included,
   in which case it would make sense to #define MPFR_NORETURN noreturn.
   But this is unlikely, as MPFR_HAVE_NORETURN should have been defined
   in such a case. So, in doubt, let us avoid any code that would use a
   noreturn macro, since it could be invalid. */
# if __MPFR_GNUC(3,0) || __MPFR_ICC(8,1,0)
#  define MPFR_NORETURN __attribute__ ((noreturn))
# elif defined(_MSC_VER) && defined(_WIN32) && (_MSC_VER >= 1200)
#  define MPFR_NORETURN __declspec (noreturn)
# endif
#endif
#ifndef MPFR_NORETURN
# define MPFR_NORETURN
#endif

#if __MPFR_GNUC(3,0) || __MPFR_ICC(8,1,0)
# define MPFR_CONST_ATTR    __attribute__ ((const))
#else
# define MPFR_CONST_ATTR
#endif

#if __MPFR_GNUC(3,0) || __MPFR_ICC(8,1,0)
# define MPFR_PURE_FUNCTION_ATTR    __attribute__ ((pure))
#else
# define MPFR_PURE_FUNCTION_ATTR
#endif

/* The hot attribute on a function is used to inform the compiler
   that the function is a hot spot of the compiled program. */
#if __MPFR_GNUC(4,3)
# define MPFR_HOT_FUNCTION_ATTR     __attribute__ ((hot))
#else
# define MPFR_HOT_FUNCTION_ATTR
#endif

/* The cold attribute on functions is used to inform the compiler
   that the function is unlikely to be executed. */
#if __MPFR_GNUC(4,3)
# define MPFR_COLD_FUNCTION_ATTR    __attribute__ ((cold))
#else
# define MPFR_COLD_FUNCTION_ATTR
#endif

/* add MPFR_MAYBE_UNUSED after a variable declaration to avoid compiler
   warnings if it is not used */
#if __MPFR_GNUC(3,4)
#define MPFR_MAYBE_UNUSED __attribute__ ((unused))
#else
#define MPFR_MAYBE_UNUSED
#endif



/******************************************************
 ************* Global Internal Variables **************
 ******************************************************/

#if defined (__cplusplus)
extern "C" {
#endif

/* Cache struct */
struct __gmpfr_cache_s {
  mpfr_t x;
  int inexact;
  int (*func)(mpfr_ptr, mpfr_rnd_t);
};
typedef struct __gmpfr_cache_s mpfr_cache_t[1];
typedef struct __gmpfr_cache_s *mpfr_cache_ptr;

__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_flags_t __gmpfr_flags;
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_exp_t   __gmpfr_emin;
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_exp_t   __gmpfr_emax;
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_prec_t  __gmpfr_default_fp_bit_precision;
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_rnd_t   __gmpfr_default_rounding_mode;
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_cache_t __gmpfr_cache_const_euler;
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_cache_t __gmpfr_cache_const_catalan;

#ifndef MPFR_USE_LOGGING
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_cache_t __gmpfr_cache_const_pi;
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_cache_t __gmpfr_cache_const_log2;
#else
/* Two constants are used by the logging functions (via mpfr_fprintf,
   then mpfr_log, for the base conversion): pi and log(2). Since the
   mpfr_cache function isn't re-entrant when working on the same cache,
   we need to define two caches for each constant. */
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_cache_t __gmpfr_normal_pi;
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_cache_t __gmpfr_normal_log2;
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_cache_t __gmpfr_logging_pi;
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_cache_t __gmpfr_logging_log2;
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_cache_ptr __gmpfr_cache_const_pi;
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_cache_ptr __gmpfr_cache_const_log2;
#endif

#define BASE_MAX 62
__MPFR_DECLSPEC extern const __mpfr_struct __gmpfr_l2b[BASE_MAX-1][2];

/* Note: do not use the following values when they can be outside the
   current exponent range, e.g. when the exponent range has not been
   extended yet; under such a condition, they can be used only in
   mpfr_cmpabs. */
__MPFR_DECLSPEC extern const mpfr_t __gmpfr_one;
__MPFR_DECLSPEC extern const mpfr_t __gmpfr_two;
__MPFR_DECLSPEC extern const mpfr_t __gmpfr_four;
__MPFR_DECLSPEC extern const mpfr_t __gmpfr_mone;
__MPFR_DECLSPEC extern const mpfr_t __gmpfr_const_log2_RNDD;
__MPFR_DECLSPEC extern const mpfr_t __gmpfr_const_log2_RNDU;


#if defined (__cplusplus)
 }
#endif

/* Replace some common functions for direct access to the global vars.
   The casts prevent these macros from being used as a lvalue (and this
   method makes sure that the expressions have the correct type). */
#define mpfr_get_emin() ((mpfr_exp_t) __gmpfr_emin)
#define mpfr_get_emax() ((mpfr_exp_t) __gmpfr_emax)
#define mpfr_get_default_rounding_mode() \
  ((mpfr_rnd_t) __gmpfr_default_rounding_mode)
#define mpfr_get_default_prec() \
  ((mpfr_prec_t) __gmpfr_default_fp_bit_precision)

/* Flags related macros. */
/* Note: Function-like macros that modify __gmpfr_flags are not defined
   because of the risk to break the sequence point rules if two such
   macros are used in the same expression (without a sequence point
   between). For instance, mpfr_sgn currently uses mpfr_set_erangeflag,
   which mustn't be implemented as a macro for this reason. */

#define mpfr_flags_test(mask) \
  (__gmpfr_flags & (mpfr_flags_t) (mask))

#if MPFR_FLAGS_ALL <= INT_MAX
#define mpfr_underflow_p() \
  ((int) (__gmpfr_flags & MPFR_FLAGS_UNDERFLOW))
#define mpfr_overflow_p() \
  ((int) (__gmpfr_flags & MPFR_FLAGS_OVERFLOW))
#define mpfr_nanflag_p() \
  ((int) (__gmpfr_flags & MPFR_FLAGS_NAN))
#define mpfr_inexflag_p() \
  ((int) (__gmpfr_flags & MPFR_FLAGS_INEXACT))
#define mpfr_erangeflag_p() \
  ((int) (__gmpfr_flags & MPFR_FLAGS_ERANGE))
#define mpfr_divby0_p() \
  ((int) (__gmpfr_flags & MPFR_FLAGS_DIVBY0))
#endif

/* Use a do-while statement for the following macros in order to prevent
   one from using them in an expression, as the sequence point rules could
   be broken if __gmpfr_flags is assigned twice in the same expression
   (via macro expansions). For instance, the mpfr_sgn macro currently uses
   mpfr_set_erangeflag, which mustn't be implemented as a function-like
   macro for this reason. It is not clear whether an expression with
   sequence points, like (void) (0, __gmpfr_flags = 0), would avoid UB. */
#define MPFR_CLEAR_FLAGS() \
  do __gmpfr_flags = 0; while (0)
#define MPFR_CLEAR_UNDERFLOW() \
  do __gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_UNDERFLOW; while (0)
#define MPFR_CLEAR_OVERFLOW() \
  do __gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_OVERFLOW; while (0)
#define MPFR_CLEAR_DIVBY0() \
  do __gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_DIVBY0; while (0)
#define MPFR_CLEAR_NANFLAG() \
  do __gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_NAN; while (0)
#define MPFR_CLEAR_INEXFLAG() \
  do __gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_INEXACT; while (0)
#define MPFR_CLEAR_ERANGEFLAG() \
  do __gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_ERANGE; while (0)
#define MPFR_SET_UNDERFLOW() \
  do __gmpfr_flags |= MPFR_FLAGS_UNDERFLOW; while (0)
#define MPFR_SET_OVERFLOW() \
  do __gmpfr_flags |= MPFR_FLAGS_OVERFLOW; while (0)
#define MPFR_SET_DIVBY0() \
  do __gmpfr_flags |= MPFR_FLAGS_DIVBY0; while (0)
#define MPFR_SET_NANFLAG() \
  do __gmpfr_flags |= MPFR_FLAGS_NAN; while (0)
#define MPFR_SET_INEXFLAG() \
  do __gmpfr_flags |= MPFR_FLAGS_INEXACT; while (0)
#define MPFR_SET_ERANGEFLAG() \
  do __gmpfr_flags |= MPFR_FLAGS_ERANGE; while (0)

/* Testing an exception flag correctly is tricky. There are mainly two
   pitfalls: First, one needs to remember to clear the corresponding
   flag, in case it was set before the function call or during some
   intermediate computations (in practice, one can clear all the flags).
   Secondly, one needs to test the flag early enough, i.e. before it
   can be modified by another function. Moreover, it is quite difficult
   (if not impossible) to reliably check problems with "make check". To
   avoid these pitfalls, it is recommended to use the following macros.
   Other use of the exception-flag predicate functions/macros will be
   detected by mpfrlint.
   Note: _op can be either a statement or an expression.
   MPFR_BLOCK_EXCEP should be used only inside a block; it is useful to
   detect some exception in order to exit the block as soon as possible. */
#define MPFR_BLOCK_DECL(_flags) mpfr_flags_t _flags
/* The (void) (_flags) makes sure that _flags is read at least once (it
   makes sense to use MPFR_BLOCK while _flags will never be read in the
   source, so that we wish to avoid the corresponding warning). */
#define MPFR_BLOCK(_flags,_op)          \
  do                                    \
    {                                   \
      MPFR_CLEAR_FLAGS ();              \
      _op;                              \
      (_flags) = __gmpfr_flags;         \
      (void) (_flags);                  \
    }                                   \
  while (0)
#define MPFR_BLOCK_TEST(_flags,_f) MPFR_UNLIKELY ((_flags) & (_f))
#define MPFR_BLOCK_EXCEP (__gmpfr_flags & (MPFR_FLAGS_UNDERFLOW | \
                                           MPFR_FLAGS_OVERFLOW | \
                                           MPFR_FLAGS_DIVBY0 | \
                                           MPFR_FLAGS_NAN))
/* Let's use a MPFR_ prefix, because e.g. OVERFLOW is defined by glibc's
   math.h, though this is not a reserved identifier! */
#define MPFR_UNDERFLOW(_flags)  MPFR_BLOCK_TEST (_flags, MPFR_FLAGS_UNDERFLOW)
#define MPFR_OVERFLOW(_flags)   MPFR_BLOCK_TEST (_flags, MPFR_FLAGS_OVERFLOW)
#define MPFR_NANFLAG(_flags)    MPFR_BLOCK_TEST (_flags, MPFR_FLAGS_NAN)
#define MPFR_INEXFLAG(_flags)   MPFR_BLOCK_TEST (_flags, MPFR_FLAGS_INEXACT)
#define MPFR_ERANGEFLAG(_flags) MPFR_BLOCK_TEST (_flags, MPFR_FLAGS_ERANGE)
#define MPFR_DIVBY0(_flags)     MPFR_BLOCK_TEST (_flags, MPFR_FLAGS_DIVBY0)


/******************************************************
 ******************** Assertions **********************
 ******************************************************/

/* MPFR_WANT_ASSERT can take 4 values (the default value is 0):
   -1 (or below): Do not check any assertion. Discouraged, in particular
     for a shared library (for time-critical applications, LTO with a
     static library should also be used anyway).
   0: Check normal assertions.
   1: Check debugging assertions too.
   2 (or above): Additional checks that may take time. For instance,
     some functions may be tested by using two different implementations
     and comparing the results.
*/

/* Note: do not use GMP macros ASSERT_ALWAYS and ASSERT as they are not
   expressions, and as a consequence, they cannot be used in a for(),
   with a comma operator and so on. */

/* MPFR_ASSERTN(expr): assertions that should normally be checked,
     otherwise give a hint to the compiler.
   MPFR_ASSERTD(expr): assertions that should be checked when testing,
     otherwise give a hint to the compiler.
   MPFR_DBGRES(assignment): to be used when the result is tested only
     in an MPFR_ASSERTD expression (in order to avoid a warning, e.g.
     with GCC's -Wunused-but-set-variable, in non-debug mode).
 */
#ifndef MPFR_WANT_ASSERT
# define MPFR_WANT_ASSERT 0
#endif

#if MPFR_WANT_ASSERT < 0
# undef MPFR_EXP_CHECK
# define MPFR_ASSERTN(expr)  MPFR_ASSUME (expr)
#else
# define MPFR_ASSERTN(expr)  \
  ((void) ((MPFR_LIKELY(expr)) || (ASSERT_FAIL(expr),0)))
#endif

#if MPFR_WANT_ASSERT > 0
# define MPFR_EXP_CHECK 1
# define MPFR_ASSERTD(expr)  MPFR_ASSERTN (expr)
# define MPFR_DBGRES(A)      (A)
#else
# define MPFR_ASSERTD(expr)  MPFR_ASSUME (expr)
# define MPFR_DBGRES(A)      ((void) (A))
#endif

/* MPFR_ASSUME is like assert(), but it is a hint to a compiler about a
   statement of fact in a function call free expression, which allows
   the compiler to generate better machine code.
   __builtin_unreachable has been introduced in GCC 4.5 but it works
   fine since 4.8 only (before it may generate unoptimized code if there
   are more than one decision).
   The use of __builtin_constant_p is to avoid calling any function in
   assertions: it detects if the expression can be reduced to a constant,
   i.e. if the expression has any other effect than being false or true
   which is wrong if there is any function call.
   In the development code, it is better to use MPFR_ASSERTD than
   MPFR_ASSUME, since it'll check if the condition is true in debug
   build.
*/
#if defined(MPFR_HAVE_BUILTIN_UNREACHABLE) && __MPFR_GNUC(4, 8)
# define MPFR_ASSUME(x)                                 \
    (! __builtin_constant_p (!!(x) || !(x)) || (x) ?    \
     (void) 0 : __builtin_unreachable())
#elif defined(_MSC_VER)
# define MPFR_ASSUME(x) __assume(x)
#else
# define MPFR_ASSUME(x) ((void) 0)
#endif

#include "mpfr-sassert.h"

/* Code to deal with impossible, for functions returning an int.
   The "return 0;" avoids an error with current GCC versions and
   "-Werror=return-type".
   WARNING: It doesn't use do { } while (0) for Insure++ */
#if defined(HAVE_BUILTIN_UNREACHABLE)
# define MPFR_RET_NEVER_GO_HERE() { __builtin_unreachable(); }
#else
# define MPFR_RET_NEVER_GO_HERE()  { MPFR_ASSERTN(0); return 0; }
#endif

/******************************************************
 ******************** Warnings ************************
 ******************************************************/

/* MPFR_WARNING is no longer useful, but let's keep the macro in case
   it needs to be used again in the future. */

#ifdef MPFR_USE_WARNINGS
# define MPFR_WARNING(W)                    \
  do                                        \
    {                                       \
      char *q = getenv ("MPFR_QUIET");      \
      if (q == NULL || *q == 0)             \
        fprintf (stderr, "MPFR: %s\n", W);  \
    }                                       \
  while (0)
#else
# define MPFR_WARNING(W)  ((void) 0)
#endif


/******************************************************
 ****************** double macros *********************
 ******************************************************/

/* Precision used for lower precision computations */
#define MPFR_SMALL_PRECISION 32

/* Definition of constants */
#define LOG2 0.69314718055994528622 /* log(2) rounded to zero on 53 bits */

/* MPFR_DOUBLE_SPEC = 1 if the C type 'double' corresponds to IEEE-754
   double precision, 0 if it doesn't, and undefined if one doesn't know.
   On all the tested machines, MPFR_DOUBLE_SPEC = 1. To have this macro
   defined here, #include <float.h> is needed. If need be, other values
   could be defined for other specs (once they are known). */
#if !defined(MPFR_DOUBLE_SPEC) && defined(FLT_RADIX) && \
    defined(DBL_MANT_DIG) && defined(DBL_MIN_EXP) && defined(DBL_MAX_EXP)
# if FLT_RADIX == 2 && DBL_MANT_DIG == 53 && \
     DBL_MIN_EXP == -1021 && DBL_MAX_EXP == 1024
#  define MPFR_DOUBLE_SPEC 1
# else
#  define MPFR_DOUBLE_SPEC 0
# endif
#endif

/* Debug non IEEE floats */
#ifdef XDEBUG
# undef _GMP_IEEE_FLOATS
#endif
#ifndef _GMP_IEEE_FLOATS
# define _GMP_IEEE_FLOATS 0
#endif

#ifndef IEEE_DBL_MANT_DIG
#define IEEE_DBL_MANT_DIG 53
#endif
#define MPFR_LIMBS_PER_DOUBLE ((IEEE_DBL_MANT_DIG-1)/GMP_NUMB_BITS+1)

#ifndef IEEE_FLT_MANT_DIG
#define IEEE_FLT_MANT_DIG 24
#endif
#define MPFR_LIMBS_PER_FLT ((IEEE_FLT_MANT_DIG-1)/GMP_NUMB_BITS+1)

/* Visual C++ doesn't support +1.0/0.0, -1.0/0.0 and 0.0/0.0
   at compile time.
   Clang with -fsanitize=undefined is a bit similar due to a bug:
     http://llvm.org/bugs/show_bug.cgi?id=17381
   but even without its sanitizer, it may be better to use the
   double_zero version until IEEE 754 division by zero is properly
   supported:
     http://llvm.org/bugs/show_bug.cgi?id=17000
*/
#if (defined(_MSC_VER) && defined(_WIN32) && (_MSC_VER >= 1200)) || \
    defined(__clang__)
static double double_zero = 0.0;
# define DBL_NAN (double_zero/double_zero)
# define DBL_POS_INF ((double) 1.0/double_zero)
# define DBL_NEG_INF ((double)-1.0/double_zero)
# define DBL_NEG_ZERO (-double_zero)
#else
# define DBL_POS_INF ((double) 1.0/0.0)
# define DBL_NEG_INF ((double)-1.0/0.0)
# define DBL_NAN     ((double) 0.0/0.0)
# define DBL_NEG_ZERO (-0.0)
#endif

/* Note: In the past, there was specific code for _GMP_IEEE_FLOATS, which
   was based on NaN and Inf memory representations. This code was breaking
   the aliasing rules (see ISO C99, 6.5#6 and 6.5#7 on the effective type)
   and for this reason it did not behave correctly with GCC 4.5.0 20091119.
   The code needed a memory transfer and was probably not better than the
   macros below with a good compiler (a fix based on the NaN / Inf memory
   representation would be even worse due to C limitations), and this code
   could be selected only when MPFR was built with --with-gmp-build, thus
   introducing a difference (bad for maintaining/testing MPFR); therefore
   it has been removed. The old code required that the argument x be an
   lvalue of type double. We still require that, in case one would need
   to change the macros below, e.g. for some broken compiler. But the
   LVALUE(x) condition could be removed if really necessary. */
/* Below, the &(x) == &(x) || &(x) != &(x) allows to make sure that x
   is a lvalue without (probably) any warning from the compiler.  The
   &(x) != &(x) is needed to avoid a failure under Mac OS X 10.4.11
   (with Xcode 2.4.1, i.e. the latest one). */
#define LVALUE(x) (&(x) == &(x) || &(x) != &(x))
#define DOUBLE_ISINF(x) (LVALUE(x) && ((x) > DBL_MAX || (x) < -DBL_MAX))
/* The DOUBLE_ISNAN(x) macro is also valid on long double x
   (assuming that the compiler isn't too broken). */
#ifdef MPFR_NANISNAN
/* Avoid MIPSpro / IRIX64 / gcc -ffast-math (incorrect) optimizations.
   The + must not be replaced by a ||. With gcc -ffast-math, NaN is
   regarded as a positive number or something like that; the second
   test catches this case. */
# define DOUBLE_ISNAN(x) \
    (LVALUE(x) && !((((x) >= 0.0) + ((x) <= 0.0)) && -(x)*(x) <= 0.0))
#else
# define DOUBLE_ISNAN(x) (LVALUE(x) && (x) != (x))
#endif

/******************************************************
 *************** Long double macros *******************
 ******************************************************/

/* We try to get the exact value of the precision of long double
   (provided by the implementation) in order to provide correct
   rounding in this case (not guaranteed if the C implementation
   does not have an adequate long double arithmetic). Note that
   it may be lower than the precision of some numbers that can
   be represented in a long double; e.g. on FreeBSD/x86, it is
   53 because the processor is configured to round in double
   precision (even when using the long double type -- this is a
   limitation of the x87 arithmetic), and on Mac OS X, it is 106
   because the implementation is a double-double arithmetic.
   Otherwise (e.g. in base 10), we get an upper bound of the
   precision, and correct rounding isn't currently provided.
*/
#if defined(LDBL_MANT_DIG) && FLT_RADIX == 2
# define MPFR_LDBL_MANT_DIG LDBL_MANT_DIG
#else
# define MPFR_LDBL_MANT_DIG \
  (sizeof(long double)*GMP_NUMB_BITS/sizeof(mp_limb_t))
#endif
#define MPFR_LIMBS_PER_LONG_DOUBLE \
  ((sizeof(long double)-1)/sizeof(mp_limb_t)+1)

/* this is standardized by IEEE 754 */
#define IEEE_FLOAT128_MANT_DIG 113

/* LONGDOUBLE_NAN_ACTION executes the code "action" if x is a NaN. */

/* On hppa2.0n-hp-hpux10 with the unbundled HP cc, the test x!=x on a NaN
   has been seen false, meaning NaNs are not detected.  This seemed to
   happen only after other comparisons, not sure what's really going on.  In
   any case we can pick apart the bytes to identify a NaN.  */
#ifdef HAVE_LDOUBLE_IEEE_QUAD_BIG
# define LONGDOUBLE_NAN_ACTION(x, action)                       \
  do {                                                          \
    union {                                                     \
      long double    ld;                                        \
      struct {                                                  \
        unsigned int sign : 1;                                  \
        unsigned int exp  : 15;                                 \
        unsigned int man3 : 16;                                 \
        unsigned int man2 : 32;                                 \
        unsigned int man1 : 32;                                 \
        unsigned int man0 : 32;                                 \
      } s;                                                      \
    } u;                                                        \
    u.ld = (x);                                                 \
    if (u.s.exp == 0x7FFFL                                      \
        && (u.s.man0 | u.s.man1 | u.s.man2 | u.s.man3) != 0)    \
      { action; }                                               \
  } while (0)
#endif

#ifdef HAVE_LDOUBLE_IEEE_QUAD_LITTLE
# define LONGDOUBLE_NAN_ACTION(x, action)                       \
  do {                                                          \
    union {                                                     \
      long double    ld;                                        \
      struct {                                                  \
        unsigned int man0 : 32;                                 \
        unsigned int man1 : 32;                                 \
        unsigned int man2 : 32;                                 \
        unsigned int man3 : 16;                                 \
        unsigned int exp  : 15;                                 \
        unsigned int sign : 1;                                  \
      } s;                                                      \
    } u;                                                        \
    u.ld = (x);                                                 \
    if (u.s.exp == 0x7FFFL                                      \
        && (u.s.man0 | u.s.man1 | u.s.man2 | u.s.man3) != 0)    \
      { action; }                                               \
  } while (0)
#endif

/* Under IEEE rules, NaN is not equal to anything, including itself.
   "volatile" here stops "cc" on mips64-sgi-irix6.5 from optimizing away
   x!=x. */
#ifndef LONGDOUBLE_NAN_ACTION
# define LONGDOUBLE_NAN_ACTION(x, action)               \
  do {                                                  \
    volatile long double __x = LONGDOUBLE_VOLATILE (x); \
    if ((x) != __x)                                     \
      { action; }                                       \
  } while (0)
# define WANT_LONGDOUBLE_VOLATILE 1
#endif

/* If we don't have a proper "volatile" then volatile is #defined to empty,
   in this case call through an external function to stop the compiler
   optimizing anything. */
#ifdef WANT_LONGDOUBLE_VOLATILE
# ifdef volatile
__MPFR_DECLSPEC long double __gmpfr_longdouble_volatile _MPFR_PROTO ((long double)) MPFR_CONST_ATTR;
#  define LONGDOUBLE_VOLATILE(x)  (__gmpfr_longdouble_volatile (x))
#  define WANT_GMPFR_LONGDOUBLE_VOLATILE 1
# else
#  define LONGDOUBLE_VOLATILE(x)  (x)
# endif
#endif

/* Some special case for IEEE_EXT Litle Endian */
#if HAVE_LDOUBLE_IEEE_EXT_LITTLE

typedef union {
  long double    ld;
  struct {
    unsigned int manl : 32;
    unsigned int manh : 32;
    unsigned int expl : 8 ;
    unsigned int exph : 7;
    unsigned int sign : 1;
  } s;
} mpfr_long_double_t;

/* #undef MPFR_LDBL_MANT_DIG */
#undef MPFR_LIMBS_PER_LONG_DOUBLE
/* #define MPFR_LDBL_MANT_DIG   64 */
#define MPFR_LIMBS_PER_LONG_DOUBLE ((64-1)/GMP_NUMB_BITS+1)

#endif

/******************************************************
 *************** _Decimal64 support *******************
 ******************************************************/

#ifdef MPFR_WANT_DECIMAL_FLOATS
/* to cast between binary64 and decimal64 */
union ieee_double_decimal64 { double d; _Decimal64 d64; };
#endif /* MPFR_WANT_DECIMAL_FLOATS */

/******************************************************
 **************** mpfr_t properties *******************
 ******************************************************/

#define MPFR_PREC_COND(p) ((p) >= MPFR_PREC_MIN && (p) <= MPFR_PREC_MAX)
#define MPFR_PREC_IN_RANGE(p) (MPFR_ASSERTD (MPFR_PREC_COND(p)), (p))

/* In the following macro, p is usually a mpfr_prec_t, but this macro
   works with other integer types (without integer overflow). Checking
   that p >= 1 in debug mode is useful here because this macro can be
   used on a computed precision (in particular, this formula does not
   work for a degenerate case p = 0, and could give different results
   on different platforms). But let us not use an assertion checking
   in the MPFR_LAST_LIMB() and MPFR_LIMB_SIZE() macros below to avoid
   too much expansion for assertions (in practice, this should be a
   problem just when testing MPFR with the --enable-assert configure
   option and the -ansi -pedantic-errors gcc compiler flags). */
#define MPFR_PREC2LIMBS(p) \
  (MPFR_ASSERTD ((p) >= 1), ((p) - 1) / GMP_NUMB_BITS + 1)

#define MPFR_PREC(x)      ((x)->_mpfr_prec)
#define MPFR_EXP(x)       ((x)->_mpfr_exp)
#define MPFR_MANT(x)      ((x)->_mpfr_d)
#define MPFR_GET_PREC(x)  MPFR_PREC_IN_RANGE (MPFR_PREC (x))
#define MPFR_LAST_LIMB(x) ((MPFR_GET_PREC (x) - 1) / GMP_NUMB_BITS)
#define MPFR_LIMB_SIZE(x) (MPFR_LAST_LIMB (x) + 1)


/******************************************************
 **************** exponent properties *****************
 ******************************************************/

/* Limits of the mpfr_exp_t type (NOT those of valid exponent values).
   These macros can be used in preprocessor directives. */
#if   _MPFR_EXP_FORMAT == 1
# define MPFR_EXP_MAX (SHRT_MAX)
# define MPFR_EXP_MIN (SHRT_MIN)
#elif _MPFR_EXP_FORMAT == 2
# define MPFR_EXP_MAX (INT_MAX)
# define MPFR_EXP_MIN (INT_MIN)
#elif _MPFR_EXP_FORMAT == 3
# define MPFR_EXP_MAX (LONG_MAX)
# define MPFR_EXP_MIN (LONG_MIN)
#elif _MPFR_EXP_FORMAT == 4
/* Note: MPFR_EXP_MAX and MPFR_EXP_MIN must not be used in #if directives
   if _MPFR_EXP_FORMAT == 4 and MPFR_HAVE_INTMAX_MAX is not defined. */
# define MPFR_EXP_MAX (MPFR_INTMAX_MAX)
# define MPFR_EXP_MIN (MPFR_INTMAX_MIN)
#else
# error "Invalid MPFR Exp format"
#endif

/* Before doing a cast to mpfr_uexp_t, make sure that the value is
   nonnegative. */
#define MPFR_UEXP(X) (MPFR_ASSERTD ((X) >= 0), (mpfr_uexp_t) (X))

#if _MPFR_EXP_FORMAT <= 3
typedef long int mpfr_eexp_t;
# define mpfr_get_exp_t(x,r) mpfr_get_si((x),(r))
# define mpfr_set_exp_t(x,e,r) mpfr_set_si((x),(e),(r))
# define MPFR_EXP_FSPEC "l"
#else
typedef intmax_t mpfr_eexp_t;
# define mpfr_get_exp_t(x,r) mpfr_get_sj((x),(r))
# define mpfr_set_exp_t(x,e,r) mpfr_set_sj((x),(e),(r))
# define MPFR_EXP_FSPEC "j"
#endif

/* Invalid exponent value (to track bugs...) */
#define MPFR_EXP_INVALID \
 ((mpfr_exp_t) 1 << (GMP_NUMB_BITS*sizeof(mpfr_exp_t)/sizeof(mp_limb_t)-2))

/* Definition of the exponent limits for MPFR numbers.
 * These limits are chosen so that if e is such an exponent, then 2e-1 and
 * 2e+1 are representable. This is useful for intermediate computations,
 * in particular the multiplication.
 */
#undef MPFR_EMIN_MIN
#undef MPFR_EMIN_MAX
#undef MPFR_EMAX_MIN
#undef MPFR_EMAX_MAX
#define MPFR_EMIN_MIN (1-MPFR_EXP_INVALID)
#define MPFR_EMIN_MAX (MPFR_EXP_INVALID-1)
#define MPFR_EMAX_MIN (1-MPFR_EXP_INVALID)
#define MPFR_EMAX_MAX (MPFR_EXP_INVALID-1)

/* Use MPFR_GET_EXP and MPFR_SET_EXP instead of MPFR_EXP directly,
   unless when the exponent may be out-of-range, for instance when
   setting the exponent before calling mpfr_check_range.
   MPFR_EXP_CHECK is defined when MPFR_WANT_ASSERT is defined, but if you
   don't use MPFR_WANT_ASSERT (for speed reasons), you can still define
   MPFR_EXP_CHECK by setting -DMPFR_EXP_CHECK in $CFLAGS.
   Note about MPFR_SET_EXP:
     The exp expression is required to have a signed type. To avoid spurious
     failures, we could cast (exp) to mpfr_exp_t, but this wouldn't allow us
     to detect some bugs that can occur on particular platforms. Anyway, an
     unsigned type for exp is suspicious and should be regarded as a bug.
*/

#ifdef MPFR_EXP_CHECK
# define MPFR_GET_EXP(x)          (mpfr_get_exp) (x)
# define MPFR_SET_EXP(x, exp)     (MPFR_ASSERTN ((exp) >= __gmpfr_emin && \
                                                 (exp) <= __gmpfr_emax),  \
                                   (void) (MPFR_EXP (x) = (exp)))
# define MPFR_SET_INVALID_EXP(x)  ((void) (MPFR_EXP (x) = MPFR_EXP_INVALID))
#else
# define MPFR_GET_EXP(x)          MPFR_EXP (x)
# define MPFR_SET_EXP(x, exp)     ((void) (MPFR_EXP (x) = (exp)))
# define MPFR_SET_INVALID_EXP(x)  ((void) 0)
#endif



/******************************************************
 ********** Singular Values (NAN, INF, ZERO) **********
 ******************************************************/

/* Enum special value of exponent. */
# define MPFR_EXP_ZERO (MPFR_EXP_MIN+1)
# define MPFR_EXP_NAN  (MPFR_EXP_MIN+2)
# define MPFR_EXP_INF  (MPFR_EXP_MIN+3)

#define MPFR_IS_NAN(x)   (MPFR_EXP(x) == MPFR_EXP_NAN)
#define MPFR_SET_NAN(x)  (MPFR_EXP(x) =  MPFR_EXP_NAN)
#define MPFR_IS_INF(x)   (MPFR_EXP(x) == MPFR_EXP_INF)
#define MPFR_SET_INF(x)  (MPFR_EXP(x) =  MPFR_EXP_INF)
#define MPFR_IS_ZERO(x)  (MPFR_EXP(x) == MPFR_EXP_ZERO)
#define MPFR_SET_ZERO(x) (MPFR_EXP(x) =  MPFR_EXP_ZERO)
#define MPFR_NOTZERO(x)  (MPFR_EXP(x) != MPFR_EXP_ZERO)

#define MPFR_IS_NORMALIZED(x) \
  (MPFR_LIMB_MSB (MPFR_MANT(x)[MPFR_LAST_LIMB(x)]) != 0)

#define MPFR_IS_FP(x)       (!MPFR_IS_NAN(x) && !MPFR_IS_INF(x))
#define MPFR_IS_SINGULAR(x) (MPFR_EXP(x) <= MPFR_EXP_INF)
#define MPFR_IS_PURE_FP(x) \
  (!MPFR_IS_SINGULAR(x) && (MPFR_ASSERTD (MPFR_IS_NORMALIZED (x)), 1))

#define MPFR_ARE_SINGULAR(x,y) \
  (MPFR_UNLIKELY(MPFR_IS_SINGULAR(x)) || MPFR_UNLIKELY(MPFR_IS_SINGULAR(y)))

#define MPFR_IS_POWER_OF_2(x) \
  (mpfr_cmp_ui_2exp ((x), 1, MPFR_GET_EXP (x) - 1) == 0)


/******************************************************
 ********************* Sign Macros ********************
 ******************************************************/

#define MPFR_SIGN_POS (1)
#define MPFR_SIGN_NEG (-1)

#define MPFR_IS_STRICTPOS(x)  (MPFR_NOTZERO (x) && MPFR_IS_POS (x))
#define MPFR_IS_STRICTNEG(x)  (MPFR_NOTZERO (x) && MPFR_IS_NEG (x))

#define MPFR_IS_NEG(x) (MPFR_SIGN(x) < 0)
#define MPFR_IS_POS(x) (MPFR_SIGN(x) > 0)

#define MPFR_SET_POS(x) (MPFR_SIGN(x) = MPFR_SIGN_POS)
#define MPFR_SET_NEG(x) (MPFR_SIGN(x) = MPFR_SIGN_NEG)

#define MPFR_CHANGE_SIGN(x) (MPFR_SIGN(x) = -MPFR_SIGN(x))
#define MPFR_SET_SAME_SIGN(x, y) (MPFR_SIGN(x) = MPFR_SIGN(y))
#define MPFR_SET_OPPOSITE_SIGN(x, y) (MPFR_SIGN(x) = -MPFR_SIGN(y))
#define MPFR_ASSERT_SIGN(s) \
 (MPFR_ASSERTD((s) == MPFR_SIGN_POS || (s) == MPFR_SIGN_NEG))
#define MPFR_SET_SIGN(x, s) \
  (MPFR_ASSERT_SIGN(s), MPFR_SIGN(x) = s)
#define MPFR_IS_POS_SIGN(s1) ((s1) > 0)
#define MPFR_IS_NEG_SIGN(s1) ((s1) < 0)
#define MPFR_MULT_SIGN(s1, s2) ((s1) * (s2))
/* Transform a sign to 1 or -1 */
#define MPFR_FROM_SIGN_TO_INT(s) (s)
#define MPFR_INT_SIGN(x) MPFR_FROM_SIGN_TO_INT(MPFR_SIGN(x))



/******************************************************
 ***************** Ternary Value Macros ***************
 ******************************************************/

/* Special inexact value */
#define MPFR_EVEN_INEX 2

/* Macros for functions returning two inexact values in an 'int' */
#define INEXPOS(y) ((y) == 0 ? 0 : (((y) > 0) ? 1 : 2))
#define INEX(y,z) (INEXPOS(y) | (INEXPOS(z) << 2))

/* When returning the ternary inexact value, ALWAYS use one of the
   following two macros, unless the flag comes from another function
   returning the ternary inexact value */
#define MPFR_RET(I) return \
  (I) ? ((__gmpfr_flags |= MPFR_FLAGS_INEXACT), (I)) : 0
#define MPFR_RET_NAN return (__gmpfr_flags |= MPFR_FLAGS_NAN), 0

#define SIGN(I) ((I) < 0 ? -1 : (I) > 0)
#define SAME_SIGN(I1,I2) (SIGN (I1) == SIGN (I2))



/******************************************************
 ************** Rounding mode macros  *****************
 ******************************************************/

/* MPFR_RND_MAX gives the number of supported rounding modes by all functions.
 * Once faithful rounding is implemented, MPFR_RNDA should be changed
 * to MPFR_RNDF. But this will also require more changes in the tests.
 */
#define MPFR_RND_MAX ((mpfr_rnd_t)((MPFR_RNDA)+1))

/* We want to test this :
 *  (rnd == MPFR_RNDU && test) || (rnd == RNDD && !test)
 * ie it transforms RNDU or RNDD to Away or Zero according to the sign */
#define MPFR_IS_RNDUTEST_OR_RNDDNOTTEST(rnd, test) \
  (((rnd) + (test)) == MPFR_RNDD)

/* We want to test if rnd = Zero, or Away.
   'test' is 1 if negative, and 0 if positive. */
#define MPFR_IS_LIKE_RNDZ(rnd, test) \
  ((rnd) == MPFR_RNDZ || MPFR_IS_RNDUTEST_OR_RNDDNOTTEST (rnd, test))

#define MPFR_IS_LIKE_RNDU(rnd, sign)                    \
  (((rnd) == MPFR_RNDU) ||                              \
   ((rnd) == MPFR_RNDZ && MPFR_IS_NEG_SIGN (sign)) ||   \
   ((rnd) == MPFR_RNDA && MPFR_IS_POS_SIGN (sign)))

#define MPFR_IS_LIKE_RNDD(rnd, sign)                    \
  (((rnd) == MPFR_RNDD) ||                              \
   ((rnd) == MPFR_RNDZ && MPFR_IS_POS_SIGN (sign)) ||   \
   ((rnd) == MPFR_RNDA && MPFR_IS_NEG_SIGN (sign)))

/* Invert a rounding mode, RNDZ and RNDA are unchanged */
#define MPFR_INVERT_RND(rnd) ((rnd) == MPFR_RNDU ? MPFR_RNDD :          \
                              (rnd) == MPFR_RNDD ? MPFR_RNDU : (rnd))

/* Transform RNDU and RNDD to RNDZ according to test */
#define MPFR_UPDATE_RND_MODE(rnd, test)                             \
  do {                                                              \
    if (MPFR_UNLIKELY(MPFR_IS_RNDUTEST_OR_RNDDNOTTEST(rnd, test)))  \
      rnd = MPFR_RNDZ;                                              \
  } while (0)

/* Transform RNDU and RNDD to RNDZ or RNDA according to sign,
   leave the other modes unchanged */
#define MPFR_UPDATE2_RND_MODE(rnd, sign)                       \
  do {                                                         \
    if (rnd == MPFR_RNDU)                                      \
      rnd = MPFR_IS_POS_SIGN (sign) ? MPFR_RNDA : MPFR_RNDZ;   \
    else if (rnd == MPFR_RNDD)                                 \
      rnd = MPFR_IS_NEG_SIGN (sign) ? MPFR_RNDA : MPFR_RNDZ;   \
  } while (0)


/******************************************************
 ******************* Limb Macros **********************
 ******************************************************/

/* Definition of simple mp_limb_t constants */
#define MPFR_LIMB_ZERO    ((mp_limb_t) 0)
#define MPFR_LIMB_ONE     ((mp_limb_t) 1)
#define MPFR_LIMB_HIGHBIT (MPFR_LIMB_ONE << (GMP_NUMB_BITS - 1))
#define MPFR_LIMB_MAX     ((mp_limb_t) -1)

/* Mask to get the Most Significant Bit of a limb */
#define MPFR_LIMB_MSB(l) ((l) & MPFR_LIMB_HIGHBIT)

/* Mask for the low 's' bits of a limb */
#define MPFR_LIMB_MASK(s) ((MPFR_LIMB_ONE << (s)) - MPFR_LIMB_ONE)



/******************************************************
 ********************** Memory ************************
 ******************************************************/

/* Heap memory handling */
typedef union { mp_size_t s; mp_limb_t l; } mpfr_size_limb_t;
#define MPFR_GET_ALLOC_SIZE(x) \
 ( ((mp_size_t*) MPFR_MANT(x))[-1] + 0)
#define MPFR_SET_ALLOC_SIZE(x, n) \
 ( ((mp_size_t*) MPFR_MANT(x))[-1] = n)
#define MPFR_MALLOC_SIZE(s) \
  ( sizeof(mpfr_size_limb_t) + MPFR_BYTES_PER_MP_LIMB * ((size_t) s) )
#define MPFR_SET_MANT_PTR(x,p) \
   (MPFR_MANT(x) = (mp_limb_t*) ((mpfr_size_limb_t*) p + 1))
#define MPFR_GET_REAL_PTR(x) \
   ((mp_limb_t*) ((mpfr_size_limb_t*) MPFR_MANT(x) - 1))

/* Temporary memory handling */
#ifndef TMP_SALLOC
/* GMP 4.1.x or below or internals */
#define MPFR_TMP_DECL TMP_DECL
#define MPFR_TMP_MARK TMP_MARK
#define MPFR_TMP_ALLOC TMP_ALLOC
#define MPFR_TMP_FREE TMP_FREE
#else
#define MPFR_TMP_DECL(x) TMP_DECL
#define MPFR_TMP_MARK(x) TMP_MARK
#define MPFR_TMP_ALLOC(s) TMP_ALLOC(s)
#define MPFR_TMP_FREE(x) TMP_FREE
#endif

#define MPFR_TMP_LIMBS_ALLOC(N) \
  ((mp_limb_t *) MPFR_TMP_ALLOC ((size_t) (N) * MPFR_BYTES_PER_MP_LIMB))

/* temporary allocate 1 limb at xp, and initialize mpfr variable x */
/* The temporary var doesn't have any size field, but it doesn't matter
 * since only functions dealing with the Heap care about it */
#define MPFR_TMP_INIT1(xp, x, p)                                     \
 ( MPFR_PREC(x) = (p),                                               \
   MPFR_MANT(x) = (xp),                                              \
   MPFR_SET_POS(x),                                                  \
   MPFR_SET_INVALID_EXP(x))

#define MPFR_TMP_INIT(xp, x, p, s)                                   \
  (xp = MPFR_TMP_LIMBS_ALLOC(s),                                     \
   MPFR_TMP_INIT1(xp, x, p))

#define MPFR_TMP_INIT_ABS(d, s)                                      \
 ( MPFR_PREC(d) = MPFR_PREC(s),                                      \
   MPFR_MANT(d) = MPFR_MANT(s),                                      \
   MPFR_SET_POS(d),                                                  \
   MPFR_EXP(d)  = MPFR_EXP(s))



/******************************************************
 *****************  Cache macros **********************
 ******************************************************/

#define mpfr_const_pi(_d,_r)    mpfr_cache(_d, __gmpfr_cache_const_pi,_r)
#define mpfr_const_log2(_d,_r)  mpfr_cache(_d, __gmpfr_cache_const_log2, _r)
#define mpfr_const_euler(_d,_r) mpfr_cache(_d, __gmpfr_cache_const_euler, _r)
#define mpfr_const_catalan(_d,_r) mpfr_cache(_d,__gmpfr_cache_const_catalan,_r)

#define MPFR_DECL_INIT_CACHE(_cache,_func)                           \
  MPFR_THREAD_ATTR mpfr_cache_t _cache =                             \
    {{{{0,MPFR_SIGN_POS,0,(mp_limb_t*)0}},0,_func}}



/******************************************************
 *******************  Threshold ***********************
 ******************************************************/

#include "mparam.h"

/******************************************************
 *****************  Useful macros *********************
 ******************************************************/

/* Theses macros help the compiler to determine if a test is
   likely or unlikely. The !! is necessary in case x is larger
   than a long. */
#if defined MPFR_DEBUG_PREDICTION && __MPFR_GNUC(3,0)

/* Code to debug branch prediction, based on Ulrich Drepper's paper
 * "What Every Programmer Should Know About Memory":
 *   http://people.freebsd.org/~lstewart/articles/cpumemory.pdf
 */
asm (".section predict_data, \"aw\"; .previous\n"
     ".section predict_line, \"a\"; .previous\n"
     ".section predict_file, \"a\"; .previous");
# if defined __x86_64__
#  define MPFR_DEBUGPRED__(e,E)                                         \
  ({ long int _e = !!(e);                                               \
    asm volatile (".pushsection predict_data\n"                         \
                  "..predictcnt%=: .quad 0; .quad 0\n"                  \
                  ".section predict_line; .quad %c1\n"                  \
                  ".section predict_file; .quad %c2; .popsection\n"     \
                  "addq $1,..predictcnt%=(,%0,8)"                       \
                  : : "r" (_e == E), "i" (__LINE__), "i" (__FILE__));   \
    __builtin_expect (_e, E);                                           \
  })
# elif defined __i386__
#  define MPFR_DEBUGPRED__(e,E)                                         \
  ({ long int _e = !!(e);                                               \
    asm volatile (".pushsection predict_data\n"                         \
                  "..predictcnt%=: .long 0; .long 0\n"                  \
                  ".section predict_line; .long %c1\n"                  \
                  ".section predict_file; .long %c2; .popsection\n"     \
                  "incl ..predictcnt%=(,%0,4)"                          \
                  : : "r" (_e == E), "i" (__LINE__), "i" (__FILE__));   \
    __builtin_expect (_e, E);                                           \
  })
# else
#  error "MPFR_DEBUGPRED__ definition missing"
# endif

# define MPFR_LIKELY(x) MPFR_DEBUGPRED__ ((x), 1)
# define MPFR_UNLIKELY(x) MPFR_DEBUGPRED__ ((x), 0)

#elif __MPFR_GNUC(3,0) || __MPFR_ICC(8,1,0)

# define MPFR_LIKELY(x) (__builtin_expect(!!(x), 1))
# define MPFR_UNLIKELY(x) (__builtin_expect(!!(x), 0))

#else

# define MPFR_LIKELY(x) (x)
# define MPFR_UNLIKELY(x) (x)

#endif

/* Declare that some variable is initialized before being used (without a
   dummy initialization) in order to avoid some compiler warnings. Use the
   VAR = VAR trick (see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=36296#c3)
   only with gcc as this is undefined behavior, and we don't know what other
   compilers do (they may also be smarter). This self-initialization trick
   could be disabled with future gcc versions.
   However, for clang (which defines __GNUC__), this trick must not be used
   as it currently generates a warning, at least with:
     Debian clang version 3.0-6.2 (tags/RELEASE_30/final) (based on LLVM 3.0)
     __VERSION__ "4.2.1 Compatible Debian Clang 3.0 (tags/RELEASE_30/final)"
     __clang__ 1
     __clang_major__ 3
     __clang_minor__ 0
     __clang_patchlevel__ 0
     __clang_version__ "3.0 (tags/RELEASE_30/final)"
   (see https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=705583 for this
   problem with clang). */
#if defined(__GNUC__) && !defined(__clang__)
# define INITIALIZED(VAR) VAR = VAR
#else
# define INITIALIZED(VAR) VAR
#endif

/* Ceil log 2: If GCC, uses a GCC extension, otherwise calls a function */
/* Warning:
 *   Needs to define MPFR_NEED_LONGLONG.
 *   Computes ceil(log2(x)) only for x integer (unsigned long)
 *   Undefined if x is 0 */
#if __MPFR_GNUC(2,95) || __MPFR_ICC(8,1,0)
/* Note: This macro MPFR_INT_CEIL_LOG2 shouldn't be used in an MPFR_ASSERT*
   macro, either directly or indirectly via other macros, otherwise it can
   yield an error due to a too large stringized expression in ASSERT_FAIL.
   A static inline function could be a better solution than this macro. */
# define MPFR_INT_CEIL_LOG2(x)                            \
    (MPFR_UNLIKELY ((x) == 1) ? 0 :                       \
     __extension__ ({ int _b; mp_limb_t _limb;            \
      MPFR_ASSERTN ((x) > 1);                             \
      _limb = (x) - 1;                                    \
      MPFR_ASSERTN (_limb == (x) - 1);                    \
      count_leading_zeros (_b, _limb);                    \
      (GMP_NUMB_BITS - _b); }))
#else
# define MPFR_INT_CEIL_LOG2(x) (__gmpfr_int_ceil_log2(x))
#endif

/* Add two integers with overflow handling */
/* Example: MPFR_SADD_OVERFLOW (c, a, b, long, unsigned long,
 *                              LONG_MIN, LONG_MAX,
 *                              goto overflow, goto underflow); */
#define MPFR_UADD_OVERFLOW(c,a,b,ACTION_IF_OVERFLOW)                  \
 do {                                                                 \
  (c) = (a) + (b);                                                    \
  if ((c) < (a)) ACTION_IF_OVERFLOW;                                  \
 } while (0)

#define MPFR_SADD_OVERFLOW(c,a,b,STYPE,UTYPE,MIN,MAX,ACTION_IF_POS_OVERFLOW,ACTION_IF_NEG_OVERFLOW) \
  do {                                                                \
  if ((a) >= 0 && (b) >= 0) {                                         \
         UTYPE uc,ua,ub;                                              \
         ua = (UTYPE) (a); ub = (UTYPE) (b);                          \
         MPFR_UADD_OVERFLOW (uc, ua, ub, ACTION_IF_POS_OVERFLOW);     \
         if (uc > (UTYPE)(MAX)) ACTION_IF_POS_OVERFLOW;               \
         else (c) = (STYPE) uc;                                       \
  } else if ((a) < 0 && (b) < 0) {                                    \
         UTYPE uc,ua,ub;                                              \
         ua = -(UTYPE) (a); ub = -(UTYPE) (b);                        \
         MPFR_UADD_OVERFLOW (uc, ua, ub, ACTION_IF_NEG_OVERFLOW);     \
         if (uc >= -(UTYPE)(MIN) || uc > (UTYPE)(MAX)) {              \
           if (uc ==  -(UTYPE)(MIN)) (c) = (MIN);                     \
           else  ACTION_IF_NEG_OVERFLOW; }                            \
         else (c) = -(STYPE) uc;                                      \
  } else (c) = (a) + (b);                                             \
 } while (0)


/* Set a number to 1 (Fast) - It doesn't check if 1 is in the exponent range */
#define MPFR_SET_ONE(x)                                               \
do {                                                                  \
  mp_size_t _size = MPFR_LAST_LIMB(x);                                \
  MPFR_SET_POS(x);                                                    \
  MPFR_EXP(x) = 1;                                                    \
  MPN_ZERO ( MPFR_MANT(x), _size);                                    \
  MPFR_MANT(x)[_size] = MPFR_LIMB_HIGHBIT;                            \
} while (0)

/* Compute s = (-a) % GMP_NUMB_BITS as unsigned */
#define MPFR_UNSIGNED_MINUS_MODULO(s, a)                              \
  do                                                                  \
    {                                                                 \
      if (IS_POW2 (GMP_NUMB_BITS))                                    \
        (s) = (- (unsigned int) (a)) % GMP_NUMB_BITS;                 \
      else                                                            \
        {                                                             \
          (s) = (a) % GMP_NUMB_BITS;                                  \
          if ((s) != 0)                                               \
            (s) = GMP_NUMB_BITS - (s);                                \
        }                                                             \
      MPFR_ASSERTD ((s) >= 0 && (s) < GMP_NUMB_BITS);                 \
    }                                                                 \
  while (0)

/* Use it only for debug reasons */
/*   MPFR_TRACE (operation) : execute operation iff DEBUG flag is set */
/*   MPFR_DUMP (x) : print x (a mpfr_t) on stdout */
#ifdef DEBUG
# define MPFR_TRACE(x) x
#else
# define MPFR_TRACE(x) (void) 0
#endif
#define MPFR_DUMP(x) ( printf(#x"="), mpfr_dump(x) )

/* Test if X (positive) is a power of 2 */
#define IS_POW2(X) (((X) & ((X) - 1)) == 0)
#define NOT_POW2(X) (((X) & ((X) - 1)) != 0)

/* Safe absolute value and difference (to avoid possible integer overflow) */
/* type is the target (unsigned) type */
#define SAFE_ABS(type,x) ((x) >= 0 ? (type)(x) : -(type)(x))
#define SAFE_DIFF(type,x,y) (MPFR_ASSERTD((x) >= (y)), (type)(x) - (type)(y))

#define mpfr_get_d1(x) mpfr_get_d(x,__gmpfr_default_rounding_mode)

/* Store in r the size in bits of the mpz_t z */
#define MPFR_MPZ_SIZEINBASE2(r, z)              \
  do {                                          \
   int _cnt;                                    \
   mp_size_t _size;                             \
   MPFR_ASSERTD (mpz_sgn (z) != 0);             \
   _size = ABSIZ(z);                            \
   MPFR_ASSERTD (_size >= 1);                   \
   count_leading_zeros (_cnt, PTR(z)[_size-1]); \
   (r) = _size * GMP_NUMB_BITS - _cnt;          \
  } while (0)

/* MPFR_LCONV_DPTS can also be forced to 0 or 1 by the user. */
#ifndef MPFR_LCONV_DPTS
# if defined(HAVE_LOCALE_H) && \
     defined(HAVE_STRUCT_LCONV_DECIMAL_POINT) && \
     defined(HAVE_STRUCT_LCONV_THOUSANDS_SEP)
#  define MPFR_LCONV_DPTS 1
# else
#  define MPFR_LCONV_DPTS 0
# endif
#endif

#if MPFR_LCONV_DPTS
#include <locale.h>
/* Warning! In case of signed char, the value of MPFR_DECIMAL_POINT may
   be negative (the ISO C99 does not seem to forbid negative values). */
#define MPFR_DECIMAL_POINT (localeconv()->decimal_point[0])
#define MPFR_THOUSANDS_SEPARATOR (localeconv()->thousands_sep[0])
#else
#define MPFR_DECIMAL_POINT ((char) '.')
#define MPFR_THOUSANDS_SEPARATOR ('\0')
#endif


/* Set y to s*significand(x)*2^e, for example MPFR_ALIAS(y,x,1,MPFR_EXP(x))
   sets y to |x|, and MPFR_ALIAS(y,x,MPFR_SIGN(x),0) sets y to x*2^f such
   that 1/2 <= |y| < 1. Does not check y is in the valid exponent range.
   WARNING! x and y share the same mantissa. So, some operations are
   not valid if x has been provided via an argument, e.g., trying to
   modify the mantissa of y, even temporarily, or calling mpfr_clear on y.
*/
#define MPFR_ALIAS(y,x,s,e)                     \
  do                                            \
    {                                           \
      MPFR_PREC(y) = MPFR_PREC(x);              \
      MPFR_SIGN(y) = (s);                       \
      MPFR_EXP(y) = (e);                        \
      MPFR_MANT(y) = MPFR_MANT(x);              \
    } while (0)


/******************************************************
 **************  Save exponent macros  ****************
 ******************************************************/

/* See README.dev for details on how to use the macros.
   They are used to set the exponent range to the maximum
   temporarily */

typedef struct {
  mpfr_flags_t saved_flags;
  mpfr_exp_t saved_emin;
  mpfr_exp_t saved_emax;
} mpfr_save_expo_t;

/* Minimum and maximum exponents of the extended exponent range. */
#define MPFR_EXT_EMIN MPFR_EMIN_MIN
#define MPFR_EXT_EMAX MPFR_EMAX_MAX

#define MPFR_SAVE_EXPO_DECL(x) mpfr_save_expo_t x
#define MPFR_SAVE_EXPO_MARK(x)     \
 ((x).saved_flags = __gmpfr_flags, \
  (x).saved_emin = __gmpfr_emin,   \
  (x).saved_emax = __gmpfr_emax,   \
  __gmpfr_emin = MPFR_EXT_EMIN,    \
  __gmpfr_emax = MPFR_EXT_EMAX)
#define MPFR_SAVE_EXPO_FREE(x)     \
 (__gmpfr_flags = (x).saved_flags, \
  __gmpfr_emin = (x).saved_emin,   \
  __gmpfr_emax = (x).saved_emax)
#define MPFR_SAVE_EXPO_UPDATE_FLAGS(x, flags)  \
  (x).saved_flags |= (flags)

/* Speed up final checking */
#define mpfr_check_range(x,t,r) \
 (MPFR_LIKELY (MPFR_EXP (x) >= __gmpfr_emin && MPFR_EXP (x) <= __gmpfr_emax) \
  ? ((t) ? (__gmpfr_flags |= MPFR_FLAGS_INEXACT, (t)) : 0)                   \
  : mpfr_check_range(x,t,r))


/******************************************************
 *****************  Inline Rounding *******************
 ******************************************************/

/*
 * Note: due to the labels, one cannot use a macro MPFR_RNDRAW* more than
 * once in a function (otherwise these labels would not be unique).
 */

/*
 * Round mantissa (srcp, sprec) to mpfr_t dest using rounding mode rnd
 * assuming dest's sign is sign.
 * In rounding to nearest mode, execute MIDDLE_HANDLER when the value
 * is the middle of two consecutive numbers in dest precision.
 * Execute OVERFLOW_HANDLER in case of overflow when rounding.
 */
#define MPFR_RNDRAW_GEN(inexact, dest, srcp, sprec, rnd, sign,              \
                        MIDDLE_HANDLER, OVERFLOW_HANDLER)                   \
  do {                                                                      \
    mp_size_t _dests, _srcs;                                                \
    mp_limb_t *_destp;                                                      \
    mpfr_prec_t _destprec, _srcprec;                                        \
                                                                            \
    /* Check Trivial Case when Dest Mantissa has more bits than source */   \
    _srcprec = (sprec);                                                     \
    _destprec = MPFR_PREC (dest);                                           \
    MPFR_ASSERTD (_srcprec >= MPFR_PREC_MIN);                               \
    MPFR_ASSERTD (_destprec >= MPFR_PREC_MIN);                              \
    _destp = MPFR_MANT (dest);                                              \
    if (MPFR_UNLIKELY (_destprec >= _srcprec))                              \
      {                                                                     \
        _srcs  = MPFR_PREC2LIMBS (_srcprec);                                \
        _dests = MPFR_PREC2LIMBS (_destprec) - _srcs;                       \
        MPN_COPY (_destp + _dests, srcp, _srcs);                            \
        MPN_ZERO (_destp, _dests);                                          \
        inexact = 0;                                                        \
      }                                                                     \
    else                                                                    \
      {                                                                     \
        /* Non trivial case: rounding needed */                             \
        mpfr_prec_t _sh;                                                    \
        mp_limb_t *_sp;                                                     \
        mp_limb_t _rb, _sb, _ulp;                                           \
                                                                            \
        /* Compute Position and shift */                                    \
        _srcs  = MPFR_PREC2LIMBS (_srcprec);                                \
        _dests = MPFR_PREC2LIMBS (_destprec);                               \
        MPFR_UNSIGNED_MINUS_MODULO (_sh, _destprec);                        \
        _sp = (srcp) + _srcs - _dests;                                      \
                                                                            \
        /* General case when prec % GMP_NUMB_BITS != 0 */                   \
        if (MPFR_LIKELY (_sh != 0))                                         \
          {                                                                 \
            mp_limb_t _mask;                                                \
            /* Compute Rounding Bit and Sticky Bit */                       \
            /* Note: in directed rounding modes, if the rounding bit */     \
            /* is 1, the behavior does not depend on the sticky bit; */     \
            /* thus we will not try to compute it in this case (this */     \
            /* can be much faster and avoids to read uninitialized   */     \
            /* data in the current mpfr_mul implementation). We just */     \
            /* make sure that _sb is initialized.                    */     \
            _mask = MPFR_LIMB_ONE << (_sh - 1);                             \
            _rb = _sp[0] & _mask;                                           \
            _sb = _sp[0] & (_mask - 1);                                     \
            if (MPFR_UNLIKELY (_sb == 0) &&                                 \
                ((rnd) == MPFR_RNDN || _rb == 0))                           \
              { /* TODO: Improve it */                                      \
                mp_limb_t *_tmp;                                            \
                mp_size_t _n;                                               \
                for (_tmp = _sp, _n = _srcs - _dests ;                      \
                     _n != 0 && _sb == 0 ; _n--)                            \
                  _sb = *--_tmp;                                            \
              }                                                             \
            _ulp = 2 * _mask;                                               \
          }                                                                 \
        else /* _sh == 0 */                                                 \
          {                                                                 \
            MPFR_ASSERTD (_dests < _srcs);                                  \
            /* Compute Rounding Bit and Sticky Bit - see note above */      \
            _rb = _sp[-1] & MPFR_LIMB_HIGHBIT;                              \
            _sb = _sp[-1] & (MPFR_LIMB_HIGHBIT-1);                          \
            if (MPFR_UNLIKELY (_sb == 0) &&                                 \
                ((rnd) == MPFR_RNDN || _rb == 0))                           \
              {                                                             \
                mp_limb_t *_tmp;                                            \
                mp_size_t _n;                                               \
                for (_tmp = _sp - 1, _n = _srcs - _dests - 1 ;              \
                     _n != 0 && _sb == 0 ; _n--)                            \
                  _sb = *--_tmp;                                            \
              }                                                             \
            _ulp = MPFR_LIMB_ONE;                                           \
          }                                                                 \
        /* Rounding */                                                      \
        if (MPFR_LIKELY (rnd == MPFR_RNDN))                                 \
          {                                                                 \
            if (_rb == 0)                                                   \
              {                                                             \
              trunc:                                                        \
                inexact = MPFR_LIKELY ((_sb | _rb) != 0) ? -sign : 0;       \
              trunc_doit:                                                   \
                MPN_COPY (_destp, _sp, _dests);                             \
                _destp[0] &= ~(_ulp - 1);                                   \
              }                                                             \
            else if (MPFR_UNLIKELY (_sb == 0))                              \
              { /* Middle of two consecutive representable numbers */       \
                MIDDLE_HANDLER;                                             \
              }                                                             \
            else                                                            \
              {                                                             \
                if (0)                                                      \
                  goto addoneulp_doit; /* dummy code to avoid warning */    \
              addoneulp:                                                    \
                inexact = sign;                                             \
              addoneulp_doit:                                               \
                if (MPFR_UNLIKELY (mpn_add_1 (_destp, _sp, _dests, _ulp)))  \
                  {                                                         \
                    _destp[_dests - 1] = MPFR_LIMB_HIGHBIT;                 \
                    OVERFLOW_HANDLER;                                       \
                  }                                                         \
                _destp[0] &= ~(_ulp - 1);                                   \
              }                                                             \
          }                                                                 \
        else                                                                \
          { /* Directed rounding mode */                                    \
            if (MPFR_LIKELY (MPFR_IS_LIKE_RNDZ (rnd,                        \
                                                MPFR_IS_NEG_SIGN (sign))))  \
              goto trunc;                                                   \
             else if (MPFR_UNLIKELY ((_sb | _rb) == 0))                     \
               {                                                            \
                 inexact = 0;                                               \
                 goto trunc_doit;                                           \
               }                                                            \
             else                                                           \
              goto addoneulp;                                               \
          }                                                                 \
      }                                                                     \
  } while (0)

/*
 * Round mantissa (srcp, sprec) to mpfr_t dest using rounding mode rnd
 * assuming dest's sign is sign.
 * Execute OVERFLOW_HANDLER in case of overflow when rounding.
 */
#define MPFR_RNDRAW(inexact, dest, srcp, sprec, rnd, sign, OVERFLOW_HANDLER) \
   MPFR_RNDRAW_GEN (inexact, dest, srcp, sprec, rnd, sign,                   \
     if ((_sp[0] & _ulp) == 0)                                               \
       {                                                                     \
         inexact = -sign;                                                    \
         goto trunc_doit;                                                    \
       }                                                                     \
     else                                                                    \
       goto addoneulp;                                                       \
     , OVERFLOW_HANDLER)

/*
 * Round mantissa (srcp, sprec) to mpfr_t dest using rounding mode rnd
 * assuming dest's sign is sign.
 * Execute OVERFLOW_HANDLER in case of overflow when rounding.
 * Set inexact to +/- MPFR_EVEN_INEX in case of even rounding.
 */
#define MPFR_RNDRAW_EVEN(inexact, dest, srcp, sprec, rnd, sign, \
                         OVERFLOW_HANDLER)                      \
   MPFR_RNDRAW_GEN (inexact, dest, srcp, sprec, rnd, sign,      \
     if ((_sp[0] & _ulp) == 0)                                  \
       {                                                        \
         inexact = -MPFR_EVEN_INEX * sign;                      \
         goto trunc_doit;                                       \
       }                                                        \
     else                                                       \
       {                                                        \
         inexact = MPFR_EVEN_INEX * sign;                       \
         goto addoneulp_doit;                                   \
       }                                                        \
     , OVERFLOW_HANDLER)

/* Return TRUE if b is non singular and we can round it to precision 'prec'
   and determine the ternary value, with rounding mode 'rnd', and with
   error at most 'error' */
#define MPFR_CAN_ROUND(b,err,prec,rnd)                                       \
 (!MPFR_IS_SINGULAR (b) && mpfr_round_p (MPFR_MANT (b), MPFR_LIMB_SIZE (b),  \
                                         (err), (prec) + ((rnd)==MPFR_RNDN)))

/* Copy the sign and the significand, and handle the exponent in exp. */
#define MPFR_SETRAW(inexact,dest,src,exp,rnd)                           \
  if (MPFR_UNLIKELY (dest != src))                                      \
    {                                                                   \
      MPFR_SET_SIGN (dest, MPFR_SIGN (src));                            \
      if (MPFR_LIKELY (MPFR_PREC (dest) == MPFR_PREC (src)))            \
        {                                                               \
          MPN_COPY (MPFR_MANT (dest), MPFR_MANT (src),                  \
                    MPFR_LIMB_SIZE (src));                              \
          inexact = 0;                                                  \
        }                                                               \
      else                                                              \
        {                                                               \
          MPFR_RNDRAW (inexact, dest, MPFR_MANT (src), MPFR_PREC (src), \
                       rnd, MPFR_SIGN (src), exp++);                    \
        }                                                               \
    }                                                                   \
  else                                                                  \
    inexact = 0;

/* TODO: fix this description (see round_near_x.c). */
/* Assuming that the function has a Taylor expansion which looks like:
    y=o(f(x)) = o(v + g(x)) with |g(x)| <= 2^(EXP(v)-err)
   we can quickly set y to v if x is small (ie err > prec(y)+1) in most
   cases. It assumes that f(x) is not representable exactly as a FP number.
   v must not be a singular value (NAN, INF or ZERO); usual values are
   v=1 or v=x.

   y is the destination (a mpfr_t), v the value to set (a mpfr_t),
   err1+err2 with err2 <= 3 the error term (mpfr_exp_t's), dir (an int) is
   the direction of the committed error (if dir = 0, it rounds toward 0,
   if dir=1, it rounds away from 0), rnd the rounding mode.

   It returns from the function a ternary value in case of success.
   If you want to free something, you must fill the "extra" field
   in consequences, otherwise put nothing in it.

   The test is less restrictive than necessary, but the function
   will finish the check itself.

   Note: err1 + err2 is allowed to overflow as mpfr_exp_t, but it must give
   its real value as mpfr_uexp_t.
*/
#define MPFR_FAST_COMPUTE_IF_SMALL_INPUT(y,v,err1,err2,dir,rnd,extra)   \
  do {                                                                  \
    mpfr_ptr _y = (y);                                                  \
    mpfr_exp_t _err1 = (err1);                                          \
    mpfr_exp_t _err2 = (err2);                                          \
    if (_err1 > 0)                                                      \
      {                                                                 \
        mpfr_uexp_t _err = (mpfr_uexp_t) _err1 + _err2;                 \
        if (MPFR_UNLIKELY (_err > MPFR_PREC (_y) + 1))                  \
          {                                                             \
            int _inexact = mpfr_round_near_x (_y,(v),_err,(dir),(rnd)); \
            if (_inexact != 0)                                          \
              {                                                         \
                extra;                                                  \
                return _inexact;                                        \
              }                                                         \
          }                                                             \
      }                                                                 \
  } while (0)

/* Variant, to be called somewhere after MPFR_SAVE_EXPO_MARK. This variant
   is needed when there are some computations before or when some non-zero
   real constant is used, such as __gmpfr_one for mpfr_cos. */
#define MPFR_SMALL_INPUT_AFTER_SAVE_EXPO(y,v,err1,err2,dir,rnd,expo,extra) \
  do {                                                                  \
    mpfr_ptr _y = (y);                                                  \
    mpfr_exp_t _err1 = (err1);                                          \
    mpfr_exp_t _err2 = (err2);                                          \
    if (_err1 > 0)                                                      \
      {                                                                 \
        mpfr_uexp_t _err = (mpfr_uexp_t) _err1 + _err2;                 \
        if (MPFR_UNLIKELY (_err > MPFR_PREC (_y) + 1))                  \
          {                                                             \
            int _inexact;                                               \
            MPFR_CLEAR_FLAGS ();                                        \
            _inexact = mpfr_round_near_x (_y,(v),_err,(dir),(rnd));     \
            if (_inexact != 0)                                          \
              {                                                         \
                extra;                                                  \
                MPFR_SAVE_EXPO_UPDATE_FLAGS (expo, __gmpfr_flags);      \
                MPFR_SAVE_EXPO_FREE (expo);                             \
                return mpfr_check_range (_y, _inexact, (rnd));          \
              }                                                         \
          }                                                             \
      }                                                                 \
  } while (0)

/******************************************************
 ***************  Ziv Loop Macro  *********************
 ******************************************************/

/* To safely increase some precision, detecting integer overflows.
   This macro is particularly useful when determining the initial
   working precision before Ziv's loop. P is a precision, X is an
   arbitrary nonnegative integer.
   Note: On 2012-02-23, the MPFR_PREC_MAX value has been decreased
   by 256 from the maximum value representable in the mpfr_prec_t
   type, in order to avoid some integer overflows when this macro
   is not used (if the result is larger than MPFR_PREC_MAX, this
   should be detected with a later assertion, e.g. in mpfr_init2).
   But this change is mainly for existing code that has not been
   updated yet. So, it is advised to always use MPFR_ADD_PREC if
   the result can be larger than MPFR_PREC_MAX. */
#define MPFR_ADD_PREC(P,X) \
  (MPFR_ASSERTN ((X) <= MPFR_PREC_MAX - (P)), (P) + (X))

#ifndef MPFR_USE_LOGGING

#define MPFR_ZIV_DECL(_x) mpfr_prec_t _x
#define MPFR_ZIV_INIT(_x, _p) (_x) = GMP_NUMB_BITS
#define MPFR_ZIV_NEXT(_x, _p) ((_p) = MPFR_ADD_PREC (_p, _x), (_x) = (_p)/2)
#define MPFR_ZIV_FREE(x)

#else

/* The following test on glibc is there mainly for Darwin (Mac OS X), to
   obtain a better error message. The real test should have been a test
   concerning nested functions in gcc, which are disabled by default on
   Darwin; but it is not possible to do that without a configure test. */
# if defined (__cplusplus) || !(__MPFR_GNUC(3,0) && __MPFR_GLIBC(2,0))
#  error "Logging not supported (needs gcc >= 3.0 and GNU C Library >= 2.0)."
# endif

/* Use LOGGING */

/* Note: the mpfr_log_level >= 0 below avoids to take into account
   Ziv loops used by the MPFR functions called by the mpfr_fprintf
   in LOG_PRINT. */

#define MPFR_ZIV_DECL(_x)                                               \
  mpfr_prec_t _x;                                                       \
  int _x ## _cpt = 1;                                                   \
  static unsigned long  _x ## _loop = 0, _x ## _bad = 0;                \
  static const char *_x ## _fname = __func__;                           \
  auto void __attribute__ ((destructor)) x ## _f  (void);               \
  void __attribute__ ((destructor)) x ## _f  (void) {                   \
    if (_x ## _loop != 0 && (MPFR_LOG_STAT_F & mpfr_log_type)) {        \
      fprintf (mpfr_log_file,                                           \
               "%s: Ziv failed %2.2f%% (%lu bad cases / %lu calls)\n",  \
               _x ## _fname, (double) 100.0 * _x ## _bad / _x ## _loop, \
               _x ## _bad, _x ## _loop );                               \
      if (mpfr_log_flush)                                               \
        fflush (mpfr_log_file);                                         \
    }                                                                   \
  }

#define MPFR_ZIV_INIT(_x, _p)                                           \
  do                                                                    \
    {                                                                   \
      (_x) = GMP_NUMB_BITS;                                             \
      if (mpfr_log_level >= 0)                                          \
        _x ## _loop ++;                                                 \
      LOG_PRINT (MPFR_LOG_BADCASE_F, "%s:ZIV 1st prec=%Pd\n",           \
                 __func__, (mpfr_prec_t) (_p));                         \
    }                                                                   \
  while (0)

#define MPFR_ZIV_NEXT(_x, _p)                                           \
  do                                                                    \
    {                                                                   \
      (_p) = MPFR_ADD_PREC (_p, _x);                                    \
      (_x) = (_p) / 2;                                                  \
      if (mpfr_log_level >= 0)                                          \
        _x ## _bad += (_x ## _cpt == 1);                                \
      _x ## _cpt ++;                                                    \
      LOG_PRINT (MPFR_LOG_BADCASE_F, "%s:ZIV new prec=%Pd\n",           \
                 __func__, (mpfr_prec_t) (_p));                         \
    }                                                                   \
  while (0)

#define MPFR_ZIV_FREE(_x)                                               \
  do                                                                    \
    if (_x ## _cpt > 1)                                                 \
      LOG_PRINT (MPFR_LOG_BADCASE_F, "%s:ZIV %d loops\n",               \
                 __func__, _x ## _cpt);                                 \
  while (0)

#endif


/******************************************************
 ***************  Logging Macros  *********************
 ******************************************************/

/* The different kind of LOG */
#define MPFR_LOG_INPUT_F    1
#define MPFR_LOG_OUTPUT_F   2
#define MPFR_LOG_INTERNAL_F 4
#define MPFR_LOG_TIME_F     8
#define MPFR_LOG_BADCASE_F  16
#define MPFR_LOG_MSG_F      32
#define MPFR_LOG_STAT_F     64

#ifdef MPFR_USE_LOGGING

/* Check if we can support this feature */
# ifdef MPFR_USE_THREAD_SAFE
#  error "Enable either `Logging' or `thread-safe', not both"
# endif
# if !__MPFR_GNUC(3,0)
#  error "Logging not supported (GCC >= 3.0)"
# endif

#if defined (__cplusplus)
extern "C" {
#endif

__MPFR_DECLSPEC extern FILE *mpfr_log_file;
__MPFR_DECLSPEC extern int   mpfr_log_flush;
__MPFR_DECLSPEC extern int   mpfr_log_type;
__MPFR_DECLSPEC extern int   mpfr_log_level;
__MPFR_DECLSPEC extern int   mpfr_log_current;
__MPFR_DECLSPEC extern mpfr_prec_t mpfr_log_prec;

#if defined (__cplusplus)
 }
#endif

/* LOG_PRINT calls mpfr_fprintf on mpfr_log_file with logging disabled
   (recursive logging is not wanted and freezes MPFR). */
#define LOG_PRINT(type, format, ...)                                    \
  do                                                                    \
    if ((mpfr_log_type & (type)) && mpfr_log_current <= mpfr_log_level) \
      {                                                                 \
        int old_level = mpfr_log_level;                                 \
        mpfr_log_level = -1;  /* disable logging in mpfr_fprintf */     \
        __gmpfr_cache_const_pi = __gmpfr_logging_pi;                    \
        __gmpfr_cache_const_log2 = __gmpfr_logging_log2;                \
        mpfr_fprintf (mpfr_log_file, format, __VA_ARGS__);              \
        if (mpfr_log_flush)                                             \
          fflush (mpfr_log_file);                                       \
        mpfr_log_level = old_level;                                     \
        __gmpfr_cache_const_pi = __gmpfr_normal_pi;                     \
        __gmpfr_cache_const_log2 = __gmpfr_normal_log2;                 \
      }                                                                 \
  while (0)

#define MPFR_LOG_VAR(x)                                                 \
  LOG_PRINT (MPFR_LOG_INTERNAL_F, "%s.%d:%s[%#Pu]=%.*Rg\n", __func__,   \
             __LINE__, #x, mpfr_get_prec (x), mpfr_log_prec, x)

#define MPFR_LOG_MSG2(format, ...)                                      \
  LOG_PRINT (MPFR_LOG_MSG_F, "%s.%d: "format, __func__, __LINE__, __VA_ARGS__)
#define MPFR_LOG_MSG(x) MPFR_LOG_MSG2 x

#define MPFR_LOG_BEGIN2(format, ...)                                    \
  mpfr_log_current ++;                                                  \
  LOG_PRINT (MPFR_LOG_INPUT_F, "%s:IN  "format"\n", __func__, __VA_ARGS__); \
  if ((MPFR_LOG_TIME_F & mpfr_log_type) &&                              \
      (mpfr_log_current <= mpfr_log_level))                             \
    __gmpfr_log_time = mpfr_get_cputime ();
#define MPFR_LOG_BEGIN(x)                                               \
  int __gmpfr_log_time = 0;                                             \
  MPFR_LOG_BEGIN2 x

#define MPFR_LOG_END2(format, ...)                                      \
  LOG_PRINT (MPFR_LOG_TIME_F, "%s:TIM %dms\n", __mpfr_log_fname,        \
             mpfr_get_cputime () - __gmpfr_log_time);                   \
  LOG_PRINT (MPFR_LOG_OUTPUT_F, "%s:OUT "format"\n", __mpfr_log_fname,  \
             __VA_ARGS__);                                              \
  mpfr_log_current --;
#define MPFR_LOG_END(x)                                                 \
  static const char *__mpfr_log_fname = __func__;                       \
  MPFR_LOG_END2 x

#define MPFR_LOG_FUNC(begin,end)                                        \
  static const char *__mpfr_log_fname = __func__;                       \
  auto void __mpfr_log_cleanup (int *time);                             \
  void __mpfr_log_cleanup (int *time) {                                 \
    int __gmpfr_log_time = *time;                                       \
    MPFR_LOG_END2 end; }                                                \
  int __gmpfr_log_time __attribute__ ((cleanup (__mpfr_log_cleanup)));  \
  __gmpfr_log_time = 0;                                                 \
  MPFR_LOG_BEGIN2 begin

#else /* MPFR_USE_LOGGING */

/* Define void macro for logging */

#define MPFR_LOG_VAR(x)
#define MPFR_LOG_BEGIN(x)
#define MPFR_LOG_END(x)
#define MPFR_LOG_MSG(x)
#define MPFR_LOG_FUNC(x,y)

#endif /* MPFR_USE_LOGGING */


/**************************************************************
 ************  Group Initialize Functions Macros  *************
 **************************************************************/

#ifndef MPFR_GROUP_STATIC_SIZE
# define MPFR_GROUP_STATIC_SIZE 16
#endif

struct mpfr_group_t {
  size_t     alloc;
  mp_limb_t *mant;
  mp_limb_t  tab[MPFR_GROUP_STATIC_SIZE];
};

#define MPFR_GROUP_DECL(g) struct mpfr_group_t g
#define MPFR_GROUP_CLEAR(g) do {                                 \
 MPFR_LOG_MSG (("GROUP_CLEAR: ptr = 0x%lX, size = %lu\n",        \
                (unsigned long) (g).mant,                        \
                (unsigned long) (g).alloc));                     \
 if (MPFR_UNLIKELY ((g).alloc != 0)) {                           \
   MPFR_ASSERTD ((g).mant != (g).tab);                           \
   (*__gmp_free_func) ((g).mant, (g).alloc);                     \
 }} while (0)

#define MPFR_GROUP_INIT_TEMPLATE(g, prec, num, handler) do {            \
 mpfr_prec_t _prec = (prec);                                            \
 mp_size_t _size;                                                       \
 MPFR_ASSERTD (_prec >= MPFR_PREC_MIN);                                 \
 if (MPFR_UNLIKELY (_prec > MPFR_PREC_MAX))                             \
   mpfr_abort_prec_max ();                                              \
 _size = MPFR_PREC2LIMBS (_prec);                                       \
 if (MPFR_UNLIKELY (_size * (num) > MPFR_GROUP_STATIC_SIZE))            \
   {                                                                    \
     (g).alloc = (num) * _size * sizeof (mp_limb_t);                    \
     (g).mant = (mp_limb_t *) (*__gmp_allocate_func) ((g).alloc);       \
   }                                                                    \
 else                                                                   \
   {                                                                    \
     (g).alloc = 0;                                                     \
     (g).mant = (g).tab;                                                \
   }                                                                    \
 MPFR_LOG_MSG (("GROUP_INIT: ptr = 0x%lX, size = %lu\n",                \
                (unsigned long) (g).mant, (unsigned long) (g).alloc));  \
 handler;                                                               \
 } while (0)
#define MPFR_GROUP_TINIT(g, n, x)                       \
  MPFR_TMP_INIT1 ((g).mant + _size * (n), x, _prec)

#define MPFR_GROUP_INIT_1(g, prec, x)                            \
 MPFR_GROUP_INIT_TEMPLATE(g, prec, 1, MPFR_GROUP_TINIT(g, 0, x))
#define MPFR_GROUP_INIT_2(g, prec, x, y)                         \
 MPFR_GROUP_INIT_TEMPLATE(g, prec, 2,                            \
   MPFR_GROUP_TINIT(g, 0, x);MPFR_GROUP_TINIT(g, 1, y))
#define MPFR_GROUP_INIT_3(g, prec, x, y, z)                      \
 MPFR_GROUP_INIT_TEMPLATE(g, prec, 3,                            \
   MPFR_GROUP_TINIT(g, 0, x);MPFR_GROUP_TINIT(g, 1, y);          \
   MPFR_GROUP_TINIT(g, 2, z))
#define MPFR_GROUP_INIT_4(g, prec, x, y, z, t)                   \
 MPFR_GROUP_INIT_TEMPLATE(g, prec, 4,                            \
   MPFR_GROUP_TINIT(g, 0, x);MPFR_GROUP_TINIT(g, 1, y);          \
   MPFR_GROUP_TINIT(g, 2, z);MPFR_GROUP_TINIT(g, 3, t))
#define MPFR_GROUP_INIT_5(g, prec, x, y, z, t, a)                \
 MPFR_GROUP_INIT_TEMPLATE(g, prec, 5,                            \
   MPFR_GROUP_TINIT(g, 0, x);MPFR_GROUP_TINIT(g, 1, y);          \
   MPFR_GROUP_TINIT(g, 2, z);MPFR_GROUP_TINIT(g, 3, t);          \
   MPFR_GROUP_TINIT(g, 4, a))
#define MPFR_GROUP_INIT_6(g, prec, x, y, z, t, a, b)             \
 MPFR_GROUP_INIT_TEMPLATE(g, prec, 6,                            \
   MPFR_GROUP_TINIT(g, 0, x);MPFR_GROUP_TINIT(g, 1, y);          \
   MPFR_GROUP_TINIT(g, 2, z);MPFR_GROUP_TINIT(g, 3, t);          \
   MPFR_GROUP_TINIT(g, 4, a);MPFR_GROUP_TINIT(g, 5, b))

#define MPFR_GROUP_REPREC_TEMPLATE(g, prec, num, handler) do {          \
 mpfr_prec_t _prec = (prec);                                            \
 size_t    _oalloc = (g).alloc;                                         \
 mp_size_t _size;                                                       \
 MPFR_LOG_MSG (("GROUP_REPREC: oldptr = 0x%lX, oldsize = %lu\n",        \
                (unsigned long) (g).mant, (unsigned long) _oalloc));    \
 MPFR_ASSERTD (_prec >= MPFR_PREC_MIN);                                 \
 if (MPFR_UNLIKELY (_prec > MPFR_PREC_MAX))                             \
   mpfr_abort_prec_max ();                                              \
 _size = MPFR_PREC2LIMBS (_prec);                                       \
 (g).alloc = (num) * _size * sizeof (mp_limb_t);                        \
 if (MPFR_LIKELY (_oalloc == 0))                                        \
   (g).mant = (mp_limb_t *) (*__gmp_allocate_func) ((g).alloc);         \
 else                                                                   \
   (g).mant = (mp_limb_t *)                                             \
     (*__gmp_reallocate_func) ((g).mant, _oalloc, (g).alloc);           \
 MPFR_LOG_MSG (("GROUP_REPREC: newptr = 0x%lX, newsize = %lu\n",        \
                (unsigned long) (g).mant, (unsigned long) (g).alloc));  \
 handler;                                                               \
 } while (0)

#define MPFR_GROUP_REPREC_1(g, prec, x)                          \
 MPFR_GROUP_REPREC_TEMPLATE(g, prec, 1, MPFR_GROUP_TINIT(g, 0, x))
#define MPFR_GROUP_REPREC_2(g, prec, x, y)                       \
 MPFR_GROUP_REPREC_TEMPLATE(g, prec, 2,                          \
   MPFR_GROUP_TINIT(g, 0, x);MPFR_GROUP_TINIT(g, 1, y))
#define MPFR_GROUP_REPREC_3(g, prec, x, y, z)                    \
 MPFR_GROUP_REPREC_TEMPLATE(g, prec, 3,                          \
   MPFR_GROUP_TINIT(g, 0, x);MPFR_GROUP_TINIT(g, 1, y);          \
   MPFR_GROUP_TINIT(g, 2, z))
#define MPFR_GROUP_REPREC_4(g, prec, x, y, z, t)                 \
 MPFR_GROUP_REPREC_TEMPLATE(g, prec, 4,                          \
   MPFR_GROUP_TINIT(g, 0, x);MPFR_GROUP_TINIT(g, 1, y);          \
   MPFR_GROUP_TINIT(g, 2, z);MPFR_GROUP_TINIT(g, 3, t))
#define MPFR_GROUP_REPREC_5(g, prec, x, y, z, t, a)              \
 MPFR_GROUP_REPREC_TEMPLATE(g, prec, 5,                          \
   MPFR_GROUP_TINIT(g, 0, x);MPFR_GROUP_TINIT(g, 1, y);          \
   MPFR_GROUP_TINIT(g, 2, z);MPFR_GROUP_TINIT(g, 3, t);          \
   MPFR_GROUP_TINIT(g, 4, a))
#define MPFR_GROUP_REPREC_6(g, prec, x, y, z, t, a, b)           \
 MPFR_GROUP_REPREC_TEMPLATE(g, prec, 6,                          \
   MPFR_GROUP_TINIT(g, 0, x);MPFR_GROUP_TINIT(g, 1, y);          \
   MPFR_GROUP_TINIT(g, 2, z);MPFR_GROUP_TINIT(g, 3, t);          \
   MPFR_GROUP_TINIT(g, 4, a);MPFR_GROUP_TINIT(g, 5, b))


/******************************************************
 ***************  Internal Functions  *****************
 ******************************************************/

#if defined (__cplusplus)
extern "C" {
#endif

MPFR_COLD_FUNCTION_ATTR __MPFR_DECLSPEC int
  mpfr_underflow _MPFR_PROTO ((mpfr_ptr, mpfr_rnd_t, int));
MPFR_COLD_FUNCTION_ATTR __MPFR_DECLSPEC int
  mpfr_overflow _MPFR_PROTO ((mpfr_ptr, mpfr_rnd_t, int));

__MPFR_DECLSPEC int mpfr_add1 _MPFR_PROTO ((mpfr_ptr, mpfr_srcptr,
                                            mpfr_srcptr, mpfr_rnd_t));
__MPFR_DECLSPEC int mpfr_sub1 _MPFR_PROTO ((mpfr_ptr, mpfr_srcptr,
                                            mpfr_srcptr, mpfr_rnd_t));
__MPFR_DECLSPEC int mpfr_add1sp _MPFR_PROTO ((mpfr_ptr, mpfr_srcptr,
                                              mpfr_srcptr, mpfr_rnd_t));
__MPFR_DECLSPEC int mpfr_sub1sp _MPFR_PROTO ((mpfr_ptr, mpfr_srcptr,
                                              mpfr_srcptr, mpfr_rnd_t));
__MPFR_DECLSPEC int mpfr_can_round_raw _MPFR_PROTO ((const mp_limb_t *,
             mp_size_t, int, mpfr_exp_t, mpfr_rnd_t, mpfr_rnd_t, mpfr_prec_t));

__MPFR_DECLSPEC int mpfr_cmp2 _MPFR_PROTO ((mpfr_srcptr, mpfr_srcptr,
                                            mpfr_prec_t *));

__MPFR_DECLSPEC long          __gmpfr_ceil_log2     _MPFR_PROTO ((double));
__MPFR_DECLSPEC long          __gmpfr_floor_log2    _MPFR_PROTO ((double));
__MPFR_DECLSPEC double        __gmpfr_ceil_exp2     _MPFR_PROTO ((double));
__MPFR_DECLSPEC unsigned long __gmpfr_isqrt     _MPFR_PROTO ((unsigned long));
__MPFR_DECLSPEC unsigned long __gmpfr_cuberoot  _MPFR_PROTO ((unsigned long));
__MPFR_DECLSPEC int       __gmpfr_int_ceil_log2 _MPFR_PROTO ((unsigned long));

__MPFR_DECLSPEC mpfr_exp_t mpfr_ceil_mul _MPFR_PROTO ((mpfr_exp_t, int, int));

__MPFR_DECLSPEC int mpfr_exp_2 _MPFR_PROTO ((mpfr_ptr, mpfr_srcptr,mpfr_rnd_t));
__MPFR_DECLSPEC int mpfr_exp_3 _MPFR_PROTO ((mpfr_ptr, mpfr_srcptr,mpfr_rnd_t));
__MPFR_DECLSPEC int mpfr_powerof2_raw _MPFR_PROTO ((mpfr_srcptr));

__MPFR_DECLSPEC int mpfr_pow_general _MPFR_PROTO ((mpfr_ptr, mpfr_srcptr,
                           mpfr_srcptr, mpfr_rnd_t, int, mpfr_save_expo_t *));

__MPFR_DECLSPEC void mpfr_setmax _MPFR_PROTO ((mpfr_ptr, mpfr_exp_t));
__MPFR_DECLSPEC void mpfr_setmin _MPFR_PROTO ((mpfr_ptr, mpfr_exp_t));

__MPFR_DECLSPEC long mpfr_mpn_exp _MPFR_PROTO ((mp_limb_t *, mpfr_exp_t *, int,
                                                mpfr_exp_t, size_t));

#ifdef _MPFR_H_HAVE_FILE
__MPFR_DECLSPEC void mpfr_fprint_binary _MPFR_PROTO ((FILE *, mpfr_srcptr));
#endif
__MPFR_DECLSPEC void mpfr_print_binary _MPFR_PROTO ((mpfr_srcptr));
__MPFR_DECLSPEC void mpfr_print_mant_binary _MPFR_PROTO ((const char*,
                                          const mp_limb_t*, mpfr_prec_t));
__MPFR_DECLSPEC void mpfr_set_str_binary _MPFR_PROTO((mpfr_ptr, const char*));

__MPFR_DECLSPEC int mpfr_round_raw _MPFR_PROTO ((mp_limb_t *,
       const mp_limb_t *, mpfr_prec_t, int, mpfr_prec_t, mpfr_rnd_t, int *));
__MPFR_DECLSPEC int mpfr_round_raw_2 _MPFR_PROTO ((const mp_limb_t *,
             mpfr_prec_t, int, mpfr_prec_t, mpfr_rnd_t));
/* No longer defined (see round_prec.c).
   Uncomment if it needs to be defined again.
__MPFR_DECLSPEC int mpfr_round_raw_3 _MPFR_PROTO ((const mp_limb_t *,
             mpfr_prec_t, int, mpfr_prec_t, mpfr_rnd_t, int *));
*/
__MPFR_DECLSPEC int mpfr_round_raw_4 _MPFR_PROTO ((mp_limb_t *,
       const mp_limb_t *, mpfr_prec_t, int, mpfr_prec_t, mpfr_rnd_t));

#define mpfr_round_raw2(xp, xn, neg, r, prec) \
  mpfr_round_raw_2((xp),(xn)*GMP_NUMB_BITS,(neg),(prec),(r))

__MPFR_DECLSPEC int mpfr_check _MPFR_PROTO ((mpfr_srcptr));

__MPFR_DECLSPEC int mpfr_sum_sort _MPFR_PROTO ((mpfr_srcptr *const,
                                                unsigned long, mpfr_srcptr *,
                                                mpfr_prec_t *));

__MPFR_DECLSPEC int mpfr_get_cputime _MPFR_PROTO ((void));

__MPFR_DECLSPEC void mpfr_nexttozero _MPFR_PROTO ((mpfr_ptr));
__MPFR_DECLSPEC void mpfr_nexttoinf _MPFR_PROTO ((mpfr_ptr));

__MPFR_DECLSPEC int mpfr_const_pi_internal _MPFR_PROTO ((mpfr_ptr,mpfr_rnd_t));
__MPFR_DECLSPEC int mpfr_const_log2_internal _MPFR_PROTO((mpfr_ptr,mpfr_rnd_t));
__MPFR_DECLSPEC int mpfr_const_euler_internal _MPFR_PROTO((mpfr_ptr, mpfr_rnd_t));
__MPFR_DECLSPEC int mpfr_const_catalan_internal _MPFR_PROTO((mpfr_ptr, mpfr_rnd_t));

#if 0
__MPFR_DECLSPEC void mpfr_init_cache _MPFR_PROTO ((mpfr_cache_t,
                                           int(*)(mpfr_ptr,mpfr_rnd_t)));
#endif
__MPFR_DECLSPEC void mpfr_clear_cache _MPFR_PROTO ((mpfr_cache_t));
__MPFR_DECLSPEC int  mpfr_cache _MPFR_PROTO ((mpfr_ptr, mpfr_cache_t,
                                              mpfr_rnd_t));

__MPFR_DECLSPEC void mpfr_mulhigh_n _MPFR_PROTO ((mpfr_limb_ptr,
                        mpfr_limb_srcptr, mpfr_limb_srcptr, mp_size_t));
__MPFR_DECLSPEC void mpfr_mullow_n  _MPFR_PROTO ((mpfr_limb_ptr,
                        mpfr_limb_srcptr, mpfr_limb_srcptr, mp_size_t));
__MPFR_DECLSPEC void mpfr_sqrhigh_n _MPFR_PROTO ((mpfr_limb_ptr,
                        mpfr_limb_srcptr, mp_size_t));
__MPFR_DECLSPEC mp_limb_t mpfr_divhigh_n _MPFR_PROTO ((mpfr_limb_ptr,
                        mpfr_limb_ptr, mpfr_limb_ptr, mp_size_t));

__MPFR_DECLSPEC int mpfr_round_p _MPFR_PROTO ((mp_limb_t *, mp_size_t,
                                               mpfr_exp_t, mpfr_prec_t));

__MPFR_DECLSPEC int mpfr_round_near_x _MPFR_PROTO ((mpfr_ptr, mpfr_srcptr,
                                                    mpfr_uexp_t, int,
                                                    mpfr_rnd_t));
__MPFR_DECLSPEC MPFR_COLD_FUNCTION_ATTR MPFR_NORETURN void
  mpfr_abort_prec_max _MPFR_PROTO ((void));

__MPFR_DECLSPEC void mpfr_rand_raw _MPFR_PROTO((mpfr_limb_ptr, gmp_randstate_t,
                                                mpfr_prec_t));

__MPFR_DECLSPEC mpz_srcptr mpfr_bernoulli_cache _MPFR_PROTO ((unsigned long));
__MPFR_DECLSPEC void mpfr_bernoulli_freecache _MPFR_PROTO ((void));

__MPFR_DECLSPEC int mpfr_sincos_fast _MPFR_PROTO((mpfr_t, mpfr_t,
                                                  mpfr_srcptr, mpfr_rnd_t));

__MPFR_DECLSPEC double mpfr_scale2 _MPFR_PROTO((double, int));

__MPFR_DECLSPEC void mpfr_div_ui2 _MPFR_PROTO((mpfr_ptr, mpfr_srcptr,
                                               unsigned long int, unsigned long int,
                                               mpfr_rnd_t));

__MPFR_DECLSPEC void mpfr_gamma_one_and_two_third _MPFR_PROTO((mpfr_ptr, mpfr_ptr, mpfr_prec_t));

__MPFR_DECLSPEC void mpfr_mpz_init _MPFR_PROTO((mpz_ptr));
__MPFR_DECLSPEC void mpfr_mpz_clear _MPFR_PROTO((mpz_ptr));

#if defined (__cplusplus)
}
#endif

/******************************************************
 **************  Internal mpz caching *****************
 ******************************************************/

/* Cache for mpz_t */
#if !defined(MPFR_MY_MPZ_INIT) || MPFR_MY_MPZ_INIT != 0
# undef mpz_init
# undef mpz_clear
# define mpz_init mpfr_mpz_init
# define mpz_clear mpfr_mpz_clear
#endif

/******************************************************
 ********** Compute LOG2(LOG2(MPFR_PREC_MAX))**********
 ******************************************************/
#if   _MPFR_PREC_FORMAT == 1
# define MPFR_PREC_MAX_TEMP USHRT_MAX
#elif _MPFR_PREC_FORMAT == 2
# define MPFR_PREC_MAX_TEMP UINT_MAX
#elif _MPFR_PREC_FORMAT == 3
# define MPFR_PREC_MAX_TEMP ULONG_MAX
#endif

#if MPFR_PREC_MAX_TEMP == 255UL
# define MPFR_PREC_BITS 8
# define MPFR_LOG2_PREC_BITS 3
#elif MPFR_PREC_MAX_TEMP == 65535UL
# define MPFR_PREC_BITS 16
# define MPFR_LOG2_PREC_BITS 4
#elif MPFR_PREC_MAX_TEMP == 4294967295UL
# define MPFR_PREC_BITS 32
# define MPFR_LOG2_PREC_BITS 5
#else
/* Assume MPFR_PREC_MAX_TEMP == 18446744073709551615ULL */
# define MPFR_PREC_BITS 64
# define MPFR_LOG2_PREC_BITS 6
#endif


#endif
