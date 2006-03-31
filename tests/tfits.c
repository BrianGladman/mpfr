/* Test file for:
 mpfr_fits_sint_p, mpfr_fits_slong_p, mpfr_fits_sshort_p,
 mpfr_fits_uint_p, mpfr_fits_ulong_p, mpfr_fits_ushort_p

Copyright 2004, 2005, 2006 Free Software Foundation, Inc.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#if HAVE_CONFIG_H
# include "config.h"       /* for a build within gmp */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

/* The ISO C99 standard specifies that in C++ implementations the
   INTMAX_MAX, ... macros should only be defined if explicitly requested.  */
#if defined __cplusplus
# define __STDC_LIMIT_MACROS
# define __STDC_CONSTANT_MACROS
#endif

#ifdef HAVE_STDINT_H
# include <stdint.h>
#endif
#ifdef HAVE_INTTYPES_H
# include <inttypes.h>
#endif

#include "mpfr-test.h"

#define ERROR1 { printf("Initial error for x="); mpfr_dump(x); exit(1); }
#define ERROR2 { printf("Error for x="); mpfr_dump(x); exit(1); }

static void check_intmax (void);

int
main (void)
{
  mpfr_t x;

  tests_start_mpfr ();

  mpfr_init2 (x, 256);

  /* Check NAN */
  mpfr_set_nan(x);
  if (mpfr_fits_ulong_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_slong_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_uint_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_sint_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_ushort_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_sshort_p(x, GMP_RNDN))
    ERROR1;

  /* Check INF */
  mpfr_set_inf(x, 1);
  if (mpfr_fits_ulong_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_slong_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_uint_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_sint_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_ushort_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_sshort_p(x, GMP_RNDN))
    ERROR1;

  /* Check Zero */
  MPFR_SET_ZERO(x);
  if (!mpfr_fits_ulong_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_slong_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_uint_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_sint_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_ushort_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_sshort_p(x, GMP_RNDN))
    ERROR2;

  /* Check small op */
  mpfr_set_str1 (x, "1@-1");
  if (!mpfr_fits_ulong_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_slong_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_uint_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_sint_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_ushort_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_sshort_p(x, GMP_RNDN))
    ERROR2;

  /* Check 17 */
  mpfr_set_ui (x, 17, GMP_RNDN);
  if (!mpfr_fits_ulong_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_slong_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_uint_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_sint_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_ushort_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_sshort_p(x, GMP_RNDN))
    ERROR2;

  /* Check all other values */
  mpfr_set_ui(x, ULONG_MAX, GMP_RNDN);
  mpfr_mul_2exp(x, x, 1, GMP_RNDN);
  if (mpfr_fits_ulong_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_slong_p(x, GMP_RNDN))
    ERROR1;
  mpfr_mul_2exp(x, x, 40, GMP_RNDN);
  if (mpfr_fits_ulong_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_uint_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_sint_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_ushort_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_sshort_p(x, GMP_RNDN))
    ERROR1;

  mpfr_set_ui(x, ULONG_MAX, GMP_RNDN);
  if (!mpfr_fits_ulong_p(x, GMP_RNDN))
    ERROR2;
  mpfr_set_ui(x, LONG_MAX, GMP_RNDN);
  if (!mpfr_fits_slong_p(x, GMP_RNDN))
    ERROR2;
  mpfr_set_ui(x, UINT_MAX, GMP_RNDN);
  if (!mpfr_fits_uint_p(x, GMP_RNDN))
    ERROR2;
  mpfr_set_ui(x, INT_MAX, GMP_RNDN);
  if (!mpfr_fits_sint_p(x, GMP_RNDN))
    ERROR2;
  mpfr_set_ui(x, USHRT_MAX, GMP_RNDN);
  if (!mpfr_fits_ushort_p(x, GMP_RNDN))
    ERROR2;
  mpfr_set_ui(x, SHRT_MAX, GMP_RNDN);
  if (!mpfr_fits_sshort_p(x, GMP_RNDN))
    ERROR2;

  mpfr_set_si(x, 1, GMP_RNDN);
  if (!mpfr_fits_sint_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_sshort_p(x, GMP_RNDN))
    ERROR2;

  /* Check negative value */
  mpfr_set_si (x, -1, GMP_RNDN);
  if (!mpfr_fits_sint_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_sshort_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_slong_p(x, GMP_RNDN))
    ERROR2;
  if (mpfr_fits_uint_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_ushort_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_ulong_p(x, GMP_RNDN))
    ERROR1;

  mpfr_clear (x);

  check_intmax ();

  tests_end_mpfr ();
  return 0;
}

static void check_intmax (void)
{
#ifdef _MPFR_H_HAVE_INTMAX_T
  mpfr_t x;

  mpfr_init2 (x, sizeof(uintmax_t)*CHAR_BIT);

  /* Check NAN */
  mpfr_set_nan(x);
  if (mpfr_fits_uintmax_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_intmax_p(x, GMP_RNDN))
    ERROR1;

  /* Check INF */
  mpfr_set_inf(x, 1);
  if (mpfr_fits_uintmax_p(x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_intmax_p(x, GMP_RNDN))
    ERROR1;

  /* Check Zero */
  MPFR_SET_ZERO(x);
  if (!mpfr_fits_uintmax_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_intmax_p(x, GMP_RNDN))
    ERROR2;

  /* Check small op */
  mpfr_set_str1 (x, "1@-1");
  if (!mpfr_fits_uintmax_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_intmax_p(x, GMP_RNDN))
    ERROR2;

  /* Check 17 */
  mpfr_set_ui (x, 17, GMP_RNDN);
  if (!mpfr_fits_uintmax_p(x, GMP_RNDN))
    ERROR2;
  if (!mpfr_fits_intmax_p(x, GMP_RNDN))
    ERROR2;

  /* Check hugest */
  mpfr_set_ui_2exp (x, 42, sizeof (uintmax_t) * 32, GMP_RNDN);
  if (mpfr_fits_uintmax_p (x, GMP_RNDN))
    ERROR1;
  if (mpfr_fits_intmax_p (x, GMP_RNDN))
    ERROR1;

  /* Check all other values */
  mpfr_set_uj (x, UINTMAX_MAX, GMP_RNDN);
  mpfr_add_ui (x, x, 1, GMP_RNDN);
  if (mpfr_fits_uintmax_p (x, GMP_RNDN))
    ERROR1;
  mpfr_set_uj (x, UINTMAX_MAX, GMP_RNDN);
  if (!mpfr_fits_uintmax_p (x, GMP_RNDN))
    ERROR2;
  mpfr_set_sj (x, INTMAX_MAX, GMP_RNDN);
  mpfr_add_ui (x, x, 1, GMP_RNDN);
  if (mpfr_fits_intmax_p (x, GMP_RNDN))
    ERROR1;
  mpfr_set_sj (x, INTMAX_MAX, GMP_RNDN);
  if (!mpfr_fits_intmax_p (x, GMP_RNDN))
    ERROR2;
  mpfr_set_sj (x, INTMAX_MIN, GMP_RNDN);
  if (!mpfr_fits_intmax_p (x, GMP_RNDN))
    ERROR2;
  mpfr_sub_ui (x, x, 1, GMP_RNDN);
  if (mpfr_fits_intmax_p (x, GMP_RNDN))
    ERROR1;

  /* Check negative value */
  mpfr_set_si (x, -1, GMP_RNDN);
  if (!mpfr_fits_intmax_p (x, GMP_RNDN))
    ERROR2;
  if (mpfr_fits_uintmax_p (x, GMP_RNDN))
    ERROR1;

  mpfr_clear (x);
#endif
}

