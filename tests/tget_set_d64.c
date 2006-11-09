/* Test file for mpfr_get_decimal64 and mpfr_set_decimal64.

Copyright 2006 Free Software Foundation, Inc.

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

#include <stdlib.h> /* for exit */
#include "mpfr-test.h"

/* #define DEBUG */

#if MPFR_WANT_DECIMAL_FLOATS
static void
print_decimal64 (_Decimal64 d)
{
  union ieee_double_extract x;
  union ieee_double_decimal64 y;
  unsigned int Gh, i;

  y.d64 = d;
  x.d = y.d;
  Gh = x.s.exp >> 6;
  printf ("|%d%d%d%d%d%d", x.s.sig, Gh >> 4, (Gh >> 3) & 1,
	  (Gh >> 2) & 1, (Gh >> 1) & 1, Gh & 1);
  printf ("%d%d%d%d%d%d", (x.s.exp >> 5) & 1, (x.s.exp >> 4) & 1,
	  (x.s.exp >> 3) & 1, (x.s.exp >> 2) & 1, (x.s.exp >> 1) & 1,
	  x.s.exp & 1);
  for (i = 20; i > 0; i--)
    printf ("%d", (x.s.manh >> (i - 1)) & 1);
  for (i = 32; i > 0; i--)
    printf ("%d", (x.s.manl >> (i - 1)) & 1);
  printf ("|\n");
}

static void
check_inf_nan ()
{
  mpfr_t  x, y;
  _Decimal64 d;

  mpfr_init2 (x, 123);
  mpfr_init2 (y, 123);

  mpfr_set_nan (x);
  d = mpfr_get_decimal64 (x, GMP_RNDZ);
  mpfr_set_ui (x, 1, GMP_RNDZ);
  mpfr_set_decimal64 (x, d, GMP_RNDZ);
  ASSERT_ALWAYS (mpfr_nan_p (x));

  mpfr_set_inf (x, 1);
  d = mpfr_get_decimal64 (x, GMP_RNDZ);
  mpfr_set_ui (x, 1, GMP_RNDZ);
  mpfr_set_decimal64 (x, d, GMP_RNDZ);
  ASSERT_ALWAYS (mpfr_inf_p (x) && mpfr_sgn (x) > 0);

  mpfr_set_inf (x, -1);
  d = mpfr_get_decimal64 (x, GMP_RNDZ);
  mpfr_set_ui (x, 1, GMP_RNDZ);
  mpfr_set_decimal64 (x, d, GMP_RNDZ);
  ASSERT_ALWAYS (mpfr_inf_p (x) && mpfr_sgn (x) < 0);

  mpfr_set_ui (x, 0, GMP_RNDZ);
  d = mpfr_get_decimal64 (x, GMP_RNDZ);
  mpfr_set_ui (x, 1, GMP_RNDZ);
  mpfr_set_decimal64 (x, d, GMP_RNDZ);
  ASSERT_ALWAYS (mpfr_cmp_ui (x, 0) == 0 && MPFR_SIGN (x) > 0);

  mpfr_set_ui (x, 0, GMP_RNDZ);
  mpfr_neg (x, x, GMP_RNDZ);
  d = mpfr_get_decimal64 (x, GMP_RNDZ);
  mpfr_set_ui (x, 1, GMP_RNDZ);
  mpfr_set_decimal64 (x, d, GMP_RNDZ);
  ASSERT_ALWAYS (mpfr_cmp_ui (x, 0) == 0 && MPFR_SIGN (x) < 0);

  mpfr_set_ui (x, 1, GMP_RNDZ);
  d = mpfr_get_decimal64 (x, GMP_RNDZ);
  mpfr_set_ui (x, 0, GMP_RNDZ);
  mpfr_set_decimal64 (x, d, GMP_RNDZ);
  ASSERT_ALWAYS (mpfr_cmp_ui (x, 1) == 0);

  mpfr_set_si (x, -1, GMP_RNDZ);
  d = mpfr_get_decimal64 (x, GMP_RNDZ);
  mpfr_set_ui (x, 0, GMP_RNDZ);
  mpfr_set_decimal64 (x, d, GMP_RNDZ);
  ASSERT_ALWAYS (mpfr_cmp_si (x, -1) == 0);

  mpfr_set_ui (x, 2, GMP_RNDZ);
  d = mpfr_get_decimal64 (x, GMP_RNDZ);
  mpfr_set_ui (x, 0, GMP_RNDZ);
  mpfr_set_decimal64 (x, d, GMP_RNDZ);
  ASSERT_ALWAYS (mpfr_cmp_ui (x, 2) == 0);

  mpfr_set_ui (x, 99, GMP_RNDZ);
  d = mpfr_get_decimal64 (x, GMP_RNDZ);
  mpfr_set_ui (x, 0, GMP_RNDZ);
  mpfr_set_decimal64 (x, d, GMP_RNDZ);
  ASSERT_ALWAYS (mpfr_cmp_ui (x, 99) == 0);

  mpfr_set_str (x, "9999999999999999", 10, GMP_RNDZ);
  mpfr_set (y, x, GMP_RNDZ);
  d = mpfr_get_decimal64 (x, GMP_RNDZ);
  mpfr_set_ui (x, 0, GMP_RNDZ);
  mpfr_set_decimal64 (x, d, GMP_RNDZ);
  ASSERT_ALWAYS (mpfr_cmp (x, y) == 0);

  /* smallest normal number */
  mpfr_set_str (x, "1E-383", 10, GMP_RNDU);
  mpfr_set (y, x, GMP_RNDZ);
  d = mpfr_get_decimal64 (x, GMP_RNDZ);
  mpfr_set_ui (x, 0, GMP_RNDZ);
  mpfr_set_decimal64 (x, d, GMP_RNDU);
  ASSERT_ALWAYS (mpfr_cmp (x, y) == 0);

  /* smallest subnormal number */
  mpfr_set_str (x, "1E-398", 10, GMP_RNDU);
  mpfr_set (y, x, GMP_RNDZ);
  d = mpfr_get_decimal64 (x, GMP_RNDZ);
  mpfr_set_ui (x, 0, GMP_RNDZ);
  mpfr_set_decimal64 (x, d, GMP_RNDU);
  ASSERT_ALWAYS (mpfr_cmp (x, y) == 0);

  /* subnormal number with exponent change when we round back
     from 16 digits to 1 digit */
  mpfr_set_str (x, "9.9E-398", 10, GMP_RNDN);
  d = mpfr_get_decimal64 (x, GMP_RNDU); /* should be 1E-397 */
  mpfr_set_ui (x, 0, GMP_RNDZ);
  mpfr_set_decimal64 (x, d, GMP_RNDD);
  mpfr_set_str (y, "1E-397", 10, GMP_RNDN);
  ASSERT_ALWAYS (mpfr_cmp (x, y) == 0);

  /* largest number */
  mpfr_set_str (x, "9.999999999999999E384", 10, GMP_RNDZ);
  mpfr_set (y, x, GMP_RNDZ);
  d = mpfr_get_decimal64 (x, GMP_RNDU);
  mpfr_set_ui (x, 0, GMP_RNDZ);
  mpfr_set_decimal64 (x, d, GMP_RNDZ);
  ASSERT_ALWAYS (mpfr_cmp (x, y) == 0);

  mpfr_set_prec (x, 53);
  mpfr_set_prec (y, 53);

  /* largest number */
  mpfr_set_str (x, "9.999999999999999E384", 10, GMP_RNDZ);
  d = mpfr_get_decimal64 (x, GMP_RNDZ);
  mpfr_set_decimal64 (y, d, GMP_RNDU);
  ASSERT_ALWAYS (mpfr_cmp (x, y) == 0);

  mpfr_clear (x);
  mpfr_clear (y);
}

static void
check_random (void)
{
  mpfr_t  x, y;
  _Decimal64 d;
  int i;
  
  mpfr_init2 (x, 49);
  mpfr_init2 (y, 49);

  for (i = 0; i < 100000; i++)
    {
      mpfr_random (x); /* 0 <= x < 1 */
      /* the normal decimal64 range contains [2^(-1272), 2^1278] */
      mpfr_mul_2si (x, x, (i % 2550) - 1272, GMP_RNDN);
      if (mpfr_get_exp (x) <= -1272)
	mpfr_mul_2exp (x, x, -1271 - mpfr_get_exp (x), GMP_RNDN);
      d = mpfr_get_decimal64 (x, GMP_RNDN);
      mpfr_set_decimal64 (y, d, GMP_RNDN);
      if (mpfr_cmp (x, y) != 0)
	{
	  printf ("x="); mpfr_dump (x);
	  printf ("d="); print_decimal64 (d);
	  printf ("y="); mpfr_dump (y);
	  exit (1);
	}
    }

  mpfr_clear (x);
  mpfr_clear (y);
}

/* check with native decimal formats */
static void
check_native (void)
{
  mpfr_t x;
  _Decimal64 d;

  mpfr_init2 (x, 53);

  /* check important constants are correctly converted */
  mpfr_set_ui (x, 17, GMP_RNDN);
  d = mpfr_get_decimal64 (x, GMP_RNDN);
  MPFR_ASSERTN(d == 17.0dd);

  mpfr_set_ui (x, 42, GMP_RNDN);
  d = mpfr_get_decimal64 (x, GMP_RNDN);
  MPFR_ASSERTN(d == 42.0dd);

  mpfr_set_decimal64 (x, 17.0dd, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (x, 17) == 0);

  mpfr_set_decimal64 (x, 42.0dd, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (x, 42) == 0);

  mpfr_clear (x);
}
#endif /* MPFR_WANT_DECIMAL_FLOATS */

int
main (void)
{
  tests_start_mpfr ();
  mpfr_test_init ();

#if MPFR_WANT_DECIMAL_FLOATS
#ifdef DEBUG
#ifdef DPD_FORMAT
  printf ("Using DPD format\n");
#else
  printf ("Using BID format\n");
#endif
#endif
  check_inf_nan ();
  check_random ();
  check_native ();
#endif

  tests_end_mpfr ();
  return 0;
}
