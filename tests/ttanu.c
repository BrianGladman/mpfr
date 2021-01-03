/* Test file for mpfr_tanu.

Copyright 2020-2021 Free Software Foundation, Inc.
Contributed by the AriC and Caramba projects, INRIA.

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
https://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include "mpfr-test.h"

static void
test_singular (void)
{
  mpfr_t x, y;
  int inexact;

  mpfr_init (x);
  mpfr_init (y);

  /* check u = 0 */
  mpfr_set_ui (x, 17, MPFR_RNDN);
  inexact = mpfr_tanu (y, x, 0, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (y));

  /* check x = NaN */
  mpfr_set_nan (x);
  inexact = mpfr_tanu (y, x, 1, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (y));

  /* check x = +Inf */
  mpfr_set_inf (x, 1);
  inexact = mpfr_tanu (y, x, 1, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (y));

  /* check x = -Inf */
  mpfr_set_inf (x, -1);
  inexact = mpfr_tanu (y, x, 1, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (y));

  /* check x = +0 */
  mpfr_set_zero (x, 1);
  inexact = mpfr_tanu (y, x, 1, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_zero_p (y) && mpfr_signbit (y) == 0);
  MPFR_ASSERTN(inexact == 0);

  /* check x = -0 */
  mpfr_set_zero (x, -1);
  inexact = mpfr_tanu (y, x, 1, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_zero_p (y) && mpfr_signbit (y) != 0);
  MPFR_ASSERTN(inexact == 0);

  mpfr_clear (x);
  mpfr_clear (y);
}

static void
test_exact (void)
{
  mpfr_t x, y;
  int inexact, n;

  mpfr_init2 (x, 6);
  mpfr_init2 (y, 6);

  /* check n + 0.5 for n integer */
  for (n = 0; n < 10; n++)
    {
      /* check 2n+0.5 for n>=0: +Inf and divide by 0 exception */
      mpfr_set_ui (x, 4 * n + 1, MPFR_RNDN);
      mpfr_clear_divby0 ();
      inexact = mpfr_tanu (y, x, 4, MPFR_RNDN);
      MPFR_ASSERTN(mpfr_inf_p (y) && mpfr_sgn (y) > 0);
      MPFR_ASSERTN(inexact == 0);
      MPFR_ASSERTN(mpfr_divby0_p ());

      /* check 2n+1 for n>=0: -0 */
      mpfr_set_ui (x, 4 * n + 2, MPFR_RNDN);
      mpfr_clear_divby0 ();
      inexact = mpfr_tanu (y, x, 4, MPFR_RNDN);
      MPFR_ASSERTN(mpfr_zero_p (y) && mpfr_signbit (y) != 0);
      MPFR_ASSERTN(inexact == 0);
      MPFR_ASSERTN(!mpfr_divby0_p ());

      /* check 2n+1.5 for n>=0: -Inf and divide by 0 exception */
      mpfr_set_ui (x, 4 * n + 3, MPFR_RNDN);
      mpfr_clear_divby0 ();
      inexact = mpfr_tanu (y, x, 4, MPFR_RNDN);
      MPFR_ASSERTN(mpfr_inf_p (y) && mpfr_sgn (y) < 0);
      MPFR_ASSERTN(inexact == 0);
      MPFR_ASSERTN(mpfr_divby0_p ());

      /* check 2n+2 for n>=0: +0 */
      mpfr_set_ui (x, 4 * n + 4, MPFR_RNDN);
      mpfr_clear_divby0 ();
      inexact = mpfr_tanu (y, x, 4, MPFR_RNDN);
      MPFR_ASSERTN(mpfr_zero_p (y) && mpfr_signbit (y) == 0);
      MPFR_ASSERTN(inexact == 0);
      MPFR_ASSERTN(!mpfr_divby0_p ());

      /* check -2n-0.5 for n>=0: -Inf and divide by 0 exception */
      mpfr_set_si (x, -4 * n - 1, MPFR_RNDN);
      mpfr_clear_divby0 ();
      inexact = mpfr_tanu (y, x, 4, MPFR_RNDN);
      MPFR_ASSERTN(mpfr_inf_p (y) && mpfr_sgn (y) < 0);
      MPFR_ASSERTN(inexact == 0);
      MPFR_ASSERTN(mpfr_divby0_p ());

      /* check -2n-1 for n>=0: +0 */
      mpfr_set_si (x, -4 * n - 2, MPFR_RNDN);
      mpfr_clear_divby0 ();
      inexact = mpfr_tanu (y, x, 4, MPFR_RNDN);
      MPFR_ASSERTN(mpfr_zero_p (y) && mpfr_signbit (y) == 0);
      MPFR_ASSERTN(inexact == 0);
      MPFR_ASSERTN(!mpfr_divby0_p ());

      /* check -2n-1.5 for n>=0: +Inf and divide by 0 exception */
      mpfr_set_si (x, -4 * n - 3, MPFR_RNDN);
      mpfr_clear_divby0 ();
      inexact = mpfr_tanu (y, x, 4, MPFR_RNDN);
      MPFR_ASSERTN(mpfr_inf_p (y) && mpfr_sgn (y) > 0);
      MPFR_ASSERTN(inexact == 0);
      MPFR_ASSERTN(mpfr_divby0_p ());

      /* check -2n-2 for n>=0: -0 */
      mpfr_set_si (x, -4 * n - 4, MPFR_RNDN);
      mpfr_clear_divby0 ();
      inexact = mpfr_tanu (y, x, 4, MPFR_RNDN);
      MPFR_ASSERTN(mpfr_zero_p (y) && mpfr_signbit (y) != 0);
      MPFR_ASSERTN(inexact == 0);
      MPFR_ASSERTN(!mpfr_divby0_p ());
    }

  /* check 2*pi*x/u = pi/4 thus x/u = 1/8, for example x=1 and u=8 */
  mpfr_set_ui (x, 1, MPFR_RNDN);
  inexact = mpfr_tanu (y, x, 8, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (y, 1) == 0 && inexact == 0);

  /* check 2*pi*x/u = 3*pi/4 thus x/u = 3/8, for example x=3 and u=8 */
  mpfr_set_ui (x, 3, MPFR_RNDN);
  inexact = mpfr_tanu (y, x, 8, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_cmp_si (y, -1) == 0 && inexact == 0);

  /* check 2*pi*x/u = 5*pi/4 thus x/u = 5/8, for example x=5 and u=8 */
  mpfr_set_ui (x, 5, MPFR_RNDN);
  inexact = mpfr_tanu (y, x, 8, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (y, 1) == 0 && inexact == 0);

  /* check 2*pi*x/u = 7*pi/4 thus x/u = 7/8, for example x=7 and u=8 */
  mpfr_set_ui (x, 7, MPFR_RNDN);
  inexact = mpfr_tanu (y, x, 8, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_cmp_si (y, -1) == 0 && inexact == 0);

  mpfr_clear (x);
  mpfr_clear (y);
}

static void
test_regular (void)
{
  mpfr_t x, y, z;
  int inexact;

  mpfr_init2 (x, 53);
  mpfr_init2 (y, 53);
  mpfr_init2 (z, 53);

  mpfr_set_ui (x, 17, MPFR_RNDN);
  inexact = mpfr_tanu (y, x, 42, MPFR_RNDN);
  /* y should be tan(2*17*pi/42) rounded to nearest */
  mpfr_set_str (z, "-0xa.e89b03074638p-4", 16, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_equal_p (y, z));
  MPFR_ASSERTN(inexact > 0);

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
}

/* FIXME[VL]: For mpfr_tanu, the range reduction should not be expensive.
   If I'm not mistaken, this is linear in the bitsize of the exponent
   since one just needs to compute the argument modulo the integer u. */
#define TEST_FUNCTION mpfr_tanu
#define ULONG_ARG2
#ifndef MPFR_USE_MINI_GMP
#define REDUCE_EMAX 262143 /* otherwise arg. reduction is too expensive */
#else
#define REDUCE_EMAX 16383  /* reduce further since mini-gmp works in O(n^2) */
#endif
#include "tgeneric.c"

static int
mpfr_tan2pi (mpfr_ptr y, mpfr_srcptr x, mpfr_rnd_t r)
{
  return mpfr_tanu (y, x, 1, r);
}

int
main (void)
{
  tests_start_mpfr ();

  test_singular ();
  test_exact ();
  test_regular ();

  test_generic (MPFR_PREC_MIN, 100, 1);

  data_check ("data/tan2pi", mpfr_tan2pi, "mpfr_tan2pi");

  tests_end_mpfr ();
  return 0;
}
