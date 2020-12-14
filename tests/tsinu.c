/* Test file for mpfr_sinu.

Copyright 2020 Free Software Foundation, Inc.
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
  inexact = mpfr_sinu (y, x, 0, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (y));

  /* check x = NaN */
  mpfr_set_nan (x);
  inexact = mpfr_sinu (y, x, 1, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (y));

  /* check x = +Inf */
  mpfr_set_inf (x, 1);
  inexact = mpfr_sinu (y, x, 1, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (y));

  /* check x = -Inf */
  mpfr_set_inf (x, -1);
  inexact = mpfr_sinu (y, x, 1, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (y));

  /* check x = +0 */
  mpfr_set_zero (x, 1);
  inexact = mpfr_sinu (y, x, 1, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_zero_p (y) && mpfr_signbit (y) == 0);
  MPFR_ASSERTN(inexact == 0);

  /* check x = -0 */
  mpfr_set_zero (x, -1);
  inexact = mpfr_sinu (y, x, 1, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_zero_p (y) && mpfr_signbit (y) != 0);
  MPFR_ASSERTN(inexact == 0);

  mpfr_clear (x);
  mpfr_clear (y);
}

static void
test_exact (void)
{
  mpfr_t x, y;
  int inexact;

  mpfr_init (x);
  mpfr_init (y);

  /* check 2*pi*x/u = pi/2 thus x/u = pi/4 for x=1 and u=4 */
  mpfr_set_ui (x, 1, MPFR_RNDN);
  inexact = mpfr_sinu (y, x, 4, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (y, 1) == 0 && inexact == 0);

  /* check 2*pi*x/u = pi thus x/u=pi/2 for x=2 and u=4 */
  mpfr_set_ui (x, 2, MPFR_RNDN);
  inexact = mpfr_sinu (y, x, 4, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_zero_p (y) && mpfr_signbit (y) == 0);
  MPFR_ASSERTN(inexact == 0);

  /* check 2*pi*x/u = 3*pi/2 thus x/u = 3*pi/4 for x=3 and u=4 */
  mpfr_set_ui (x, 3, MPFR_RNDN);
  inexact = mpfr_sinu (y, x, 4, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_cmp_si (y, -1) == 0 && inexact == 0);

  /* check 2*pi*x/u = 2*pi thus x/u = pi for x=4 and u=4 */
  mpfr_set_ui (x, 4, MPFR_RNDN);
  inexact = mpfr_sinu (y, x, 4, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_zero_p (y) && mpfr_signbit (y) == 0);
  MPFR_ASSERTN(inexact == 0);

  /* check 2*pi*x/u = -pi/2 thus x/u = pi/4 for x=-1 and u=4 */
  mpfr_set_si (x, -1, MPFR_RNDN);
  inexact = mpfr_sinu (y, x, 4, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_cmp_si (y, -1) == 0 && inexact == 0);

  /* check 2*pi*x/u = -pi thus x/u=pi/2 for x=-2 and u=4 */
  mpfr_set_si (x, -2, MPFR_RNDN);
  inexact = mpfr_sinu (y, x, 4, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_zero_p (y) && mpfr_signbit (y) != 0);
  MPFR_ASSERTN(inexact == 0);

  /* check 2*pi*x/u = -3*pi/2 thus x/u = -3*pi/4 for x=3 and u=4 */
  mpfr_set_si (x, -3, MPFR_RNDN);
  inexact = mpfr_sinu (y, x, 4, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (y, 1) == 0 && inexact == 0);

  /* check 2*pi*x/u = -2*pi thus x/u = pi for x=4 and u=4 */
  mpfr_set_si (x, -4, MPFR_RNDN);
  inexact = mpfr_sinu (y, x, 4, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_zero_p (y) && mpfr_signbit (y) != 0);
  MPFR_ASSERTN(inexact == 0);

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
  inexact = mpfr_sinu (y, x, 42, MPFR_RNDN);
  /* y should be sin(2*17*pi/42) rounded to nearest */
  mpfr_set_str (z, "0x9.035be4a906768p-4", 16, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_equal_p (y, z));
  MPFR_ASSERTN(inexact > 0);

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
}

/* FIXME[VL]: For mpfr_sinu, the range reduction should not be expensive.
   If I'm not mistaken, this is linear in the bitsize of the exponent
   since one just needs to compute the argument modulo the integer u. */
#define TEST_FUNCTION mpfr_sinu
#define ULONG_ARG2
#ifndef MPFR_USE_MINI_GMP
#define REDUCE_EMAX 262143 /* otherwise arg. reduction is too expensive */
#else
#define REDUCE_EMAX 16383  /* reduce further since mini-gmp works in O(n^2) */
#endif
#include "tgeneric.c"

int
main (void)
{
  tests_start_mpfr ();

  test_singular ();
  test_exact ();
  test_regular ();

  test_generic (MPFR_PREC_MIN, 100, 1);

  tests_end_mpfr ();
  return 0;
}
