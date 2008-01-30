/* Test file for mpfr_rec_sqrt.

Copyright 2008 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

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

#include <stdio.h>
#include <stdlib.h>

#include "mpfr-test.h"

#define TEST_FUNCTION mpfr_rec_sqrt
#define TEST_RANDOM_POS 8 /* 8/512 = 1/64 of the tested numbers are negative */
#include "tgeneric.c"

static void
special (void)
{
  mpfr_t x, y;
  int inex;

  mpfr_init (x);
  mpfr_init (y);

  /* rec_sqrt(NaN) = NaN */
  mpfr_set_nan (x);
  inex = mpfr_rec_sqrt (x, x, GMP_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (x) && inex == 0);

  /* rec_sqrt(+Inf) = +0 */
  mpfr_set_inf (x, 1);
  inex = mpfr_rec_sqrt (x, x, GMP_RNDN);
  MPFR_ASSERTN(mpfr_zero_p (x) && MPFR_IS_POS(x) && inex == 0);

  /* rec_sqrt(-Inf) = NaN */
  mpfr_set_inf (x, -1);
  inex = mpfr_rec_sqrt (x, x, GMP_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (x) && inex == 0);

  /* rec_sqrt(+0) = +Inf */
  mpfr_set_ui (x, 0, GMP_RNDN);
  inex = mpfr_rec_sqrt (x, x, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (x) && MPFR_IS_POS(x) && inex == 0);

  /* rec_sqrt(-0) = +Inf */
  mpfr_set_ui (x, 0, GMP_RNDN);
  mpfr_neg (x, x, GMP_RNDN);
  inex = mpfr_rec_sqrt (x, x, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (x) && MPFR_IS_POS(x) && inex == 0);

  /* rec_sqrt(-1) = NaN */
  mpfr_set_si (x, -1, GMP_RNDN);
  inex = mpfr_rec_sqrt (x, x, GMP_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (x) && inex == 0);

  /* rec_sqrt(1) = 1 */
  mpfr_set_ui (x, 1, GMP_RNDN);
  inex = mpfr_rec_sqrt (x, x, GMP_RNDN);
  MPFR_ASSERTN((mpfr_cmp_ui (x, 1) == 0) && (inex == 0));

  mpfr_set_prec (x, 23);
  mpfr_set_prec (y, 33);
  mpfr_set_str_binary (x, "1.0001110110101001010100e-1");
  inex = mpfr_rec_sqrt (y, x, GMP_RNDU);
  mpfr_set_prec (x, 33);
  mpfr_set_str_binary (x, "1.01010110101110100100100101011");
  MPFR_ASSERTN (inex > 0 && mpfr_cmp (x, y) == 0);

  mpfr_clear (x);
  mpfr_clear (y);
}

int
main (void)
{
  tests_start_mpfr ();

  special ();
  test_generic (2, 300, 15);

  tests_end_mpfr ();
  return 0;
}
