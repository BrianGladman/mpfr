/* Test file for mpfr_total_order.

Copyright 2018 Free Software Foundation, Inc.
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
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include "mpfr-test.h"

int
main (void)
{
  mpfr_t x, y;

  tests_start_mpfr ();

  mpfr_init (x);
  mpfr_init (y);

  /* If x < y, totalOrder(x, y) is true */
  mpfr_set_ui (x, 1, MPFR_RNDN);
  mpfr_set_ui (y, 2, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_total_order (x, y) != 0);

  /* If x > y, totalOrder(x, y) is false */
  mpfr_set_ui (x, 2, MPFR_RNDN);
  mpfr_set_ui (y, 1, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_total_order (x, y) == 0);

  /* totalOrder(-0, +0) is true */
  mpfr_set_ui (x, 0, MPFR_RNDN);
  mpfr_neg (x, x, MPFR_RNDN);
  mpfr_set_ui (y, 0, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_total_order (x, y) != 0);

  /* totalOrder(+0, -0) is false */
  mpfr_set_ui (x, 0, MPFR_RNDN);
  mpfr_set_ui (y, 0, MPFR_RNDN);
  mpfr_neg (y, y, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_total_order (x, y) == 0);

  /* Case x = y */
  mpfr_set_si (x, -1, MPFR_RNDN);
  mpfr_set_si (y, -1, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_total_order (x, y) != 0);

  mpfr_set_si (x, 1, MPFR_RNDN);
  mpfr_set_si (y, 1, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_total_order (x, y) != 0);

  /* Case x = -NaN */
  mpfr_set_nan (x);
  /* warning: the sign of x is unspecified in mpfr_set_nan */
  if (!mpfr_signbit (x))
    mpfr_neg (x, x, MPFR_RNDN);
  /* now x = -NaN */
  mpfr_set_si (y, 1, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_total_order (x, y) != 0);
  mpfr_set_inf (y, -1);
  MPFR_ASSERTN(mpfr_total_order (x, y) != 0);
  mpfr_set_inf (y, +1);
  MPFR_ASSERTN(mpfr_total_order (x, y) != 0);

  /* Case y = +NaN */
  mpfr_set_nan (y);
  if (mpfr_signbit (y))
    mpfr_neg (y, y, MPFR_RNDN);
  /* now y = +NaN */
  mpfr_set_si (x, 1, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_total_order (x, y) != 0);
  mpfr_set_inf (x, -1);
  MPFR_ASSERTN(mpfr_total_order (x, y) != 0);
  mpfr_set_inf (x, +1);
  MPFR_ASSERTN(mpfr_total_order (x, y) != 0);

  /* Both x and y are NaN */
  mpfr_set_nan (x);
  if (!mpfr_signbit (x))
    mpfr_neg (x, x, MPFR_RNDN);
  /* now x = -NaN */
  mpfr_set_nan (y);
  if (mpfr_signbit (y))
    mpfr_neg (y, y, MPFR_RNDN);
  /* now y = +NaN */
  MPFR_ASSERTN(mpfr_total_order (x, y) != 0);
  mpfr_neg (y, y, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_total_order (x, y) != 0);

  mpfr_set_nan (x);
  if (mpfr_signbit (x))
    mpfr_neg (x, x, MPFR_RNDN);
  /* now x = +NaN */
  mpfr_set_nan (y);
  if (mpfr_signbit (y))
    mpfr_neg (y, y, MPFR_RNDN);
  /* now y = +NaN */
  MPFR_ASSERTN(mpfr_total_order (x, y) != 0);
  mpfr_neg (y, y, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_total_order (x, y) == 0);

  mpfr_clear (x);
  mpfr_clear (y);

  tests_end_mpfr ();
  return 0;
}
