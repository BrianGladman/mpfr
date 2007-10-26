/* mpfr_tli2 -- test file for dilogarithm function

Copyright 2007 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

This file is part of the MPFR Library, and was contributed by Mathieu Dutour.

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

#define TEST_FUNCTION mpfr_li2
#include "tgeneric.c"

static void
special (void)
{
  mpfr_t x, y;
  mpfr_init (x);
  mpfr_init (y);

  mpfr_set_nan (x);
  mpfr_li2 (y, x, GMP_RNDN);
  if (!mpfr_nan_p (y))
    {
      printf ("Error for li2(NaN)\n");
      exit (1);
    }

  mpfr_set_inf (x, -1);
  mpfr_li2 (y, x, GMP_RNDN);
  if (!mpfr_inf_p (y) || mpfr_sgn (y) > 0)
    {
      printf ("Error for li2(-Inf)\n");
      exit (1);
    }

  mpfr_set_inf (x, 1);
  mpfr_li2 (y, x, GMP_RNDN);
  if (!mpfr_inf_p (y) || mpfr_sgn (y) > 0)
    {
      printf ("Error for li2(+Inf)\n");
      exit (1);
    }

  mpfr_set_ui (x, 0, GMP_RNDN);
  mpfr_li2 (y, x, GMP_RNDN);
  if (!mpfr_zero_p (y) || mpfr_sgn (y) < 0)
    {
      printf ("Error for li2(+0)\n");
      exit (1);
    }

  mpfr_set_ui (x, 0, GMP_RNDN);
  mpfr_neg (x, x, GMP_RNDN);
  mpfr_li2 (y, x, GMP_RNDN);
  if (!mpfr_zero_p (y) || mpfr_sgn (y) > 0)
    {
      printf ("Error for li2(-0)\n");
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);
}

static void
normal (void)
{
  mpfr_t x, y;
  mpfr_init (x);
  mpfr_init (y);

  /* x1 = 2^-3 */
  mpfr_set_str (x, "1e-2", 2, GMP_RNDD);
  mpfr_set_str (y, "1.0000100001111010011110101001111001000010000101000001e-3",
		2, GMP_RNDD);
  mpfr_li2 (x, x, GMP_RNDN);
  if (!mpfr_equal_p (x, y))
    {
      printf ("Error for li2(x1)\n");
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);
}

int
main (int argc, char *argv[])
{
  MPFR_TEST_USE_RANDS ();
  tests_start_mpfr ();

  special ();

  test_generic (2, 100, 2);

  data_check ("data/li2", mpfr_li2, "mpfr_li2");

  tests_end_mpfr ();
  return 0;
}
