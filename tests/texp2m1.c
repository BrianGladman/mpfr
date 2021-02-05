/* Test file for mpfr_exp2m1.

Copyright 2001-2021 Free Software Foundation, Inc.
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

#define TEST_FUNCTION mpfr_exp2m1
#define TEST_RANDOM_EMIN -36
#define TEST_RANDOM_EMAX 36
#include "tgeneric.c"

static void
special (void)
{
  mpfr_t x, y;
  int i;

  mpfr_init2 (x, 2);
  mpfr_init2 (y, 2);

  mpfr_set_nan (x);
  mpfr_exp2m1 (y, x, MPFR_RNDN);
  if (!mpfr_nan_p (y))
    {
      printf ("Error for exp2m1(NaN)\n");
      exit (1);
    }

  mpfr_set_inf (x, 1);
  mpfr_exp2m1 (y, x, MPFR_RNDN);
  if (!mpfr_inf_p (y) || mpfr_sgn (y) < 0)
    {
      printf ("Error for exp2m1(+Inf)\n");
      exit (1);
    }

  mpfr_set_inf (x, -1);
  mpfr_exp2m1 (y, x, MPFR_RNDN);
  if (mpfr_cmp_si (y, -1) != 0)
    {
      printf ("Error for exp2m1(-Inf)\n");
      exit (1);
    }

  mpfr_set_ui (x, 0, MPFR_RNDN);
  mpfr_exp2m1 (y, x, MPFR_RNDN);
  if (MPFR_NOTZERO (y) || MPFR_IS_NEG (y))
    {
      printf ("Error for exp2m1(+0)\n");
      exit (1);
    }

  mpfr_neg (x, x, MPFR_RNDN);
  mpfr_exp2m1 (y, x, MPFR_RNDN);
  if (MPFR_NOTZERO (y) || MPFR_IS_POS (y))
    {
      printf ("Error for exp2m1(-0)\n");
      exit (1);
    }

  /* Check overflow of exp2m1(x) */
  mpfr_clear_flags ();
  mpfr_set_str_binary (x, "1.1E1000000000");
  i = mpfr_exp2m1 (x, x, MPFR_RNDN);
  MPFR_ASSERTN (MPFR_IS_INF (x) && MPFR_IS_POS (x));
  MPFR_ASSERTN (mpfr_overflow_p ());
  MPFR_ASSERTN (i > 0);

  mpfr_clear_flags ();
  mpfr_set_str_binary (x, "1.1E1000000000");
  i = mpfr_exp2m1 (x, x, MPFR_RNDU);
  MPFR_ASSERTN (MPFR_IS_INF (x) && MPFR_IS_POS (x));
  MPFR_ASSERTN (mpfr_overflow_p ());
  MPFR_ASSERTN (i > 0);

  mpfr_clear_flags ();
  mpfr_set_str_binary (x, "1.1E1000000000");
  i = mpfr_exp2m1 (x, x, MPFR_RNDD);
  MPFR_ASSERTN (!MPFR_IS_INF (x) && MPFR_IS_POS (x));
  MPFR_ASSERTN (mpfr_overflow_p ());
  MPFR_ASSERTN (i < 0);

  /* Check internal underflow of exp2m1 (x) */
  mpfr_set_prec (x, 2);
  mpfr_clear_flags ();
  mpfr_set_str_binary (x, "-1.1E1000000000");
  i = mpfr_exp2m1 (x, x, MPFR_RNDN);
  MPFR_ASSERTN (mpfr_cmp_si (x, -1) == 0);
  MPFR_ASSERTN (!mpfr_overflow_p () && !mpfr_underflow_p ());
  MPFR_ASSERTN (i < 0);

  mpfr_set_str_binary (x, "-1.1E1000000000");
  i = mpfr_exp2m1 (x, x, MPFR_RNDD);
  MPFR_ASSERTN (mpfr_cmp_si (x, -1) == 0);
  MPFR_ASSERTN (!mpfr_overflow_p () && !mpfr_underflow_p ());
  MPFR_ASSERTN (i < 0);

  mpfr_set_str_binary (x, "-1.1E1000000000");
  i = mpfr_exp2m1 (x, x, MPFR_RNDZ);
  MPFR_ASSERTN (mpfr_cmp_str (x, "-0.11", 2, MPFR_RNDN) == 0);
  MPFR_ASSERTN (!mpfr_overflow_p () && !mpfr_underflow_p ());
  MPFR_ASSERTN (i > 0);

  mpfr_set_str_binary (x, "-1.1E1000000000");
  i = mpfr_exp2m1 (x, x, MPFR_RNDU);
  MPFR_ASSERTN (mpfr_cmp_str (x, "-0.11", 2, MPFR_RNDN) == 0);
  MPFR_ASSERTN (!mpfr_overflow_p () && !mpfr_underflow_p ());
  MPFR_ASSERTN (i > 0);

  /* hard-coded test */
  mpfr_set_prec (x, 32);
  mpfr_set_prec (y, 32);
  mpfr_set_ui_2exp (x, 1686629713UL, -29, MPFR_RNDN); /* approximates Pi */
  i = mpfr_exp2m1 (y, x, MPFR_RNDN);
  MPFR_ASSERTN (mpfr_cmp_ui_2exp (y, 2100501491UL, -28) == 0);
  MPFR_ASSERTN (i < 0);

  mpfr_clear (x);
  mpfr_clear (y);
}

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  special ();

  test_generic (MPFR_PREC_MIN, 100, 100);

  tests_end_mpfr ();
  return 0;
}
