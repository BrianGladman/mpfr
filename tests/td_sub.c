/* Test file for mpfr_d_sub

Copyright 2007 Free Software Foundation, Inc.
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
#include <float.h>

#include "mpfr-test.h"

static void
check_nans (void)
{
  mpfr_t  x, y;

  mpfr_init2 (x, 123);
  mpfr_init2 (y, 123);

  /* 1.0 - nan is nan */
  mpfr_set_nan (x);
  mpfr_d_sub (y, 1.0, x, GMP_RNDN);
  MPFR_ASSERTN (mpfr_nan_p (y));

  /* 1.0 - +inf == -inf */
  mpfr_set_inf (x, 1);
  mpfr_d_sub (y, 1.0, x, GMP_RNDN);
  MPFR_ASSERTN (mpfr_inf_p (y));
  MPFR_ASSERTN (mpfr_sgn (y) < 0);

  /* 1.0 - -inf == +inf */
  mpfr_set_inf (x, -1);
  mpfr_d_sub (y, 1.0, x, GMP_RNDN);
  MPFR_ASSERTN (mpfr_inf_p (y));
  MPFR_ASSERTN (mpfr_sgn (y) > 0);

  mpfr_clear (x);
  mpfr_clear (y);
}

#define TEST_FUNCTION mpfr_d_sub
#define DOUBLE_ARG1
#define RAND_FUNCTION(x) mpfr_random2(x, MPFR_LIMB_SIZE (x), 1)
#include "tgeneric.c"

int
main (int argc, char *argv[])
{
  mpfr_t x, y, z;
  double d;
  int inexact;

  MPFR_TEST_USE_RANDS ();
  tests_start_mpfr ();

  /* check with enough precision */
  mpfr_init2 (x, IEEE_DBL_MANT_DIG);
  mpfr_init2 (y, IEEE_DBL_MANT_DIG);
  mpfr_init2 (z, IEEE_DBL_MANT_DIG);
  
  mpfr_set_str (y, "4096", 10, GMP_RNDN);
  d = 0.125;
  mpfr_clear_flags ();
  inexact = mpfr_d_sub (x, d, y, GMP_RNDN);
  if (inexact != 0)
    {
      printf ("Inexact flag error in mpfr_d_sub\n");
      exit (1);
    }
  mpfr_set_str (z, "-4095.875", 10, GMP_RNDN);
  if (mpfr_cmp (z, x))
    {
      printf ("Error in mpfr_d_sub (");
      mpfr_out_str (stdout, 10, 7, y, GMP_RNDN);
      printf (" + %.20g)\nexpected ", d);
      mpfr_out_str (stdout, 10, 0, z, GMP_RNDN);
      printf ("\ngot     ");
      mpfr_out_str (stdout, 10, 0, x, GMP_RNDN);
      printf ("\n");
      exit (1);
    }
  mpfr_clears (x, y, z, (void *)0);

  check_nans ();

  test_generic (2, 1000, 100);

  tests_end_mpfr ();
  return 0;
}
