/* Test file for mpfr_sech.

Copyright 2005, 2006, 2007 Free Software Foundation, Inc.
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

#define TEST_FUNCTION mpfr_sech
#include "tgeneric.c"

static void
check_specials (void)
{
  mpfr_t  x, y;

  mpfr_init2 (x, 123L);
  mpfr_init2 (y, 123L);

  mpfr_set_nan (x);
  mpfr_sech (y, x, GMP_RNDN);
  if (! mpfr_nan_p (y))
    {
      printf ("Error: sech(NaN) != NaN\n");
      exit (1);
    }

  mpfr_set_inf (x, 1);
  mpfr_sech (y, x, GMP_RNDN);
  if (! (MPFR_IS_ZERO (y) && MPFR_SIGN (y) > 0))
    {
      printf ("Error: sech(+Inf) != +0\n");
      exit (1);
    }

  mpfr_set_inf (x, -1);
  mpfr_sech (y, x, GMP_RNDN);
  if (! (MPFR_IS_ZERO (y) && MPFR_SIGN (y) > 0))
    {
      printf ("Error: sech(-Inf) != +0\n");
      exit (1);
    }

  /* sec(+/-0) = 1 */
  mpfr_set_ui (x, 0, GMP_RNDN);
  mpfr_sech (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 1))
    {
      printf ("Error: sech(+0) != 1\n");
      exit (1);
    }
  mpfr_neg (x, x, GMP_RNDN);
  mpfr_sech (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 1))
    {
      printf ("Error: sech(-0) != 1\n");
      exit (1);
    }

  /* check huge x */
  mpfr_set_str (x, "8e8", 10, GMP_RNDN);
  mpfr_sech (y, x, GMP_RNDN);
  if (! (mpfr_zero_p (y) && MPFR_SIGN (y) > 0))
    {
      printf ("Error: sech(8e8) != +0\n");
      exit (1);
    }
  mpfr_set_str (x, "-8e8", 10, GMP_RNDN);
  mpfr_sech (y, x, GMP_RNDN);
  if (! (mpfr_zero_p (y) && MPFR_SIGN (y) > 0))
    {
      printf ("Error: sech(-8e8) != +0\n");
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);
}

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  check_specials ();
  test_generic (2, 200, 10);

  tests_end_mpfr ();
  return 0;
}
