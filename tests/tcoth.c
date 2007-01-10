/* Test file for mpfr_coth.

Copyright 2005, 2006, 2007 Free Software Foundation, Inc.

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

#define TEST_FUNCTION mpfr_coth
#include "tgeneric.c"

static void
check_specials (void)
{
  mpfr_t  x, y;

  mpfr_init2 (x, 123L);
  mpfr_init2 (y, 123L);

  mpfr_set_nan (x);
  mpfr_coth (y, x, GMP_RNDN);
  if (! mpfr_nan_p (y))
    {
      printf ("Error: coth(NaN) != NaN\n");
      exit (1);
    }

  mpfr_set_inf (x, 1);
  mpfr_coth (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 1))
    {
      printf ("Error: coth(Inf) != 1\n");
      exit (1);
    }

  mpfr_set_inf (x, -1);
  mpfr_coth (y, x, GMP_RNDN);
  if (mpfr_cmp_si (y, -1))
    {
      printf ("Error: coth(-Inf) != -1\n");
      exit (1);
    }

  /* cot(+/-0) = +/-0 */
  mpfr_set_ui (x, 0, GMP_RNDN);
  mpfr_coth (y, x, GMP_RNDN);
  if (! (mpfr_zero_p (y) && MPFR_SIGN (y) > 0))
    {
      printf ("Error: coth(+0) != +0\n");
      exit (1);
    }
  mpfr_neg (x, x, GMP_RNDN);
  mpfr_coth (y, x, GMP_RNDN);
  if (! (mpfr_zero_p (y) && MPFR_SIGN (y) < 0))
    {
      printf ("Error: coth(-0) != -0\n");
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);
}

static void
check_bugs (void)
{
  mpfr_t x, y;
  
  mpfr_init (x);
  mpfr_init (y);

  /* bug found by Rob (Sisyphus) on 16 Sep 2005 */
  mpfr_set_ui (x, 2, GMP_RNDN);
  mpfr_set_prec (y, 2);
  mpfr_coth (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 1))
    {
      printf ("Error for coth(2), expected 1, got ");
      mpfr_dump (y);
      exit (1);
    }

  mpfr_set_prec (x, 53);
  mpfr_set_prec (y, 53);

  mpfr_set_str (x, "18.368400284838550", 10, GMP_RNDN);
  mpfr_set_str (y, "1.0000000000000002", 10, GMP_RNDN);
  mpfr_coth (x, x, GMP_RNDN);
  if (mpfr_cmp (x, y) != 0)
    {
      printf ("Error for coth(18.368400284838550)\n");
      exit (1);
    }

  mpfr_set_str (x, "18.714973875118520", 10, GMP_RNDN);
  mpfr_coth (x, x, GMP_RNDN);
  if (mpfr_cmp (x, y) != 0)
    {
      printf ("Error for coth(18.714973875118520)\n");
      exit (1);
    }

  mpfr_set_str (x, "18.714973875118524", 10, GMP_RNDN);
  mpfr_coth (x, x, GMP_RNDN);
  if (mpfr_cmp_ui (x, 1) != 0)
    {
      printf ("Error for coth(18.714973875118524)\n");
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
  check_bugs ();
  test_generic (2, 200, 10);

  tests_end_mpfr ();
  return 0;
}
