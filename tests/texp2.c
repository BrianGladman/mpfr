/* Test file for mpfr_exp2.

Copyright 2001, 2002, 2003, 2004 Free Software Foundation.
Adapted from tarctan.c.

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
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include <stdio.h>
#include <stdlib.h>

#include "mpfr-test.h"

#define TEST_FUNCTION mpfr_exp2
#include "tgeneric.c"

int
main (int argc, char *argv[])
{
  mpfr_t x, y;

  tests_start_mpfr ();

  mpfr_init (x);
  mpfr_init (y);

  mpfr_set_ui (x, 4, GMP_RNDN);
  mpfr_exp2 (y, x, GMP_RNDN);
  mpfr_exp2 (y, x, GMP_RNDD);
  mpfr_exp2 (y, x, GMP_RNDU);
  if (mpfr_cmp_ui (y, 16) != 0)
    {
      printf ("Error for 2^4\n");
      exit (1);
    }

  mpfr_set_si (x, -4, GMP_RNDN);
  mpfr_exp2 (y, x, GMP_RNDN);
  mpfr_exp2 (y, x, GMP_RNDD);
  mpfr_exp2 (y, x, GMP_RNDU);
  if (mpfr_cmp_ui_2exp (y, 1, -4) != 0)
    {
      printf ("Error for 2^(-4)\n");
      exit (1);
    }

  mpfr_set_prec (x, 53);
  mpfr_set_prec (y, 53);
  mpfr_set_str (x, /*-1683977482443233.0 / 2199023255552.0*/
		"-7.6578429909351734750089235603809357e2", 10,GMP_RNDN);
  mpfr_exp2 (y, x, GMP_RNDN);
  if (mpfr_cmp_str1 (y, "2.991959870867646566478e-231"))
    {
      printf ("Error for x=-1683977482443233/2^41\n");
      exit (1);
    }

  MPFR_SET_INF(x);
  MPFR_SET_POS(x);
  mpfr_exp2 (y, x, GMP_RNDN);
  if(!MPFR_IS_INF(y))
    {
      printf ("evaluation of function in INF does not return INF\n");
      exit (1);
    }

  MPFR_CHANGE_SIGN(x);
  mpfr_exp2 (y, x, GMP_RNDN);
  if(!MPFR_IS_ZERO(y))
    {
      printf ("evaluation of function in -INF does not return 0\n");
      exit (1);
    }

  MPFR_SET_NAN(x);
  mpfr_exp2 (y, x, GMP_RNDN);
  if(!MPFR_IS_NAN(y))
    {
      printf ("evaluation of function in NaN does not return NaN\n");
      exit (1);
    }

  test_generic (2, 100, 100);

  mpfr_clear (x);
  mpfr_clear (y);

  tests_end_mpfr ();
  return 0;
}
