/* Test file for mpfr_cosh.

Copyright 2001, 2002, 2004 Free Software Foundation.
Adapted from tatan.c.

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

#define TEST_FUNCTION mpfr_cosh
#include "tgeneric.c"

static void
special (void)
{
  mpfr_t  x, y;

  mpfr_init (x);
  mpfr_init (y);

  mpfr_set_nan (x);
  mpfr_cosh (y, x, GMP_RNDN);
  if (!mpfr_nan_p (y))
    {
      printf ("Error: cosh(NaN) != NaN\n");
      exit (1);
    }

  mpfr_set_inf (x, 1);
  mpfr_cosh (y, x, GMP_RNDN);
  if (!mpfr_inf_p (y) || mpfr_sgn (y) < 0)
    {
      printf ("Error: cosh(+Inf) != +Inf\n");
      exit (1);
    }

  mpfr_set_inf (x, -1);
  mpfr_cosh (y, x, GMP_RNDN);
  if (!mpfr_inf_p (y) || mpfr_sgn (y) < 0)
    {
      printf ("Error: cosh(-Inf) != +Inf\n");
      exit (1);
    }

  /* cosh(+/-0) = 1 */
  mpfr_set_ui (x, 0, GMP_RNDN);
  mpfr_cosh (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 1))
    {
      printf ("Error: cosh(+0) != 1\n");
      exit (1);
    }
  mpfr_neg (x, x, GMP_RNDN);
  mpfr_cosh (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 1))
    {
      printf ("Error: cosh(-0) != 1\n");
      exit (1);
    }

  mpfr_set_prec (x, 32);
  mpfr_set_prec (y, 32);

  mpfr_set_str_binary (x, "0.1101110111111111001011101000101");
  mpfr_set_str_binary (y, "1.0110011001110000101100011001001");
  mpfr_cosh (x, x, GMP_RNDN);
  if (mpfr_cmp (x, y))
    {
      printf ("Error: mpfr_cosh for prec=32 (1)\n");
      exit (1);
    }

  mpfr_set_str_binary (x, "-0.1110111000011101010111100000101E-1");
  mpfr_set_str_binary (y, "1.0001110000101111111111100110101");
  mpfr_cosh (x, x, GMP_RNDN);
  if (mpfr_cmp (x, y))
    {
      printf ("Error: mpfr_cosh for prec=32 (2)\n");
      exit (1);
    }

  

  mpfr_clear (x);
  mpfr_clear (y);
}

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  special ();

  test_generic (2, 100, 100);

  tests_end_mpfr ();
  return 0;
}
