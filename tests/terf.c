/* Test file for mpfr_erf.

Copyright 2001, 2002, 2003, 2004 Free Software Foundation, Inc.
Contributed by Ludovic Meunier and Paul Zimmermann.

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mpfr-test.h"

#define TEST_FUNCTION mpfr_erf
#define test_generic test_generic_erf
#include "tgeneric.c"

#define TEST_FUNCTION mpfr_erfc
#define test_generic test_generic_erfc
#include "tgeneric.c"

static void
special_erf (void)
{
  mpfr_t x, y;
  int inex;

  mpfr_init2 (x, 53);
  mpfr_init2 (y, 53);

  /* erf(NaN) = NaN */
  mpfr_set_nan (x);
  mpfr_erf (y, x, GMP_RNDN);
  if (!mpfr_nan_p (y))
    {
      printf ("mpfr_erf failed for x=NaN\n");
      exit (1);
    }

  /* erf(+Inf) = 1 */
  mpfr_set_inf (x, 1);
  mpfr_erf (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 1))
    {
      printf ("mpfr_erf failed for x=+Inf\n");
      printf ("expected 1.0, got ");
      mpfr_out_str (stdout, 2, 0, y, GMP_RNDN);
      printf ("\n");
      exit (1);
    }

  /* erf(-Inf) = -1 */
  mpfr_set_inf (x, -1);
  mpfr_erf (y, x, GMP_RNDN);
  if (mpfr_cmp_si (y, -1))
    {
      printf ("mpfr_erf failed for x=-Inf\n");
      exit (1);
    }

  /* erf(+0) = +0 */
  mpfr_set_ui (x, 0, GMP_RNDN);
  mpfr_erf (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 0) || mpfr_sgn (y) < 0)
    {
      printf ("mpfr_erf failed for x=+0\n");
      exit (1);
    }

  /* erf(-0) = -0 */
  mpfr_neg (x, x, GMP_RNDN);
  mpfr_erf (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 0) || mpfr_sgn (y) > 0)
    {
      printf ("mpfr_erf failed for x=-0\n");
      exit (1);
    }

  mpfr_set_ui (x, 1, GMP_RNDN);
  mpfr_erf (x, x, GMP_RNDN);
  mpfr_set_str_binary (y, "0.11010111101110110011110100111010000010000100010001011");
  if (mpfr_cmp (x, y))
    {
      printf ("mpfr_erf failed for x=1.0, rnd=GMP_RNDN\n");
      printf ("expected ");
      mpfr_out_str (stdout, 2, 0, y, GMP_RNDN);
      printf ("\n");
      printf ("got      ");
      mpfr_out_str (stdout, 2, 0, x, GMP_RNDN);
      printf ("\n");
      exit (1);
    }

  mpfr_set_str (x, "6.6", 10, GMP_RNDN);
  mpfr_erf (x, x, GMP_RNDN);
  if (mpfr_cmp_ui (x, 1))
    {
      printf ("mpfr_erf failed for x=6.6, rnd=GMP_RNDN\n");
      printf ("expected 1\n");
      printf ("got      ");
      mpfr_out_str (stdout, 2, 0, x, GMP_RNDN);
      printf ("\n");
      exit (1);
    }

  mpfr_set_str (x, "-6.6", 10, GMP_RNDN);
  mpfr_erf (x, x, GMP_RNDN);
  if (mpfr_cmp_si (x, -1))
    {
      printf ("mpfr_erf failed for x=-6.6, rnd=GMP_RNDN\n");
      printf ("expected -1\n");
      printf ("got      ");
      mpfr_out_str (stdout, 2, 0, x, GMP_RNDN);
      printf ("\n");
      exit (1);
    }

  mpfr_set_str (x, "6.6", 10, GMP_RNDN);
  mpfr_erf (x, x, GMP_RNDZ);
  mpfr_set_str_binary (y, "0.11111111111111111111111111111111111111111111111111111");
  if (mpfr_cmp (x, y))
    {
      printf ("mpfr_erf failed for x=6.6, rnd=GMP_RNDZ\n");
      printf ("expected ");
      mpfr_out_str (stdout, 2, 0, y, GMP_RNDN);
      printf ("\n");
      printf ("got      ");
      mpfr_out_str (stdout, 2, 0, x, GMP_RNDN);
      printf ("\n");
      exit (1);
    }

  mpfr_set_str (x, "4.5", 10, GMP_RNDN);
  mpfr_erf (x, x, GMP_RNDN);
  mpfr_set_str_binary (y, "0.1111111111111111111111111111111100100111110100011");
  if (mpfr_cmp (x, y))
    {
      printf ("mpfr_erf failed for x=4.5, rnd=GMP_RNDN\n");
      printf ("expected ");
      mpfr_out_str (stdout, 2, 0, y, GMP_RNDN);
      printf ("\n");
      printf ("got      ");
      mpfr_out_str (stdout, 2, 0, x, GMP_RNDN);
      printf ("\n");
      exit (1);
    }

  mpfr_set_prec (x, 120);
  mpfr_set_prec (y, 120);
  mpfr_set_str_binary (x, "0.110100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011E3");
  mpfr_erf (x, x, GMP_RNDN);
  mpfr_set_str_binary (y, "0.11111111111111111111111111111111111111111111111111111111111111111100111111000100111011111011010000110101111100011001101");
  if (mpfr_cmp (x, y))
    {
      printf ("mpfr_erf failed for x=6.6, rnd=GMP_RNDN\n");
      printf ("expected ");
      mpfr_out_str (stdout, 2, 0, y, GMP_RNDN);
      printf ("\n");
      printf ("got      ");
      mpfr_out_str (stdout, 2, 0, x, GMP_RNDN);
      printf ("\n");
      exit (1);
    }

  mpfr_set_prec (x, 8);
  mpfr_set_prec (y, 8);
  mpfr_set_ui (x, 50, GMP_RNDN);
  inex = mpfr_erf (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 1))
    {
      printf ("mpfr_erf failed for x=50, rnd=GMP_RNDN\n");
      printf ("expected 1, got ");
      mpfr_out_str (stdout, 2, 0, y, GMP_RNDN);
      printf ("\n");
      exit (1);
    }
  if (inex <= 0)
    {
      printf ("mpfr_erf failed for x=50, rnd=GMP_RNDN: wrong ternary value\n"
              "expected positive, got %d\n", inex);
      exit (1);
    }
  inex = mpfr_erf (x, x, GMP_RNDZ);
  mpfr_nextbelow (y);
  if (mpfr_cmp (x, y))
    {
      printf ("mpfr_erf failed for x=50, rnd=GMP_RNDZ\n");
      printf ("expected ");
      mpfr_out_str (stdout, 2, 0, y, GMP_RNDN);
      printf ("\n");
      printf ("got      ");
      mpfr_out_str (stdout, 2, 0, x, GMP_RNDN);
      printf ("\n");
      exit (1);
    }
  if (inex >= 0)
    {
      printf ("mpfr_erf failed for x=50, rnd=GMP_RNDN: wrong ternary value\n"
              "expected negative, got %d\n", inex);
      exit (1);
    }

  mpfr_set_prec (x, 32);
  mpfr_set_prec (y, 32);

  mpfr_set_str_binary (x, "0.1010100100111011001111100101E-1");
  mpfr_set_str_binary (y, "0.10111000001110011010110001101011E-1");
  mpfr_erf (x, x, GMP_RNDN);
  if (mpfr_cmp (x, y))
    {
      printf ("Error: erf for prec=32 (1)\n");
      exit (1);
    }

  mpfr_set_str_binary (x, "-0.10110011011010111110010001100001");
  mpfr_set_str_binary (y, "-0.1010110110101011100010111000111");
  mpfr_erf (x, x, GMP_RNDN);
  if (mpfr_cmp (x, y))
    {
      printf ("Error: erf for prec=32 (2)\n");
      mpfr_print_binary (x); printf ("\n");
      exit (1);
    }

  mpfr_set_str_binary (x, "100.10001110011110100000110000111");
  mpfr_set_str_binary (y, "0.11111111111111111111111111111111");
  mpfr_erf (x, x, GMP_RNDN);
  if (mpfr_cmp (x, y))
    {
      printf ("Error: erf for prec=32 (3)\n");
      exit (1);
    }
  mpfr_set_str_binary (x, "100.10001110011110100000110000111");
  mpfr_erf (x, x, GMP_RNDZ);
  if (mpfr_cmp (x, y))
    {
      printf ("Error: erf for prec=32 (4)\n");
      exit (1);
    }
  mpfr_set_str_binary (x, "100.10001110011110100000110000111");
  mpfr_erf (x, x, GMP_RNDU);
  if (mpfr_cmp_ui (x, 1))
    {
      printf ("Error: erf for prec=32 (5)\n");
      exit (1);
    }

  mpfr_set_str_binary (x, "100.10001110011110100000110001000");
  mpfr_erf (x, x, GMP_RNDN);
  if (mpfr_cmp_ui (x, 1))
    {
      printf ("Error: erf for prec=32 (6)\n");
      exit (1);
    }
  mpfr_set_str_binary (x, "100.10001110011110100000110001000");
  mpfr_set_str_binary (y, "0.11111111111111111111111111111111");
  mpfr_erf (x, x, GMP_RNDZ);
  if (mpfr_cmp (x, y))
    {
      printf ("Error: erf for prec=32 (7)\n");
      exit (1);
    }
  mpfr_set_str_binary (x, "100.10001110011110100000110001000");
  mpfr_erf (x, x, GMP_RNDU);
  if (mpfr_cmp_ui (x, 1))
    {
      printf ("Error: erf for prec=32 (8)\n");
      exit (1);
    }

  mpfr_set_ui (x, 5, GMP_RNDN);
  mpfr_erf (x, x, GMP_RNDN);
  if (mpfr_cmp_ui (x, 1))
    {
      printf ("Error: erf for prec=32 (9)\n");
      exit (1);
    }
  mpfr_set_ui (x, 5, GMP_RNDN);
  mpfr_erf (x, x, GMP_RNDU);
  if (mpfr_cmp_ui (x, 1))
    {
      printf ("Error: erf for prec=32 (10)\n");
      exit (1);
    }
  mpfr_set_ui (x, 5, GMP_RNDN);
  mpfr_erf (x, x, GMP_RNDZ);
  mpfr_set_str_binary (y, "0.11111111111111111111111111111111");
  if (mpfr_cmp (x, y))
    {
      printf ("Error: erf for prec=32 (11)\n");
      exit (1);
    }
  mpfr_set_ui (x, 5, GMP_RNDN);
  mpfr_erf (x, x, GMP_RNDD);
  mpfr_set_str_binary (y, "0.11111111111111111111111111111111");
  if (mpfr_cmp (x, y))
    {
      printf ("Error: erf for prec=32 (12)\n");
      exit (1);
    }

  mpfr_set_prec (x, 43);
  mpfr_set_prec (y, 64);
  mpfr_set_str_binary (x, "-0.1101110110101111100101011101110101101001001e3");
  mpfr_erf (y, x, GMP_RNDU);
  mpfr_set_prec (x, 64);
  mpfr_set_str_binary (x, "-0.1111111111111111111111111111111111111111111111111111111111111111");
  if (mpfr_cmp (x, y))
    {
      printf ("Error: erf for prec=43,64 (13)\n");
      exit (1);
    }
  

  mpfr_clear (x);
  mpfr_clear (y);
}

static void
special_erfc (void)
{
  mpfr_t x, y;

  mpfr_inits (x, y, NULL);

  /* erfc (NaN) = NaN */
  mpfr_set_nan (x);
  mpfr_erfc (y, x, GMP_RNDN);
  if (!mpfr_nan_p (y))
    {
      printf ("mpfr_erfc failed for x=NaN\n");
      exit (1);
    }
  /* erfc(+Inf) = 0+ */
  mpfr_set_inf (x, 1);
  mpfr_erfc (y, x, GMP_RNDN);
  if (!MPFR_IS_ZERO (y) || !MPFR_IS_POS (y))
    {
      printf ("mpfr_erf failed for x=+Inf\n");
      printf ("expected 0+, got ");
      mpfr_dump (y);
      exit (1);
    }
  /* erfc(-Inf) = 2 */
  mpfr_set_inf (x, -1);
  mpfr_erfc (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 2))
    {
      printf ("mpfr_erf failed for x=-Inf\n");
      printf ("expected 2, got ");
      mpfr_dump (y);
      exit (1);
    }
  /* erf(+0) = 1 */
  mpfr_set_ui (x, 0, GMP_RNDN);
  mpfr_erfc (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 1))
    {
      printf ("mpfr_erf failed for x=+0\n");
      printf ("expected 1, got ");
      mpfr_dump (y);
      exit (1);
    }

  mpfr_clears (x, y, NULL);
}

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  /* special_erf (); */
  /* special_erfc (); */

  /* test_generic_erf (2, 100, 15); */
  test_generic_erfc (2, 100, 15);

  tests_end_mpfr ();
  return 0;
}
