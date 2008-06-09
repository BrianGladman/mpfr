/* Test file for mpfr_hypot.

Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 Free Software Foundation, Inc.
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
#include <limits.h>
#include <stdlib.h>

#include "mpfr-test.h"

static void
special (void)
{
  mpfr_t x, y, z;

  mpfr_init (x);
  mpfr_init (y);
  mpfr_init (z);

  mpfr_set_nan (x);
  mpfr_hypot (z, x, y, GMP_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (z));

  mpfr_set_inf (x, 1);
  mpfr_set_inf (y, -1);
  mpfr_hypot (z, x, y, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (z) && mpfr_sgn (z) > 0);

  mpfr_set_inf (x, -1);
  mpfr_set_nan (y);
  mpfr_hypot (z, x, y, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (z) && mpfr_sgn (z) > 0);

  mpfr_set_nan (x);
  mpfr_set_inf (y, -1);
  mpfr_hypot (z, x, y, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (z) && mpfr_sgn (z) > 0);

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
}

static void
test_large (void)
{
  mpfr_t x, y, z, t;

  mpfr_init (x);
  mpfr_init (y);
  mpfr_init (z);
  mpfr_init (t);

  mpfr_set_ui (x, 21, GMP_RNDN);
  mpfr_set_ui (y, 28, GMP_RNDN);
  mpfr_set_ui (z, 35, GMP_RNDN);

  mpfr_mul_2ui (x, x, MPFR_EMAX_DEFAULT-6, GMP_RNDN);
  mpfr_mul_2ui (y, y, MPFR_EMAX_DEFAULT-6, GMP_RNDN);
  mpfr_mul_2ui (z, z, MPFR_EMAX_DEFAULT-6, GMP_RNDN);

  mpfr_hypot (t, x, y, GMP_RNDN);
  if (mpfr_cmp (z, t))
    {
      printf ("Error in test_large: got\n");
      mpfr_out_str (stdout, 2, 0, t, GMP_RNDN);
      printf ("\ninstead of\n");
      mpfr_out_str (stdout, 2, 0, z, GMP_RNDN);
      printf ("\n");
      exit (1);
    }

  mpfr_set_prec (x, 53);
  mpfr_set_prec (t, 53);
  mpfr_set_prec (y, 53);
  mpfr_set_str_binary (x, "0.11101100011110000011101000010101010011001101000001100E-1021");
  mpfr_set_str_binary (t, "0.11111001010011000001110110001101011100001000010010100E-1021");
  mpfr_hypot (y, x, t, GMP_RNDN);

  mpfr_set_prec (x, 240);
  mpfr_set_prec (y, 22);
  mpfr_set_prec (t, 2);
  mpfr_set_str_binary (x, "0.100111011010010010110100000100000001100010011100110101101111111101011110111011011101010110100101111000111100010100110000100101011110111011100110100110100101110101101100011000001100000001111101110100100100011011011010110111100110010101000111e-7");
  mpfr_set_str_binary (y, "0.1111000010000011000111e-10");
  mpfr_hypot (t, x, y, GMP_RNDN);
  mpfr_set_str_binary (y, "0.11E-7");
  if (mpfr_cmp (t, y))
    {
      printf ("Error in mpfr_hypot (1)\n");
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
  mpfr_clear (t);
}

static void
test_large_small (void)
{
  mpfr_t x, y, z;
  int inexact;

  mpfr_init2 (x, 3);
  mpfr_init2 (y, 2);
  mpfr_init2 (z, 2);

  mpfr_set_ui_2exp (x, 1, mpfr_get_emax () / 2, GMP_RNDN);
  mpfr_set_ui_2exp (y, 1, -1, GMP_RNDN);
  inexact = mpfr_hypot (z, x, y, GMP_RNDN);
  if (inexact >= 0 || mpfr_cmp (x, z))
    {
      printf ("Error 1 in test_large_small\n");
      exit (1);
    }

  mpfr_mul_ui (x, x, 5, GMP_RNDN);
  inexact = mpfr_hypot (z, x, y, GMP_RNDN);
  if (mpfr_cmp (x, z) >= 0)
    {
      printf ("Error 2 in test_large_small\n");
      printf ("x = ");
      mpfr_out_str (stdout, 2, 0, x, GMP_RNDN);
      printf ("\n");
      printf ("y = ");
      mpfr_out_str (stdout, 2, 0, y, GMP_RNDN);
      printf ("\n");
      printf ("z = ");
      mpfr_out_str (stdout, 2, 0, z, GMP_RNDN);
      printf (" (in precision 2) instead of\n    ");
      mpfr_out_str (stdout, 2, 2, x, GMP_RNDU);
      printf ("\n");
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
}

static void
check_overflow (void)
{
  mpfr_t x, y;
  int inex, r;

  mpfr_inits2 (8, x, y, (mpfr_ptr) 0);
  mpfr_set_ui (x, 1, GMP_RNDN);
  mpfr_setmax (x, mpfr_get_emax ());

  RND_LOOP(r)
    {
      mpfr_clear_overflow ();
      inex = mpfr_hypot (y, x, x, (mp_rnd_t) r);
      if (!mpfr_overflow_p ())
        {
          printf ("No overflow in check_overflow for %s\n",
                  mpfr_print_rnd_mode ((mp_rnd_t) r));
          exit (1);
        }
      MPFR_ASSERTN (MPFR_IS_POS (y));
      if (r == GMP_RNDZ || r == GMP_RNDD)
        {
          MPFR_ASSERTN (inex < 0);
          MPFR_ASSERTN (!mpfr_inf_p (y));
          mpfr_nexttoinf (y);
        }
      else
        {
          MPFR_ASSERTN (inex > 0);
        }
      MPFR_ASSERTN (mpfr_inf_p (y));
    }

  mpfr_clears (x, y, (mpfr_ptr) 0);
}

#define TWO_ARGS
#define TEST_FUNCTION mpfr_hypot
#include "tgeneric.c"

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  special ();

  test_large ();
  test_large_small ();
  check_overflow ();

  test_generic (2, 100, 10);

  tests_end_mpfr ();
  return 0;
}
