/* Test file for exceptions.

Copyright 2001, 2002, 2003, 2004 Free Software Foundation.

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

/* Test default rounding mode */
static void
check_default_rnd(void)
{
  mp_rnd_t r, t;
  for(r = 0 ; r < 4 ; r++)
    {
      mpfr_set_default_rounding_mode (r);
      t = mpfr_get_default_rounding_mode();
      if (r !=t)
	{
	  printf("ERROR in setting / getting default rounding mode (1)\n");
	  exit(1);
	}
    }
  mpfr_set_default_rounding_mode(4);
  if (mpfr_get_default_rounding_mode() != GMP_RNDD)
    {
      printf("ERROR in setting / getting default rounding mode (2)\n");
      exit(1);
    }
  mpfr_set_default_rounding_mode(-1);
  if (mpfr_get_default_rounding_mode() != GMP_RNDD)
    {
      printf("ERROR in setting / getting default rounding mode (3)\n");
      exit(1);
    }
}

static void
mpfr_set_double_range (void)
{
  mpfr_set_default_prec (53);

  /* in double precision format, the unbiased exponent is between 0 and
     2047, where 0 is used for subnormal numbers, and 2047 for special
     numbers (infinities, NaN), and the bias is 1023, thus "normal" numbers
     have an exponent between -1022 and 1023, corresponding to numbers
     between 2^(-1022) and previous(2^(1024)).
     (The smallest subnormal number is 0.(0^51)1*2^(-1022)= 2^(-1074).)

     The smallest normal power of two is 1.0*2^(-1022).
     The largest normal power of two is 2^1023.
     (We have to add one for mpfr since mantissa are between 1/2 and 1.)
  */

  mpfr_set_emin (-1021);
  mpfr_set_emax (1024);
}

static void
test_set_underflow (void)
{
  /* static to allow non-constant initialiers in r */
  mpfr_t x, zero, min;
  mpfr_ptr r[4];
  int t[4] = { 1, -1, 1, -1 };
  mp_rnd_t i;
  int s;

  mpfr_inits (x, zero, min, (mpfr_ptr) 0);
  mpfr_set_ui (zero, 0, GMP_RNDN);
  mpfr_set_ui (min, 0, GMP_RNDN);
  mpfr_nextabove (min);
  r[0] = r[2] = min;
  r[1] = r[3] = zero;
  for (s = 1; s > 0; s = -1)
    {
      for (i = 0; i < 4; i++)
        {
          int j;
          int inex;

          j = s < 0 && i > 1 ? 5 - i : i;
          inex = mpfr_set_underflow (x, i, s);
          if (mpfr_cmp (x, r[j]) || inex * t[j] <= 0)
            {
              printf ("Error in test_set_underflow, sign = %d,"
                      " rnd_mode = %s\n", s, mpfr_print_rnd_mode (i));
              printf ("Got\n");
              mpfr_out_str (stdout, 2, 0, x, GMP_RNDN);
              printf (", inex = %d\ninstead of\n", inex);
              mpfr_out_str (stdout, 2, 0, r[j], GMP_RNDN);
              printf (", inex = %d\n", t[j]);
              exit (1);
            }
        }
      mpfr_neg (zero, zero, GMP_RNDN);
      mpfr_neg (min, min, GMP_RNDN);
    }
  mpfr_clears (x, zero, min, (mpfr_ptr) 0);
}

static void
test_set_overflow (void)
{
  /* static to allow non-constant initialiers in r */
  mpfr_t x, inf, max;
  mpfr_ptr r[4];
  int t[4] = { 1, -1, 1, -1 };
  mp_rnd_t i;
  int s;

  mpfr_inits2 (32, x, inf, max, (mpfr_ptr) 0);
  mpfr_set_inf (inf, 1);
  mpfr_set_inf (max, 1);
  mpfr_nextbelow (max);
  r[0] = r[2] = inf;
  r[1] = r[3] = max;
  for (s = 1; s > 0; s = -1)
    {
      for (i = 0; i < 4; i++)
        {
          int j;
          int inex;

          j = s < 0 && i > 1 ? 5 - i : i;
          inex = mpfr_set_overflow (x, i, s);
          if (mpfr_cmp (x, r[j]) || inex * t[j] <= 0)
            {
              printf ("Error in test_set_overflow, sign = %d,"
                      " rnd_mode = %s\n", s, mpfr_print_rnd_mode (i));
              printf ("Got\n");
              mpfr_out_str (stdout, 2, 0, x, GMP_RNDN);
              printf (", inex = %d\ninstead of\n", inex);
              mpfr_out_str (stdout, 2, 0, r[j], GMP_RNDN);
              printf (", inex = %d\n", t[j]);
              exit (1);
            }
        }
      mpfr_neg (inf, inf, GMP_RNDN);
      mpfr_neg (max, max, GMP_RNDN);
    }
  mpfr_clears (x, inf, max, (mpfr_ptr) 0);
}

int
main (int argc, char *argv[])
{
  mpfr_t x, y;
  mp_exp_t emin, emax;

  tests_start_mpfr ();

  test_set_underflow ();
  test_set_overflow ();
  check_default_rnd();

  mpfr_init (x);
  mpfr_init (y);

  emin = mpfr_get_emin ();
  emax = mpfr_get_emax ();
  if (emin >= emax)
    {
      printf ("Error: emin >= emax\n");
      exit (1);
    }

  mpfr_set_ui (x, 1, GMP_RNDN);
  mpfr_mul_2exp (x, x, 1024, GMP_RNDN);
  mpfr_set_double_range ();
  mpfr_check_range (x, 0, GMP_RNDN);
  if (!mpfr_inf_p (x) || (mpfr_sgn(x) <= 0))
    {
      printf ("Error: 2^1024 rounded to nearest should give +Inf\n");
      exit (1);
    }

  mpfr_set_emax (1025);
  mpfr_set_ui (x, 1, GMP_RNDN);
  mpfr_mul_2exp (x, x, 1024, GMP_RNDN);
  mpfr_set_double_range ();
  mpfr_check_range (x, 0, GMP_RNDD);
  if (!mpfr_number_p (x))
    {
      printf ("Error: 2^1024 rounded down should give a normal number\n");
      exit (1);
    }

  mpfr_set_ui (x, 1, GMP_RNDN);
  mpfr_mul_2exp (x, x, 1023, GMP_RNDN);
  mpfr_add (x, x, x, GMP_RNDN);
  if (!mpfr_inf_p (x) || (mpfr_sgn(x) <= 0))
    {
      printf ("Error: x+x rounded to nearest for x=2^1023 should give +Inf\n");
      printf ("emax = %ld\n", mpfr_get_emax ());
      printf ("got "); mpfr_print_binary (x); puts ("");
      exit (1);
    }

  mpfr_set_ui (x, 1, GMP_RNDN);
  mpfr_mul_2exp (x, x, 1023, GMP_RNDN);
  mpfr_add (x, x, x, GMP_RNDD);
  if (!mpfr_number_p (x))
    {
      printf ("Error: x+x rounded down for x=2^1023 should give"
              " a normal number\n");
      exit (1);
    }

  mpfr_set_ui (x, 1, GMP_RNDN);
  mpfr_div_2exp (x, x, 1022, GMP_RNDN);
  mpfr_set_str_binary (y, "1.1e-1022"); /* y = 3/2*x */
  mpfr_sub (y, y, x, GMP_RNDZ);
  if (mpfr_cmp_ui (y, 0))
    {
      printf ("Error: y-x rounded to zero should give 0"
              " for y=3/2*2^(-1022), x=2^(-1022)\n");
      printf ("y="); mpfr_print_binary (y); puts ("");
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);

  tests_end_mpfr ();
  return 0;
}
