/* Test file for mpfr_ui_sub.

Copyright 2000, 2001, 2002, 2003, 2004 Free Software Foundation.

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
#include <float.h>
#include <time.h>

#include "mpfr-test.h"

static void
special (void)
{
  mpfr_t x, y, res;
  int inexact;

  mpfr_init (x);
  mpfr_init (y);
  mpfr_init (res);

  mpfr_set_prec (x, 24);
  mpfr_set_prec (y, 24);
  mpfr_set_str_binary (y, "0.111100110001011010111");
  inexact = mpfr_ui_sub (x, 1, y, GMP_RNDN);
  if (inexact)
    {
      printf ("Wrong inexact flag: got %d, expected 0\n", inexact);
      exit (1);
    }

  mpfr_set_prec (x, 24);
  mpfr_set_prec (y, 24);
  mpfr_set_str_binary (y, "0.111100110001011010111");
  if ((inexact = mpfr_ui_sub (x, 38181761, y, GMP_RNDN)) >= 0)
    {
      printf ("Wrong inexact flag: got %d, expected -1\n", inexact);
      exit (1);
    }

  mpfr_set_prec (x, 63);
  mpfr_set_prec (y, 63);
  mpfr_set_str_binary (y, "0.111110010010100100110101101010001001100101110001000101110111111E-1");
  if ((inexact = mpfr_ui_sub (x, 1541116494, y, GMP_RNDN)) <= 0)
    {
      printf ("Wrong inexact flag: got %d, expected +1\n", inexact);
      exit (1);
    }

  mpfr_set_prec (x, 32);
  mpfr_set_prec (y, 32);
  mpfr_set_str_binary (y, "0.11011000110111010001011100011100E-1");
  if ((inexact = mpfr_ui_sub (x, 2000375416, y, GMP_RNDN)) >= 0)
    {
      printf ("Wrong inexact flag: got %d, expected -1\n", inexact);
      exit (1);
    }

  mpfr_set_prec (x, 24);
  mpfr_set_prec (y, 24);
  mpfr_set_str_binary (y, "0.110011011001010011110111E-2");
  if ((inexact = mpfr_ui_sub (x, 927694848, y, GMP_RNDN)) <= 0)
    {
      printf ("Wrong inexact flag: got %d, expected +1\n", inexact);
      exit (1);
    }

  /* bug found by Mathieu Dutour, 12 Apr 2001 */
  mpfr_set_prec (x, 5);
  mpfr_set_prec (y, 5);
  mpfr_set_prec (res, 5);
  mpfr_set_str_binary (x, "1e-12");

  mpfr_ui_sub (y, 1, x, GMP_RNDD);
  mpfr_set_str_binary (res, "0.11111");
  if (mpfr_cmp (y, res))
    {
      printf ("Error in mpfr_ui_sub (y, 1, x, GMP_RNDD) for x=2^(-12)\nexpected 1.1111e-1, got ");
      mpfr_out_str (stdout, 2, 0, y, GMP_RNDN);
      printf ("\n");
      exit (1);
    }

  mpfr_ui_sub (y, 1, x, GMP_RNDU);
  mpfr_set_str_binary (res, "1.0");
  if (mpfr_cmp (y, res))
    {
      printf ("Error in mpfr_ui_sub (y, 1, x, GMP_RNDU) for x=2^(-12)\n"
              "expected 1.0, got ");
      mpfr_out_str (stdout, 2, 0, y, GMP_RNDN);
      printf ("\n");
      exit (1);
    }

  mpfr_ui_sub (y, 1, x, GMP_RNDN);
  mpfr_set_str_binary (res, "1.0");
  if (mpfr_cmp (y, res))
    {
      printf ("Error in mpfr_ui_sub (y, 1, x, GMP_RNDN) for x=2^(-12)\n"
              "expected 1.0, got ");
      mpfr_out_str (stdout, 2, 0, y, GMP_RNDN);
      printf ("\n");
      exit (1);
    }

  mpfr_set_prec (x, 10);
  mpfr_set_prec (y, 10);
  mpfr_random (x);
  mpfr_ui_sub (y, 0, x, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmpabs (x, y) == 0 && mpfr_sgn (x) != mpfr_sgn (y));

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (res);
}

/* checks that (y-x) gives the right results with 53 bits of precision */
static void
check (unsigned long y, const char *xs, mp_rnd_t rnd_mode, const char *zs)
{
  mpfr_t xx, zz;

  mpfr_inits2 (53, xx, zz, NULL);
  mpfr_set_str1 (xx, xs);
  mpfr_ui_sub (zz, y, xx, rnd_mode);
  if (mpfr_cmp_str1 (zz, zs) )
    {
      printf ("expected difference is %s, got\n",zs);
      mpfr_out_str(stdout, 10, 0, zz, GMP_RNDN);
      printf ("mpfr_ui_sub failed for y=%lu x=%s with rnd_mode=%s\n",
              y, xs, mpfr_print_rnd_mode (rnd_mode));
      exit (1);
    }
  mpfr_clears (xx, zz, NULL);
}

/* if u = o(x-y), v = o(u-x), w = o(v+y), then x-y = u-w */
static void
check_two_sum (mp_prec_t p)
{
  unsigned int x;
  mpfr_t y, u, v, w;
  mp_rnd_t rnd;
  int inexact;

  mpfr_inits2 (p, y, u, v, w, NULL);
  do
    {
      x = randlimb ();
    }
  while (x < 1);
  mpfr_random (y);
  rnd = GMP_RNDN;
  inexact = mpfr_ui_sub (u, x, y, rnd);
  mpfr_sub_ui (v, u, x, rnd);
  mpfr_add (w, v, y, rnd);
  /* as u = (x-y) + w, we should have inexact and w of same sign */
  if (((inexact == 0) && mpfr_cmp_ui (w, 0)) ||
      ((inexact > 0) && (mpfr_cmp_ui (w, 0) <= 0)) ||
      ((inexact < 0) && (mpfr_cmp_ui (w, 0) >= 0)))
    {
      printf ("Wrong inexact flag for prec=%u, rnd=%s\n",
              (unsigned int) p, mpfr_print_rnd_mode (rnd));
      printf ("x=%u\n", x);
      printf ("y="); mpfr_print_binary(y); puts ("");
      printf ("u="); mpfr_print_binary(u); puts ("");
      printf ("v="); mpfr_print_binary(v); puts ("");
      printf ("w="); mpfr_print_binary(w); puts ("");
      printf ("inexact = %d\n", inexact);
      exit (1);
    }
  mpfr_clears (y, u, v, w, NULL);
}

static void
check_nans (void)
{
  mpfr_t  x, y;

  mpfr_init2 (x, 123L);
  mpfr_init2 (y, 123L);

  /* 1 - nan == nan */
  mpfr_set_nan (x);
  mpfr_ui_sub (y, 1L, x, GMP_RNDN);
  MPFR_ASSERTN (mpfr_nan_p (y));

  /* 1 - +inf == -inf */
  mpfr_set_inf (x, 1);
  mpfr_ui_sub (y, 1L, x, GMP_RNDN);
  MPFR_ASSERTN (mpfr_inf_p (y));
  MPFR_ASSERTN (mpfr_sgn (y) < 0);

  /* 1 - -inf == +inf */
  mpfr_set_inf (x, -1);
  mpfr_ui_sub (y, 1L, x, GMP_RNDN);
  MPFR_ASSERTN (mpfr_inf_p (y));
  MPFR_ASSERTN (mpfr_sgn (y) > 0);

  mpfr_clear (x);
  mpfr_clear (y);
}

int
main (int argc, char *argv[])
{
  mp_prec_t p;
  unsigned k;

  tests_start_mpfr ();

  check_nans ();

  special ();
  for (p=2; p<100; p++)
    for (k=0; k<100; k++)
      check_two_sum (p);

  check(1196426492, "1.4218093058435347e-3", GMP_RNDN, 
	"1.1964264919985781e9");
  check(1092583421, "-1.0880649218158844e9", GMP_RNDN, 
	"2.1806483428158845901e9");
  check(948002822, "1.22191250737771397120e+20", GMP_RNDN,
	"-1.2219125073682338611e20");
  check(832100416, "4.68311314939691330000e-215", GMP_RNDD,
	"8.3210041599999988079e8");
  check(1976245324, "1.25296395864546893357e+232", GMP_RNDZ,
	"-1.2529639586454686577e232");
  check(2128997392, "-1.08496826129284207724e+187", GMP_RNDU,
	"1.0849682612928422704e187");
  check(293607738, "-1.9967571564050541e-5", GMP_RNDU, 
	"2.9360773800002003e8");
  check(354270183, "2.9469161763489528e3", GMP_RNDN, 
	"3.5426723608382362e8");

  tests_end_mpfr ();
  return 0;
}
