/* Test file for mpfr_cos.

Copyright 2001, 2002, 2003, 2004 Free Software Foundation, Inc.

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

static void
check53 (const char *xs, const char *cos_xs, mp_rnd_t rnd_mode)
{
  mpfr_t xx, c;

  mpfr_inits2 (53, xx, c, NULL);
  mpfr_set_str1 (xx, xs); /* should be exact */
  mpfr_cos (c, xx, rnd_mode);
  if (mpfr_cmp_str1 (c, cos_xs))
    {
      printf ("mpfr_cos failed for x=%s, rnd=%s\n", 
	      xs, mpfr_print_rnd_mode (rnd_mode));
      printf ("mpfr_cos gives cos(x)=");
      mpfr_out_str(stdout, 10, 0, c, GMP_RNDN);
      printf(", expected %s\n", cos_xs);
      exit (1);
    }
  mpfr_clears (xx, c, NULL);
}

#define TEST_FUNCTION mpfr_cos
#include "tgeneric.c"

static void
check_nans (void)
{
  mpfr_t  x, y;

  mpfr_init2 (x, 123L);
  mpfr_init2 (y, 123L);

  mpfr_set_nan (x);
  mpfr_cos (y, x, GMP_RNDN);
  if (! mpfr_nan_p (y))
    {
      printf ("Error: cos(NaN) != NaN\n");
      exit (1);
    }

  mpfr_set_inf (x, 1);
  mpfr_cos (y, x, GMP_RNDN);
  if (! mpfr_nan_p (y))
    {
      printf ("Error: cos(Inf) != NaN\n");
      exit (1);
    }

  mpfr_set_inf (x, -1);
  mpfr_cos (y, x, GMP_RNDN);
  if (! mpfr_nan_p (y))
    {
      printf ("Error: cos(-Inf) != NaN\n");
      exit (1);
    }

  /* cos(+/-0) = 1 */
  mpfr_set_ui (x, 0, GMP_RNDN);
  mpfr_cos (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 1))
    {
      printf ("Error: cos(+0) != 1\n");
      exit (1);
    }
  mpfr_neg (x, x, GMP_RNDN);
  mpfr_cos (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 1))
    {
      printf ("Error: cos(-0) != 1\n");
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);
}

int
main (int argc, char *argv[])
{
  mpfr_t x, y;

  tests_start_mpfr ();

  check_nans ();

  mpfr_init (x);
  mpfr_init (y);

  mpfr_set_prec (x, 53);
  mpfr_set_prec (y, 2);
  mpfr_set_str (x, "9.81333845856942e-1", 10, GMP_RNDN);
  mpfr_cos (y, x, GMP_RNDN);

  mpfr_set_prec (x, 30);
  mpfr_set_prec (y, 30);
  mpfr_set_str_binary (x, "1.00001010001101110010100010101e-1");
  mpfr_cos (y, x, GMP_RNDU);
  mpfr_set_str_binary (x, "1.10111100010101011110101010100e-1");
  if (mpfr_cmp (y, x))
    {
      printf ("Error for prec=30, rnd=GMP_RNDU\n");
      printf ("expected "); mpfr_print_binary (x); puts ("");
      printf ("     got "); mpfr_print_binary (y); puts ("");
      exit (1);
    }

  mpfr_set_prec (x, 59);
  mpfr_set_prec (y, 59);
  mpfr_set_str_binary (x, "1.01101011101111010011111110111111111011011101100111100011e-3");
  mpfr_cos (y, x, GMP_RNDU);
  mpfr_set_str_binary (x, "1.1111011111110010001001001011100111101110100010000010010011e-1");
  if (mpfr_cmp (y, x))
    {
      printf ("Error for prec=59, rnd=GMP_RNDU\n");
      printf ("expected "); mpfr_print_binary (x); puts ("");
      printf ("     got "); mpfr_print_binary (y); puts ("");
      exit (1);
    }

  mpfr_set_prec (x, 5);
  mpfr_set_prec (y, 5);
  mpfr_set_str_binary (x, "1.1100e-2");
  mpfr_cos (y, x, GMP_RNDD);
  mpfr_set_str_binary (x, "1.1100e-1");
  if (mpfr_cmp (y, x))
    {
      printf ("Error for x=1.1100e-2, rnd=GMP_RNDD\n");
      printf ("expected 1.1100e-1, got "); mpfr_print_binary (y); puts ("");
      exit (1);
    }

  mpfr_set_prec (x, 32);
  mpfr_set_prec (y, 32);

  mpfr_set_str_binary (x, "0.10001000001001011000100001E-6");
  mpfr_set_str_binary (y, "0.1111111111111101101111001100001");
  mpfr_cos (x, x, GMP_RNDN);
  if (mpfr_cmp (x, y))
    {
      printf ("Error for prec=32 (1)\n");
      exit (1);
    }

  mpfr_set_str_binary (x, "-0.1101011110111100111010011001011E-1");
  mpfr_set_str_binary (y, "0.11101001100110111011011010100011");
  mpfr_cos (x, x, GMP_RNDN);
  if (mpfr_cmp (x, y))
    {
      printf ("Error for prec=32 (2)\n");
      exit (1);
    }

  /* huge argument reduction */
  mpfr_set_str_binary (x, "0.10000010000001101011101111001011E40");
  mpfr_set_str_binary (y, "0.10011000001111010000101011001011E-1");
  mpfr_cos (x, x, GMP_RNDN);
  if (mpfr_cmp (x, y))
    {
      printf ("Error for prec=32 (3)\n");
      exit (1);
    }

  /* worst case from PhD thesis of Vincent Lefe`vre: x=8980155785351021/2^54 */
  check53 ("4.984987858808754279e-1", "8.783012931285841817e-1", GMP_RNDN);
  check53 ("4.984987858808754279e-1", "8.783012931285840707e-1", GMP_RNDD);
  check53 ("4.984987858808754279e-1", "8.783012931285840707e-1", GMP_RNDZ);
  check53 ("4.984987858808754279e-1", "8.783012931285841817e-1", GMP_RNDU);
  check53 ("1.00031274099908640274",  "0.540039116973283217504", GMP_RNDN);
  check53 ("1.00229256850978698523",  "0.538371757797526551137", GMP_RNDZ);
  check53 ("1.00288304857059840103",  "0.537874062022526966409", GMP_RNDZ);
  check53 ("1.00591265847407274059",  "0.53531755997839769456",  GMP_RNDN);

  check53 ("1.00591265847407274059", "0.53531755997839769456",  GMP_RNDN);

  test_generic (2, 100, 100);

  mpfr_clear (x);
  mpfr_clear (y);

  tests_end_mpfr ();
  return 0;
}
