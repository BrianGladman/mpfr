/* Test file for mpfr_cos.

Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 Free Software Foundation, Inc.
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

#ifdef CHECK_EXTERNAL
static int
test_cos (mpfr_ptr a, mpfr_srcptr b, mp_rnd_t rnd_mode)
{
  int res;
  int ok = rnd_mode == GMP_RNDN && mpfr_number_p (b) && mpfr_get_prec (a)>=53;
  if (ok)
    {
      mpfr_print_raw (b);
    }
  res = mpfr_cos (a, b, rnd_mode);
  if (ok)
    {
      printf (" ");
      mpfr_print_raw (a);
      printf ("\n");
    }
  return res;
}
#else
#define test_cos mpfr_cos
#endif

static void
check53 (const char *xs, const char *cos_xs, mp_rnd_t rnd_mode)
{
  mpfr_t xx, c;

  mpfr_inits2 (53, xx, c, (void *) 0);
  mpfr_set_str1 (xx, xs); /* should be exact */
  test_cos (c, xx, rnd_mode);
  if (mpfr_cmp_str1 (c, cos_xs))
    {
      printf ("mpfr_cos failed for x=%s, rnd=%s\n",
              xs, mpfr_print_rnd_mode (rnd_mode));
      printf ("mpfr_cos gives cos(x)=");
      mpfr_out_str(stdout, 10, 0, c, GMP_RNDN);
      printf(", expected %s\n", cos_xs);
      exit (1);
    }
  mpfr_clears (xx, c, (void *) 0);
}

#define TEST_FUNCTION test_cos
#define REDUCE_EMAX 1048575 /* otherwise arg. reduction is too expensive */
#include "tgeneric.c"

static void
check_nans (void)
{
  mpfr_t  x, y;

  mpfr_init2 (x, 123L);
  mpfr_init2 (y, 123L);

  mpfr_set_nan (x);
  test_cos (y, x, GMP_RNDN);
  if (! mpfr_nan_p (y))
    {
      printf ("Error: cos(NaN) != NaN\n");
      exit (1);
    }

  mpfr_set_inf (x, 1);
  test_cos (y, x, GMP_RNDN);
  if (! mpfr_nan_p (y))
    {
      printf ("Error: cos(Inf) != NaN\n");
      exit (1);
    }

  mpfr_set_inf (x, -1);
  test_cos (y, x, GMP_RNDN);
  if (! mpfr_nan_p (y))
    {
      printf ("Error: cos(-Inf) != NaN\n");
      exit (1);
    }

  /* cos(+/-0) = 1 */
  mpfr_set_ui (x, 0, GMP_RNDN);
  test_cos (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 1))
    {
      printf ("Error: cos(+0) != 1\n");
      exit (1);
    }
  mpfr_neg (x, x, GMP_RNDN);
  test_cos (y, x, GMP_RNDN);
  if (mpfr_cmp_ui (y, 1))
    {
      printf ("Error: cos(-0) != 1\n");
      exit (1);
    }

  /* Compute ~Pi/2 to check */
  /* FIXME: Too slow!
  mpfr_set_prec (x, 20000);
  mpfr_const_pi (x, GMP_RNDD); mpfr_div_2ui (x, x, 1, GMP_RNDN);
  mpfr_set_prec (y, 24);
  test_cos (y, x, GMP_RNDN);
  if (mpfr_cmp_str (y, "0.111001010110100011000001E-20000", 2, GMP_RNDN))
    {
      printf("Error computing cos(~Pi/2)\n");
      mpfr_dump (y);
      exit (1);
      } */

  mpfr_clear (x);
  mpfr_clear (y);
}

static void
special_overflow (void)
{
  mpfr_t x, y;
  mp_exp_t emin, emax;

  emin = mpfr_get_emin ();
  emax = mpfr_get_emax ();

  mpfr_init2 (x, 24);
  mpfr_init2 (y, 73);

  /* Check special case: An overflow in const_pi could occurs! */
  set_emin (-125);
  set_emax (128);
  mpfr_set_str_binary (x, "0.111101010110110011101101E6");
  test_cos (y, x, GMP_RNDZ);
  set_emin (emin);
  set_emax (emax);

  mpfr_clear (x);
  mpfr_clear (y);
}

static void
overflowed_cos0 (void)
{
  mpfr_t x, y;
  int emax, i, inex, rnd, err = 0;
  mp_exp_t old_emax;

  old_emax = mpfr_get_emax ();

  mpfr_init2 (x, 8);
  mpfr_init2 (y, 8);

  for (emax = -1; emax <= 0; emax++)
    {
      mpfr_set_ui_2exp (y, 1, emax, GMP_RNDN);
      mpfr_nextbelow (y);
      set_emax (emax);  /* 1 is not representable. */
      /* and if emax < 0, 1 - eps is not representable either. */
      for (i = -1; i <= 1; i++)
        RND_LOOP (rnd)
        {
          mpfr_set_si_2exp (x, i, -512 * ABS (i), GMP_RNDN);
          mpfr_clear_flags ();
          inex = mpfr_cos (x, x, rnd);
          if ((i == 0 || emax < 0 || rnd == GMP_RNDN || rnd == GMP_RNDU) &&
              ! mpfr_overflow_p ())
            {
              printf ("Error in overflowed_cos0 (i = %d, rnd = %s):\n"
                      "  The overflow flag is not set.\n",
                      i, mpfr_print_rnd_mode (rnd));
              err = 1;
            }
          if (rnd == GMP_RNDZ || rnd == GMP_RNDD)
            {
              if (inex >= 0)
                {
                  printf ("Error in overflowed_cos0 (i = %d, rnd = %s):\n"
                          "  The inexact value must be negative.\n",
                          i, mpfr_print_rnd_mode (rnd));
                  err = 1;
                }
              if (! mpfr_equal_p (x, y))
                {
                  printf ("Error in overflowed_cos0 (i = %d, rnd = %s):\n"
                          "  Got ", i, mpfr_print_rnd_mode (rnd));
                  mpfr_print_binary (x);
                  printf (" instead of 0.11111111E%d.\n", emax);
                  err = 1;
                }
            }
          else
            {
              if (inex <= 0)
                {
                  printf ("Error in overflowed_cos0 (i = %d, rnd = %s):\n"
                          "  The inexact value must be positive.\n",
                          i, mpfr_print_rnd_mode (rnd));
                  err = 1;
                }
              if (! (mpfr_inf_p (x) && MPFR_SIGN (x) > 0))
                {
                  printf ("Error in overflowed_cos0 (i = %d, rnd = %s):\n"
                          "  Got ", i, mpfr_print_rnd_mode (rnd));
                  mpfr_print_binary (x);
                  printf (" instead of +Inf.\n");
                  err = 1;
                }
            }
        }
      set_emax (old_emax);
    }

  if (err)
    exit (1);
  mpfr_clear (x);
  mpfr_clear (y);
}

int
main (int argc, char *argv[])
{
  mpfr_t x, y;
  int inex;

  tests_start_mpfr ();

  special_overflow ();
  check_nans ();

  mpfr_init (x);
  mpfr_init (y);

  mpfr_set_prec (x, 53);
  mpfr_set_prec (y, 2);
  mpfr_set_str (x, "9.81333845856942e-1", 10, GMP_RNDN);
  test_cos (y, x, GMP_RNDN);

  mpfr_set_prec (x, 30);
  mpfr_set_prec (y, 30);
  mpfr_set_str_binary (x, "1.00001010001101110010100010101e-1");
  test_cos (y, x, GMP_RNDU);
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
  test_cos (y, x, GMP_RNDU);
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
  test_cos (y, x, GMP_RNDD);
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
  test_cos (x, x, GMP_RNDN);
  if (mpfr_cmp (x, y))
    {
      printf ("Error for prec=32 (1)\n");
      exit (1);
    }

  mpfr_set_str_binary (x, "-0.1101011110111100111010011001011E-1");
  mpfr_set_str_binary (y, "0.11101001100110111011011010100011");
  test_cos (x, x, GMP_RNDN);
  if (mpfr_cmp (x, y))
    {
      printf ("Error for prec=32 (2)\n");
      exit (1);
    }

  /* huge argument reduction */
  mpfr_set_str_binary (x, "0.10000010000001101011101111001011E40");
  mpfr_set_str_binary (y, "0.10011000001111010000101011001011E-1");
  test_cos (x, x, GMP_RNDN);
  if (mpfr_cmp (x, y))
    {
      printf ("Error for prec=32 (3)\n");
      exit (1);
    }

  mpfr_set_prec (x, 3);
  mpfr_set_prec (y, 3);
  mpfr_set_str_binary (x, "0.110E60");
  inex = mpfr_cos (y, x, GMP_RNDD);
  MPFR_ASSERTN(inex < 0);

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

  overflowed_cos0 ();
  test_generic (2, 100, 15);

  mpfr_clear (x);
  mpfr_clear (y);

  tests_end_mpfr ();
  return 0;
}
