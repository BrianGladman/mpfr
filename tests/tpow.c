/* Test file for mpfr_pow, mpfr_pow_ui and mpfr_pow_si.

Copyright 2000, 2001, 2002, 2003, 2004, 2005 Free Software Foundation, Inc.

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
#include <math.h>

#include "mpfr-test.h"

#ifdef CHECK_EXTERNAL
static int
test_pow (mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mp_rnd_t rnd_mode)
{
  int res;
  int ok = rnd_mode == GMP_RNDN && mpfr_number_p (b) && mpfr_number_p (c)
    && mpfr_get_prec (a) >= 53;
  if (ok)
    {
      mpfr_print_raw (b);
      printf (" ");
      mpfr_print_raw (c);
    }
  res = mpfr_pow (a, b, c, rnd_mode);
  if (ok)
    {
      printf (" ");
      mpfr_print_raw (a);
      printf ("\n");
    }
  return res;
}
#else
#define test_pow mpfr_pow
#endif

#define TEST_FUNCTION test_pow
#define TWO_ARGS
#include "tgeneric.c"

static void
check_pow_ui (void)
{
  mpfr_t a, b;
  int res;

  mpfr_init2 (a, 53);
  mpfr_init2 (b, 53);

  /* check in-place operations */
  mpfr_set_str (b, "0.6926773", 10, GMP_RNDN);
  mpfr_pow_ui (a, b, 10, GMP_RNDN);
  mpfr_pow_ui (b, b, 10, GMP_RNDN);
  if (mpfr_cmp (a, b))
    {
      printf ("Error for mpfr_pow_ui (b, b, ...)\n");
      exit (1);
    }

  /* check large exponents */
  mpfr_set_ui (b, 1, GMP_RNDN);
  mpfr_pow_ui (a, b, 4294967295UL, GMP_RNDN);

  mpfr_set_inf (a, -1);
  mpfr_pow_ui (a, a, 4049053855UL, GMP_RNDN);
  if (!mpfr_inf_p (a) || (mpfr_sgn (a) >= 0))
    {
      printf ("Error for (-Inf)^4049053855\n");
      exit (1);
    }

  mpfr_set_inf (a, -1);
  mpfr_pow_ui (a, a, (unsigned long) 30002752, GMP_RNDN);
  if (!mpfr_inf_p (a) || (mpfr_sgn (a) <= 0))
    {
      printf ("Error for (-Inf)^30002752\n");
      exit (1);
    }

  /* Check underflow */
  mpfr_set_str_binary (a, "1E-1");
  res = mpfr_pow_ui (a, a, -mpfr_get_emin (), GMP_RNDN);
  if (MPFR_GET_EXP (a) != mpfr_get_emin () + 1)
    {
      printf ("Error for (1e-1)^MPFR_EMAX_MAX\n");
      mpfr_dump (a);
      exit (1);
    }

  mpfr_set_str_binary (a, "1E-10");
  res = mpfr_pow_ui (a, a, -mpfr_get_emin (), GMP_RNDZ);
  if (!MPFR_IS_ZERO (a))
    {
      printf ("Error for (1e-10)^MPFR_EMAX_MAX\n");
      mpfr_dump (a);
      exit (1);
    }

  /* Check overflow */
  mpfr_set_str_binary (a, "1E10");
  res = mpfr_pow_ui (a, a, ULONG_MAX, GMP_RNDN);
  if (!MPFR_IS_INF (a) || MPFR_SIGN (a) < 0)
    {
      printf ("Error for (1e10)^ULONG_MAX\n");
      exit (1);
    }

  /* Check 0 */
  MPFR_SET_ZERO (a);
  MPFR_SET_POS (a);
  mpfr_set_si (b, -1, GMP_RNDN);
  res = mpfr_pow_ui (b, a, 1, GMP_RNDN);
  if (res != 0 || MPFR_IS_NEG (b))
    {
      printf ("Error for (0+)^1\n");
      exit (1);
    }
  MPFR_SET_ZERO (a);
  MPFR_SET_NEG (a);
  mpfr_set_ui (b, 1, GMP_RNDN);
  res = mpfr_pow_ui (b, a, 5, GMP_RNDN);
  if (res != 0 || MPFR_IS_POS (b))
    {
      printf ("Error for (0-)^5\n");
      exit (1);
    }
  MPFR_SET_ZERO (a);
  MPFR_SET_NEG (a);
  mpfr_set_si (b, -1, GMP_RNDN);
  res = mpfr_pow_ui (b, a, 6, GMP_RNDN);
  if (res != 0 || MPFR_IS_NEG (b))
    {
      printf ("Error for (0-)^6\n");
      exit (1);
    }

  mpfr_set_prec (a, 122);
  mpfr_set_prec (b, 122);
  mpfr_set_str_binary (a, "0.10000010010000111101001110100101101010011110011100001111000001001101000110011001001001001011001011010110110110101000111011E1");
  mpfr_set_str_binary (b, "0.11111111100101001001000001000001100011100000001110111111100011111000111011100111111111110100011000111011000100100011001011E51290375");
  mpfr_pow_ui (a, a, 2026876995UL, GMP_RNDU);
  if (mpfr_cmp (a, b) != 0)
    {
      printf ("Error for x^2026876995\n");
      exit (1);
    }

  mpfr_clear (a);
  mpfr_clear (b);
}

static void
check_pow_si (void)
{
  mpfr_t x;

  mpfr_init (x);

  mpfr_set_nan (x);
  mpfr_pow_si (x, x, -1, GMP_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (x));

  mpfr_set_inf (x, 1);
  mpfr_pow_si (x, x, -1, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (x, 0) == 0 && MPFR_IS_POS(x));

  mpfr_set_inf (x, -1);
  mpfr_pow_si (x, x, -1, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (x, 0) == 0 && MPFR_IS_NEG(x));

  mpfr_set_inf (x, -1);
  mpfr_pow_si (x, x, -2, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (x, 0) == 0 && MPFR_IS_POS(x));

  mpfr_set_ui (x, 0, GMP_RNDN);
  mpfr_pow_si (x, x, -1, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (x) && mpfr_sgn (x) > 0);

  mpfr_set_ui (x, 0, GMP_RNDN);
  mpfr_neg (x, x, GMP_RNDN);
  mpfr_pow_si (x, x, -1, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (x) && mpfr_sgn (x) < 0);

  mpfr_set_ui (x, 0, GMP_RNDN);
  mpfr_neg (x, x, GMP_RNDN);
  mpfr_pow_si (x, x, -2, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (x) && mpfr_sgn (x) > 0);

  mpfr_clear (x);
}

static void
check_special_pow_si ()
{
  mpfr_t a, b;
  mp_exp_t emin;

  mpfr_init (a);
  mpfr_init (b);
  mpfr_set_str (a, "2E100000000", 10, GMP_RNDN);
  mpfr_set_si (b, -10, GMP_RNDN);
  test_pow (b, a, b, GMP_RNDN);
  if (!MPFR_IS_ZERO(b))
    {
      printf("Pow(2E10000000, -10) failed\n");
      mpfr_dump (a);
      mpfr_dump (b);
      exit(1);
    }

  emin = mpfr_get_emin ();
  mpfr_set_emin (-10);
  mpfr_set_si (a, -2, GMP_RNDN);
  mpfr_pow_si (b, a, -10000, GMP_RNDN);
  if (!MPFR_IS_ZERO (b))
    {
      printf ("Pow_so (1, -10000) doesn't underflow if emin=-10.\n");
      mpfr_dump (a);
      mpfr_dump (b);
      exit (1);
    }
  mpfr_set_emin (emin);
  mpfr_clear (a);
  mpfr_clear (b);
}

static void
check_inexact (mp_prec_t p)
{
  mpfr_t x, y, z, t;
  unsigned long u;
  mp_prec_t q;
  int inexact, cmp;
  int rnd;

  mpfr_init2 (x, p);
  mpfr_init (y);
  mpfr_init (z);
  mpfr_init (t);
  mpfr_random (x);
  u = randlimb () % 2;
  for (q = 2; q <= p; q++)
    for (rnd = 0; rnd < GMP_RND_MAX; rnd++)
      {
        mpfr_set_prec (y, q);
        mpfr_set_prec (z, q + 10);
        mpfr_set_prec (t, q);
        inexact = mpfr_pow_ui (y, x, u, (mp_rnd_t) rnd);
        cmp = mpfr_pow_ui (z, x, u, (mp_rnd_t) rnd);
        if (mpfr_can_round (z, q + 10, (mp_rnd_t) rnd, (mp_rnd_t) rnd, q))
          {
            cmp = mpfr_set (t, z, (mp_rnd_t) rnd) || cmp;
            if (mpfr_cmp (y, t))
              {
                printf ("results differ for u=%lu rnd=%s\n",
                        u, mpfr_print_rnd_mode ((mp_rnd_t) rnd));
                printf ("x="); mpfr_print_binary (x); puts ("");
                printf ("y="); mpfr_print_binary (y); puts ("");
                printf ("t="); mpfr_print_binary (t); puts ("");
                printf ("z="); mpfr_print_binary (z); puts ("");
                exit (1);
              }
            if (((inexact == 0) && (cmp != 0)) ||
                ((inexact != 0) && (cmp == 0)))
              {
                printf ("Wrong inexact flag for p=%u, q=%u, rnd=%s\n",
                        (unsigned int) p, (unsigned int) q,
                        mpfr_print_rnd_mode ((mp_rnd_t) rnd));
                printf ("expected %d, got %d\n", cmp, inexact);
                printf ("u=%lu x=", u); mpfr_print_binary (x); puts ("");
                printf ("y="); mpfr_print_binary (y); puts ("");
                exit (1);
              }
          }
      }

  /* check exact power */
  mpfr_set_prec (x, p);
  mpfr_set_prec (y, p);
  mpfr_set_prec (z, p);
  mpfr_set_ui (x, 4, GMP_RNDN);
  mpfr_set_str (y, "0.5", 10, GMP_RNDN);
  test_pow (z, x, y, GMP_RNDZ);

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
  mpfr_clear (t);
}

static void
special ()
{
  mpfr_t x, y, z, t;

  mpfr_init2 (x, 53);
  mpfr_init2 (y, 53);
  mpfr_init2 (z, 53);
  mpfr_init2 (t, 2);

  mpfr_set_ui (x, 2, GMP_RNDN);
  mpfr_pow_si (x, x, -2, GMP_RNDN);
  if (mpfr_cmp_ui_2exp (x, 1, -2))
    {
      printf ("Error in pow_si(x,x,-2) for x=2\n");
      exit (1);
    }
  mpfr_set_ui (x, 2, GMP_RNDN);
  mpfr_set_si (y, -2, GMP_RNDN);
  test_pow (x, x, y, GMP_RNDN);
  if (mpfr_cmp_ui_2exp (x, 1, -2))
    {
      printf ("Error in pow(x,x,y) for x=2, y=-2\n");
      exit (1);
    }

  mpfr_set_prec (x, 64);
  mpfr_set_prec (y, 64);
  mpfr_set_prec (z, 64);
  mpfr_set_prec (t, 64);
  mpfr_set_str_binary (x, "0.111011000111100000111010000101010100110011010000011");
  mpfr_set_str_binary (y, "0.111110010100110000011101100011010111000010000100101");
  mpfr_set_str_binary (t, "0.1110110011110110001000110100100001001111010011111000010000011001");

  test_pow (z, x, y, GMP_RNDN);
  if (mpfr_cmp (z, t))
    {
      printf ("Error in mpfr_pow for prec=64, rnd=GMP_RNDN\n");
      exit (1);
    }

  mpfr_set_prec (x, 53);
  mpfr_set_prec (y, 53);
  mpfr_set_prec (z, 53);
  mpfr_set_str (x, "5.68824667828621954868e-01", 10, GMP_RNDN);
  mpfr_set_str (y, "9.03327850535952658895e-01", 10, GMP_RNDN);
  test_pow (z, x, y, GMP_RNDZ);
  if (mpfr_cmp_d(z, 0.60071044650456473235))
    {
      printf ("Error in mpfr_pow for prec=53, rnd=GMP_RNDZ\n");
      exit (1);
    }

  mpfr_set_prec (t, 2);
  mpfr_set_prec (x, 30);
  mpfr_set_prec (y, 30);
  mpfr_set_prec (z, 30);
  mpfr_set_str (x, "1.00000000001010111110001111011e1", 2, GMP_RNDN);
  mpfr_set_str (t, "-0.5", 10, GMP_RNDN);
  test_pow (z, x, t, GMP_RNDN);
  mpfr_set_str (y, "1.01101001111010101110000101111e-1", 2, GMP_RNDN);
  if (mpfr_cmp (z, y))
    {
      printf ("Error in mpfr_pow for prec=30, rnd=GMP_RNDN\n");
      exit (1);
    }

  mpfr_set_prec (x, 21);
  mpfr_set_prec (y, 21);
  mpfr_set_prec (z, 21);
  mpfr_set_str (x, "1.11111100100001100101", 2, GMP_RNDN);
  test_pow (z, x, t, GMP_RNDZ);
  mpfr_set_str (y, "1.01101011010001100000e-1", 2, GMP_RNDN);
  if (mpfr_cmp (z, y))
    {
      printf ("Error in mpfr_pow for prec=21, rnd=GMP_RNDZ\n");
      printf ("Expected ");
      mpfr_out_str (stdout, 2, 0, y, GMP_RNDN);
      printf ("\nGot      ");
      mpfr_out_str (stdout, 2, 0, z, GMP_RNDN);
      printf ("\n");
      exit (1);
    }

  mpfr_set_inf (x, 1);
  mpfr_set_prec (y, 2);
  mpfr_set_str_binary (y, "1E10");
  test_pow (z, x, y, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (z) && MPFR_IS_POS(z));
  mpfr_set_inf (x, -1);
  test_pow (z, x, y, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (z) && MPFR_IS_POS(z));
  mpfr_set_prec (y, 10);
  mpfr_set_str_binary (y, "1.000000001E9");
  test_pow (z, x, y, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (z) && MPFR_IS_NEG(z));
  mpfr_set_str_binary (y, "1.000000001E8");
  test_pow (z, x, y, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (z) && MPFR_IS_POS(z));

  mpfr_set_inf (x, -1);
  mpfr_set_prec (y, 2 * mp_bits_per_limb);
  mpfr_set_ui (y, 1, GMP_RNDN);
  mpfr_mul_2exp (y, y, mp_bits_per_limb - 1, GMP_RNDN);
  /* y = 2^(mp_bits_per_limb - 1) */
  test_pow (z, x, y, GMP_RNDN);
  MPFR_ASSERTN(mpfr_inf_p (z) && MPFR_IS_POS(z));
  mpfr_nextabove (y);
  test_pow (z, x, y, GMP_RNDN);
  /* y = 2^(mp_bits_per_limb - 1) + epsilon */
  MPFR_ASSERTN(mpfr_inf_p (z) && MPFR_IS_POS(z));
  mpfr_nextbelow (y);
  mpfr_div_2exp (y, y, 1, GMP_RNDN);
  mpfr_nextabove (y);
  test_pow (z, x, y, GMP_RNDN);
  /* y = 2^(mp_bits_per_limb - 2) + epsilon */
  MPFR_ASSERTN(mpfr_inf_p (z) && MPFR_IS_POS(z));

  mpfr_set_si (x, -1, GMP_RNDN);
  mpfr_set_prec (y, 2);
  mpfr_set_str_binary (y, "1E10");
  test_pow (z, x, y, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (z, 1) == 0);

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
  mpfr_clear (t);
}

void f(void);

static void
particular_cases (void)
{
  mpfr_t t[11], r;
  static const char *name[11] = {
    "NaN", "+inf", "-inf", "+0", "-0", "+1", "-1", "+2", "-2", "+0.5", "-0.5"};
  int i, j;
  int error = 0;

  for (i = 0; i < 11; i++)
    mpfr_init2 (t[i], 2);
  mpfr_init2 (r, 6);

  mpfr_set_nan (t[0]);
  mpfr_set_inf (t[1], 1);
  mpfr_set_ui (t[3], 0, GMP_RNDN);
  mpfr_set_ui (t[5], 1, GMP_RNDN);
  mpfr_set_ui (t[7], 2, GMP_RNDN);
  mpfr_div_2ui (t[9], t[5], 1, GMP_RNDN);
  for (i = 1; i < 11; i += 2)
    mpfr_neg (t[i+1], t[i], GMP_RNDN);

  for (i = 0; i < 11; i++)
    for (j = 0; j < 11; j++)
      {
        int p;
        static int q[11][11] = {
          /*          NaN +inf -inf  +0   -0   +1   -1   +2   -2  +0.5 -0.5 */
          /*  NaN */ { 0,   0,   0,  128, 128,  0,   0,   0,   0,   0,   0  },
          /* +inf */ { 0,   1,   2,  128, 128,  1,   2,   1,   2,   1,   2  },
          /* -inf */ { 0,   1,   2,  128, 128, -1,  -2,   1,   2,   1,   2  },
          /*  +0  */ { 0,   2,   1,  128, 128,  2,   1,   2,   1,   2,   1  },
          /*  -0  */ { 0,   2,   1,  128, 128, -2,  -1,   2,   1,   2,   1  },
          /*  +1  */ {128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128 },
          /*  -1  */ { 0,  128, 128, 128, 128,-128,-128, 128, 128,  0,   0  },
          /*  +2  */ { 0,   1,   2,  128, 128, 256,  64, 512,  32, 180,  90 },
          /*  -2  */ { 0,   1,   2,  128, 128,-256, -64, 512,  32,  0,   0  },
          /* +0.5 */ { 0,   2,   1,  128, 128,  64, 256,  32, 512,  90, 180 },
          /* -0.5 */ { 0,   2,   1,  128, 128, -64,-256,  32, 512,  0,   0  }
        };
        if (i == 5 && j == 1)
          f();

        test_pow (r, t[i], t[j], GMP_RNDN);
        p = mpfr_nan_p (r) ? 0 : mpfr_inf_p (r) ? 1 :
          mpfr_cmp_ui (r, 0) == 0 ? 2 :
          (int) (fabs (mpfr_get_d (r, GMP_RNDN)) * 128.0);
        if (p != 0 && MPFR_SIGN(r) < 0)
          p = -p;
        if (p != q[i][j])
          {
            printf ("Error in mpfr_pow for particular case (%s)^(%s) (%d,%d):\n"
                    "got %d instead of %d\n", name[i], name[j], i,j,p, q[i][j]);
            mpfr_dump (r);
            error = 1;
          }
      }

  for (i = 0; i < 11; i++)
    mpfr_clear (t[i]);
  mpfr_clear (r);

  if (error)
    exit (1);
}

void f (void)
{
}

static void
underflows (void)
{
  mpfr_t x, y;
  int i;

  mpfr_init2 (x, 64);
  mpfr_init2 (y, 64);

  mpfr_set_ui (x, 1, GMP_RNDN);
  mpfr_set_exp (x, mpfr_get_emin());

  for (i = 3; i < 10; i++)
    {
      mpfr_set_ui (y, i, GMP_RNDN);
      mpfr_div_2ui (y, y, 1, GMP_RNDN);
      test_pow (y, x, y, GMP_RNDN);
      if (!MPFR_IS_FP(y) || mpfr_cmp_ui (y, 0))
        {
          printf ("Error in mpfr_pow for ");
          mpfr_out_str (stdout, 2, 0, x, GMP_RNDN);
          printf (" ^ (%d/2)\nGot ", i);
          mpfr_out_str (stdout, 2, 0, y, GMP_RNDN);
          printf (" instead of 0.\n");
          exit (1);
        }
    }

  mpfr_clear (x);
  mpfr_clear (y);
}

static void
overflows (void)
{
  mpfr_t a, b;

  /* bug found by Ming J. Tsai <mingjt@delvron.us>, 4 Oct 2003 */

  mpfr_init_set_str (a, "5.1e32", 10, GMP_RNDN);
  mpfr_init (b);

  test_pow (b, a, a, GMP_RNDN);
  if (!(mpfr_inf_p (b) && mpfr_sgn (b) > 0))
    {
      printf ("Error for a^a for a=5.1e32\n");
      printf ("Expected +Inf, got ");
      mpfr_out_str (stdout, 10, 0, b, GMP_RNDN);
      printf ("\n");
      exit (1);
    }

  mpfr_clear(a);
  mpfr_clear(b);
}

int
main (void)
{
  mp_prec_t p;

  tests_start_mpfr ();

  special ();

  particular_cases ();

  check_pow_ui ();

  check_pow_si ();

  check_special_pow_si ();

  for (p = 2; p < 100; p++)
    check_inexact (p);

  underflows ();

  overflows ();

  test_generic (2, 100, 100);

  tests_end_mpfr ();
  return 0;
}
