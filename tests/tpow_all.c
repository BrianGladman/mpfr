/* Test file for the various power functions

Copyright 2008 Free Software Foundation, Inc.
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

/* Note: some tests of the other tpow* test files could be moved there.
   The main goal of this test file is to test _all_ the power functions
   on special values, to make sure that they are consistent and give the
   expected result, in particular because such special cases are handled
   in different ways in each function. */

/* Execute with at least an argument to report all the errors found by
   comparisons. */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "mpfr-test.h"

/* Behavior of cmpres (called by test_others):
 *   0: stop as soon as an error is found.
 *   1: report all errors found by test_others.
 *  -1: the 1 is changed to this value as soon as an error has been found.
 */
static int all_cmpres_errors;

static char *val[] =
  { "@NaN@", "-@Inf@", "-4", "-3", "-2", "-1.5", "-1", "-0.5", "-0",
    "0", "0.5", "1", "1.5", "2", "3", "4", "@Inf@" };

static void
err (const char *s, int i, int j, int rnd, mpfr_srcptr z, int inex)
{
  puts (s);
  printf ("x = %s, y = %s, %s\n", val[i], val[j],
          mpfr_print_rnd_mode ((mp_rnd_t) rnd));
  printf ("z = ");
  mpfr_out_str (stdout, 10, 0, z, GMP_RNDN);
  printf ("\ninex = %d\n", inex);
  exit (1);
}

static void
cmpres (int spx, const void *px, const char *sy, mp_rnd_t rnd,
        mpfr_srcptr z1, int inex1, mpfr_srcptr z2, int inex2,
        const char *s)
{
  if (MPFR_IS_NAN (z1) && MPFR_IS_NAN (z2))
    return;
  if (mpfr_equal_p (z1, z2) && SAME_SIGN (inex1, inex2))
    return;

  printf ("Error with %s\nx = ", s);
  if (spx)
    printf ("%s, ", (char *) px);
  else
    {
      mpfr_out_str (stdout, 16, 0, (mpfr_ptr) px, GMP_RNDN);
      puts (",");
    }
  printf ("y = %s, %s\n", sy, mpfr_print_rnd_mode (rnd));
  printf ("Expected ");
  mpfr_out_str (stdout, 16, 0, z1, GMP_RNDN);
  printf (", inex = %d\n", SIGN (inex1));
  printf ("Got      ");
  mpfr_out_str (stdout, 16, 0, z2, GMP_RNDN);
  printf (", inex = %d\n", SIGN (inex2));
  if (all_cmpres_errors != 0)
    all_cmpres_errors = -1;
  else
    exit (1);
}

static int
is_odd (mpfr_srcptr x)
{
  /* does not work for large integers */
  return mpfr_integer_p (x) && (mpfr_get_si (x, GMP_RNDN) & 1);
}

/* Compare the result (z1,inex1) of mpfr_pow with all flags cleared
   with those of mpfr_pow with all flags set and of the other power
   functions. Arguments x and y are the input values; sx and sy are
   their string representations (sx may be null); rnd contains the
   rounding mode. */
static void
test_others (const void *sx, const char *sy, mp_rnd_t rnd,
             mpfr_srcptr x, mpfr_srcptr y, mpfr_srcptr z1, int inex1)
{
  mpfr_t z2;
  int inex2;
  int spx = sx != NULL;

  if (!spx)
    sx = x;

  mpfr_init2 (z2, mpfr_get_prec (z1));

  __gmpfr_flags = MPFR_FLAGS_ALL;
  inex2 = mpfr_pow (z2, x, y, rnd);
  cmpres (spx, sx, sy, rnd, z1, inex1, z2, inex2, "mpfr_pow, flags set");

  /* If y is an integer that fits in an unsigned long and is not -0,
     we can test mpfr_pow_ui. */
  if (MPFR_IS_POS (y) && mpfr_integer_p (y) &&
      mpfr_fits_ulong_p (y, GMP_RNDN))
    {
      unsigned long yy = mpfr_get_ui (y, GMP_RNDN);

      mpfr_clear_flags ();
      inex2 = mpfr_pow_ui (z2, x, yy, rnd);
      cmpres (spx, sx, sy, rnd, z1, inex1, z2, inex2,
              "mpfr_pow_ui, flags cleared");
      __gmpfr_flags = MPFR_FLAGS_ALL;
      inex2 = mpfr_pow_ui (z2, x, yy, rnd);
      cmpres (spx, sx, sy, rnd, z1, inex1, z2, inex2,
              "mpfr_pow_ui, flags set");

      /* If x is an integer that fits in an unsigned long and is not -0,
         we can also test mpfr_ui_pow_ui. */
      if (MPFR_IS_POS (x) && mpfr_integer_p (x) &&
          mpfr_fits_ulong_p (x, GMP_RNDN))
        {
          unsigned long xx = mpfr_get_ui (x, GMP_RNDN);

          mpfr_clear_flags ();
          inex2 = mpfr_ui_pow_ui (z2, xx, yy, rnd);
          cmpres (spx, sx, sy, rnd, z1, inex1, z2, inex2,
                  "mpfr_ui_pow_ui, flags cleared");
          __gmpfr_flags = MPFR_FLAGS_ALL;
          inex2 = mpfr_ui_pow_ui (z2, xx, yy, rnd);
          cmpres (spx, sx, sy, rnd, z1, inex1, z2, inex2,
                  "mpfr_ui_pow_ui, flags set");
        }
    }

  /* If y is an integer but not -0, we can test mpfr_pow_z, and
     possibly mpfr_pow_si. */
  if ((MPFR_IS_POS (y) || MPFR_NOTZERO (y)) && mpfr_integer_p (y))
    {
      mpz_t yyy;

      /* If y fits in a long, we can test mpfr_pow_si. */
      if (mpfr_fits_slong_p (y, GMP_RNDN))
        {
          long yy = mpfr_get_si (y, GMP_RNDN);

          mpfr_clear_flags ();
          inex2 = mpfr_pow_si (z2, x, yy, rnd);
          cmpres (spx, sx, sy, rnd, z1, inex1, z2, inex2,
                  "mpfr_pow_si, flags cleared");
          __gmpfr_flags = MPFR_FLAGS_ALL;
          inex2 = mpfr_pow_si (z2, x, yy, rnd);
          cmpres (spx, sx, sy, rnd, z1, inex1, z2, inex2,
                  "mpfr_pow_si, flags set");
        }

      /* Test mpfr_pow_z. */
      mpz_init (yyy);
      mpfr_get_z (yyy, y, GMP_RNDN);
      mpfr_clear_flags ();
      inex2 = mpfr_pow_z (z2, x, yyy, rnd);
      cmpres (spx, sx, sy, rnd, z1, inex1, z2, inex2,
              "mpfr_pow_z, flags cleared");
      __gmpfr_flags = MPFR_FLAGS_ALL;
      inex2 = mpfr_pow_z (z2, x, yyy, rnd);
      cmpres (spx, sx, sy, rnd, z1, inex1, z2, inex2,
              "mpfr_pow_z, flags set");
      mpz_clear (yyy);
    }

  /* If x is an integer that fits in an unsigned long and is not -0,
     we can test mpfr_ui_pow. */
  if (MPFR_IS_POS (x) && mpfr_integer_p (x) &&
      mpfr_fits_ulong_p (x, GMP_RNDN))
    {
      unsigned long xx = mpfr_get_ui (x, GMP_RNDN);

      mpfr_clear_flags ();
      inex2 = mpfr_ui_pow (z2, xx, y, rnd);
      cmpres (spx, sx, sy, rnd, z1, inex1, z2, inex2,
              "mpfr_ui_pow, flags cleared");
      __gmpfr_flags = MPFR_FLAGS_ALL;
      inex2 = mpfr_ui_pow (z2, xx, y, rnd);
      cmpres (spx, sx, sy, rnd, z1, inex1, z2, inex2,
              "mpfr_ui_pow, flags set");

      /* If x = 2, we can test mpfr_exp2. */
      if (xx == 2)
        {
          mpfr_clear_flags ();
          inex2 = mpfr_exp2 (z2, y, rnd);
          cmpres (spx, sx, sy, rnd, z1, inex1, z2, inex2,
                  "mpfr_exp2, flags cleared");
          __gmpfr_flags = MPFR_FLAGS_ALL;
          inex2 = mpfr_exp2 (z2, y, rnd);
          cmpres (spx, sx, sy, rnd, z1, inex1, z2, inex2,
                  "mpfr_exp2, flags set");
        }

      /* If x = 10, we can test mpfr_exp10. */
      if (xx == 10)
        {
          mpfr_clear_flags ();
          inex2 = mpfr_exp10 (z2, y, rnd);
          cmpres (spx, sx, sy, rnd, z1, inex1, z2, inex2,
                  "mpfr_exp10, flags cleared");
          __gmpfr_flags = MPFR_FLAGS_ALL;
          inex2 = mpfr_exp10 (z2, y, rnd);
          cmpres (spx, sx, sy, rnd, z1, inex1, z2, inex2,
                  "mpfr_exp10, flags set");
        }
    }

  mpfr_clear (z2);
}

static void
tst (void)
{
  int sv = sizeof (val) / sizeof (*val);
  int i, j;
  int rnd;
  mpfr_t x, y, z, tmp;

  mpfr_inits2 (53, x, y, z, tmp, (mpfr_ptr) 0);

  for (i = 0; i < sv; i++)
    for (j = 0; j < sv; j++)
      RND_LOOP (rnd)
        {
          int exact, inex;

          if (mpfr_set_str (x, val[i], 10, GMP_RNDN) ||
              mpfr_set_str (y, val[j], 10, GMP_RNDN))
            {
              printf ("internal error for (%d,%d,%d)\n", i, j, rnd);
              exit (1);
            }
          mpfr_clear_flags ();
          inex = mpfr_pow (z, x, y, (mp_rnd_t) rnd);
          if (mpfr_underflow_p ())
            err ("got underflow", i, j, rnd, z, inex);
          if (mpfr_overflow_p ())
            err ("got overflow", i, j, rnd, z, inex);
          if (! MPFR_IS_NAN (z) && mpfr_nanflag_p ())
            err ("got NaN flag without NaN value", i, j, rnd, z, inex);
          if (MPFR_IS_NAN (z) && ! mpfr_nanflag_p ())
            err ("got NaN value without NaN flag", i, j, rnd, z, inex);
          if (inex != 0 && ! mpfr_inexflag_p ())
            err ("got non-zero ternary value without inexact flag",
                 i, j, rnd, z, inex);
          if (inex == 0 && mpfr_inexflag_p ())
            err ("got null ternary value with inexact flag",
                 i, j, rnd, z, inex);
          exact = MPFR_IS_SINGULAR (z) ||
            (mpfr_mul_2ui (tmp, z, 16, GMP_RNDN), mpfr_integer_p (tmp));
          if (exact && inex != 0)
            err ("got exact value with ternary flag different from 0",
                 i, j, rnd, z, inex);
          if (! exact && inex == 0)
            err ("got inexact value with ternary flag equal to 0",
                 i, j, rnd, z, inex);
          if (MPFR_IS_ZERO (x) && ! MPFR_IS_NAN (y) && MPFR_NOTZERO (y))
            {
              if (MPFR_IS_NEG (y) && ! MPFR_IS_INF (z))
                err ("expected an infinity", i, j, rnd, z, inex);
              if (MPFR_IS_POS (y) && ! MPFR_IS_ZERO (z))
                err ("expected a zero", i, j, rnd, z, inex);
              if ((MPFR_IS_NEG (x) && is_odd (y)) ^ MPFR_IS_NEG (z))
                err ("wrong sign", i, j, rnd, z, inex);
            }
          if (! MPFR_IS_NAN (x) && mpfr_cmp_si (x, -1) == 0)
            {
              /* x = -1 */
              if (! (MPFR_IS_INF (y) || mpfr_integer_p (y)) &&
                  ! MPFR_IS_NAN (z))
                err ("expected NaN", i, j, rnd, z, inex);
              if ((MPFR_IS_INF (y) || (mpfr_integer_p (y) && ! is_odd (y)))
                  && ! mpfr_equal_p (z, __gmpfr_one))
                err ("expected 1", i, j, rnd, z, inex);
              if (is_odd (y) &&
                  (MPFR_IS_NAN (z) || mpfr_cmp_si (z, -1) != 0))
                err ("expected -1", i, j, rnd, z, inex);
            }
          if ((mpfr_equal_p (x, __gmpfr_one) || MPFR_IS_ZERO (y)) &&
              ! mpfr_equal_p (z, __gmpfr_one))
            err ("expected 1", i, j, rnd, z, inex);
          if (MPFR_IS_PURE_FP (x) && MPFR_IS_NEG (x) &&
              MPFR_IS_FP (y) && ! mpfr_integer_p (y) &&
              ! MPFR_IS_NAN (z))
            err ("expected NaN", i, j, rnd, z, inex);
          if (MPFR_IS_INF (y) && MPFR_NOTZERO (x))
            {
              int cmpabs1 = mpfr_cmpabs (x, __gmpfr_one);

              if ((MPFR_IS_NEG (y) ? (cmpabs1 < 0) : (cmpabs1 > 0)) &&
                  ! (MPFR_IS_POS (z) && MPFR_IS_INF (z)))
                err ("expected +Inf", i, j, rnd, z, inex);
              if ((MPFR_IS_NEG (y) ? (cmpabs1 > 0) : (cmpabs1 < 0)) &&
                  ! (MPFR_IS_POS (z) && MPFR_IS_ZERO (z)))
                err ("expected +0", i, j, rnd, z, inex);
            }
          if (MPFR_IS_INF (x) && ! MPFR_IS_NAN (y) && MPFR_NOTZERO (y))
            {
              if (MPFR_IS_POS (y) && ! MPFR_IS_INF (z))
                err ("expected an infinity", i, j, rnd, z, inex);
              if (MPFR_IS_NEG (y) && ! MPFR_IS_ZERO (z))
                err ("expected a zero", i, j, rnd, z, inex);
              if ((MPFR_IS_NEG (x) && is_odd (y)) ^ MPFR_IS_NEG (z))
                err ("wrong sign", i, j, rnd, z, inex);
            }
          test_others (val[i], val[j], (mp_rnd_t) rnd, x, y, z, inex);
        }
  mpfr_clears (x, y, z, tmp, (mpfr_ptr) 0);
}

static void
underflow_up (int extended_emin)
{
  mpfr_t x, y, z, z0, eps;
  mp_exp_t n;
  int inex;
  int rnd;

  n = 1 - mpfr_get_emin ();
  MPFR_ASSERTN (n > 1);
  if (n > ULONG_MAX)
    return;

  mpfr_init2 (eps, 2);
  mpfr_set_ui_2exp (eps, 1, -1, GMP_RNDN);  /* 1/2 */
  mpfr_div_ui (eps, eps, n, GMP_RNDZ);      /* 1/(2n) rounded toward zero */

  mpfr_init2 (x, sizeof (unsigned long) * CHAR_BIT + 1);
  inex = mpfr_ui_sub (x, 1, eps, GMP_RNDN);
  MPFR_ASSERTN (inex == 0);  /* since n < 2^(size_of_long_in_bits) */
  inex = mpfr_div_2ui (x, x, 1, GMP_RNDN);  /* 1/2 - eps/2 exactly */
  MPFR_ASSERTN (inex == 0);

  mpfr_init2 (y, sizeof (unsigned long) * CHAR_BIT);
  inex = mpfr_set_ui (y, n, GMP_RNDN);
  MPFR_ASSERTN (inex == 0);

  /* 0 < eps < 1 / (2n), thus (1 - eps)^n > 1/2,
     and 1/2 (1/2)^n < (1/2 - eps/2)^n < (1/2)^n. */
  mpfr_inits2 (64, z, z0, (mpfr_ptr) 0);
  RND_LOOP (rnd)
    {
      unsigned int ufinex = MPFR_FLAGS_UNDERFLOW | MPFR_FLAGS_INEXACT;
      int expected_inex;
      char sy[256];

      mpfr_clear_flags ();
      inex = mpfr_pow (z, x, y, (mp_rnd_t) rnd);
      if (__gmpfr_flags != ufinex)
        {
          printf ("Error in underflow_up for %s",
                  mpfr_print_rnd_mode ((mp_rnd_t) rnd));
          if (extended_emin)
            printf (" and extended emin");
          printf ("\n");
          printf ("got %u instead of %u\n", __gmpfr_flags, ufinex);
          exit (1);
        }
      mpfr_set_ui (z0, 0, GMP_RNDN);
      expected_inex = rnd == GMP_RNDN || rnd == GMP_RNDU ?
        (mpfr_nextabove (z0), 1) : -1;
      sprintf (sy, "%lu", (unsigned long) n);
      cmpres (0, x, sy, (mp_rnd_t) rnd, z0, expected_inex, z, inex,
              extended_emin ? "underflow_up and extended emin" :
              "underflow_up");
      test_others (NULL, sy, (mp_rnd_t) rnd, x, y, z, inex);
    }

  mpfr_clears (x, y, z, z0, eps, (mpfr_ptr) 0);
}

static void
underflow (void)
{
  mp_exp_t emin;

  underflow_up (0);

  emin = mpfr_get_emin ();
  set_emin (MPFR_EMIN_MIN);
  if (mpfr_get_emin () != emin)
    {
      underflow_up (1);
      set_emin (emin);
    }
}

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();
  all_cmpres_errors = argc > 1;
  tst ();
  underflow ();
  tests_end_mpfr ();
  return all_cmpres_errors < 0;
}
