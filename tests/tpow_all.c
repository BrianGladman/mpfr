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

#include <stdlib.h>

#include "mpfr-test.h"

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
cmpres (int i, int j, int rnd, mpfr_srcptr z1, int inex1,
        mpfr_srcptr z2, int inex2, const char *s)
{
  if (MPFR_IS_NAN (z1) && MPFR_IS_NAN (z2))
    return;
  if (mpfr_equal_p (z1, z2) && SAME_SIGN (inex1, inex2))
    return;

  printf ("Error with %s\n", s);
  printf ("x = %s, y = %s, %s\n", val[i], val[j],
          mpfr_print_rnd_mode ((mp_rnd_t) rnd));
  printf ("Expected ");
  mpfr_out_str (stdout, 10, 0, z1, GMP_RNDN);
  printf (", inex = %d\n", SIGN (inex1));
  printf ("Got      ");
  mpfr_out_str (stdout, 10, 0, z2, GMP_RNDN);
  printf (", inex = %d\n", SIGN (inex2));
  exit (1);
}

static int
is_odd (mpfr_srcptr x)
{
  /* does not work for large integers */
  return mpfr_integer_p (x) && (mpfr_get_si (x, GMP_RNDN) & 1);
}

static void
tst (void)
{
  int sv = sizeof (val) / sizeof (*val);
  int i, j;
  int rnd;
  mpfr_t x, y, z1, z2;

  mpfr_inits2 (53, x, y, z1, z2, (mpfr_ptr) 0);

  for (i = 0; i < sv; i++)
    for (j = 0; j < sv; j++)
      RND_LOOP (rnd)
        {
          int exact, inex1, inex2;

          if (mpfr_set_str (x, val[i], 10, GMP_RNDN) ||
              mpfr_set_str (y, val[j], 10, GMP_RNDN))
            {
              printf ("internal error for (%d,%d,%d)\n", i, j, rnd);
              exit (1);
            }
          mpfr_clear_flags ();
          inex1 = mpfr_pow (z1, x, y, (mp_rnd_t) rnd);
          if (mpfr_underflow_p ())
            err ("got underflow", i, j, rnd, z1, inex1);
          if (mpfr_overflow_p ())
            err ("got overflow", i, j, rnd, z1, inex1);
          if (! MPFR_IS_NAN (z1) && mpfr_nanflag_p ())
            err ("got NaN flag without NaN value", i, j, rnd, z1, inex1);
          if (MPFR_IS_NAN (z1) && ! mpfr_nanflag_p ())
            err ("got NaN value without NaN flag", i, j, rnd, z1, inex1);
          if (inex1 != 0 && ! mpfr_inexflag_p ())
            err ("got non-zero ternary value without inexact flag",
                 i, j, rnd, z1, inex1);
          if (inex1 == 0 && mpfr_inexflag_p ())
            err ("got null ternary value with inexact flag",
                 i, j, rnd, z1, inex1);
          exact = MPFR_IS_SINGULAR (z1) ||
            (mpfr_mul_2ui (z2, z1, 16, GMP_RNDN), mpfr_integer_p (z2));
          if (exact && inex1 != 0)
            err ("got exact value with ternary flag different from 0",
                 i, j, rnd, z1, inex1);
          if (! exact && inex1 == 0)
            err ("got inexact value with ternary flag equal to 0",
                 i, j, rnd, z1, inex1);
          if (MPFR_IS_ZERO (x) && ! MPFR_IS_NAN (y) && MPFR_NOTZERO (y))
            {
              if (MPFR_IS_NEG (y) && ! MPFR_IS_INF (z1))
                err ("expected an infinity", i, j, rnd, z1, inex1);
              if (MPFR_IS_POS (y) && ! MPFR_IS_ZERO (z1))
                err ("expected a zero", i, j, rnd, z1, inex1);
              if ((MPFR_IS_NEG (x) && is_odd (y)) ^ MPFR_IS_NEG (z1))
                err ("wrong sign", i, j, rnd, z1, inex1);
            }
          if (! MPFR_IS_NAN (x) && mpfr_cmp_si (x, -1) == 0)
            {
              /* x = -1 */
              if (! (MPFR_IS_INF (y) || mpfr_integer_p (y)) &&
                  ! MPFR_IS_NAN (z1))
                err ("expected NaN", i, j, rnd, z1, inex1);
              if ((MPFR_IS_INF (y) || (mpfr_integer_p (y) && ! is_odd (y)))
                  && ! mpfr_equal_p (z1, __gmpfr_one))
                err ("expected 1", i, j, rnd, z1, inex1);
              if (is_odd (y) &&
                  (MPFR_IS_NAN (z1) || mpfr_cmp_si (z1, -1) != 0))
                err ("expected -1", i, j, rnd, z1, inex1);
            }
          if ((mpfr_equal_p (x, __gmpfr_one) || MPFR_IS_ZERO (y)) &&
              ! mpfr_equal_p (z1, __gmpfr_one))
            err ("expected 1", i, j, rnd, z1, inex1);
          if (MPFR_IS_PURE_FP (x) && MPFR_IS_NEG (x) &&
              MPFR_IS_FP (y) && ! mpfr_integer_p (y) &&
              ! MPFR_IS_NAN (z1))
            err ("expected NaN", i, j, rnd, z1, inex1);
          if (MPFR_IS_INF (y) && MPFR_NOTZERO (x))
            {
              int cmpabs1 = mpfr_cmpabs (x, __gmpfr_one);

              if ((MPFR_IS_NEG (y) ? (cmpabs1 < 0) : (cmpabs1 > 0)) &&
                  ! (MPFR_IS_POS (z1) && MPFR_IS_INF (z1)))
                err ("expected +Inf", i, j, rnd, z1, inex1);
              if ((MPFR_IS_NEG (y) ? (cmpabs1 > 0) : (cmpabs1 < 0)) &&
                  ! (MPFR_IS_POS (z1) && MPFR_IS_ZERO (z1)))
                err ("expected +0", i, j, rnd, z1, inex1);
            }
          if (MPFR_IS_INF (x) && ! MPFR_IS_NAN (y) && MPFR_NOTZERO (y))
            {
              if (MPFR_IS_POS (y) && ! MPFR_IS_INF (z1))
                err ("expected an infinity", i, j, rnd, z1, inex1);
              if (MPFR_IS_NEG (y) && ! MPFR_IS_ZERO (z1))
                err ("expected a zero", i, j, rnd, z1, inex1);
              if ((MPFR_IS_NEG (x) && is_odd (y)) ^ MPFR_IS_NEG (z1))
                err ("wrong sign", i, j, rnd, z1, inex1);
            }

          __gmpfr_flags = MPFR_FLAGS_ALL;
          inex2 = mpfr_pow (z2, x, y, (mp_rnd_t) rnd);
          cmpres (i, j, rnd, z1, inex1, z2, inex2, "mpfr_pow, flags set");

          if (mpfr_integer_p (y))
            {
              long yy = mpfr_get_si (y, GMP_RNDN);

              /* If y >= 0 and y is not -0, we can test mpfr_pow_ui. */
              if (MPFR_IS_POS (y))
                {
                  MPFR_ASSERTN (yy >= 0);
                  mpfr_clear_flags ();
                  inex2 = mpfr_pow_ui (z2, x, yy, (mp_rnd_t) rnd);
                  cmpres (i, j, rnd, z1, inex1, z2, inex2,
                          "mpfr_pow_ui, flags cleared");
                  __gmpfr_flags = MPFR_FLAGS_ALL;
                  inex2 = mpfr_pow_ui (z2, x, yy, (mp_rnd_t) rnd);
                  cmpres (i, j, rnd, z1, inex1, z2, inex2,
                          "mpfr_pow_ui, flags set");

                  /* If x >= 0 and x is not -0, we can test mpfr_ui_pow_ui. */
                  if (mpfr_integer_p (x) && MPFR_IS_POS (x))
                    {
                      unsigned long xx = mpfr_get_ui (x, GMP_RNDN);

                      mpfr_clear_flags ();
                      inex2 = mpfr_ui_pow_ui (z2, xx, yy, (mp_rnd_t) rnd);
                      cmpres (i, j, rnd, z1, inex1, z2, inex2,
                              "mpfr_ui_pow_ui, flags cleared");
                      __gmpfr_flags = MPFR_FLAGS_ALL;
                      inex2 = mpfr_ui_pow_ui (z2, xx, yy, (mp_rnd_t) rnd);
                      cmpres (i, j, rnd, z1, inex1, z2, inex2,
                              "mpfr_ui_pow_ui, flags set");
                    }
                }

              /* We can test mpfr_pow_si and mpfr_pow_z when y is not -0. */
              if (MPFR_IS_POS (y) || MPFR_NOTZERO (y))
                {
                  mpz_t yyy;

                  mpfr_clear_flags ();
                  inex2 = mpfr_pow_si (z2, x, yy, (mp_rnd_t) rnd);
                  cmpres (i, j, rnd, z1, inex1, z2, inex2,
                          "mpfr_pow_si, flags cleared");
                  __gmpfr_flags = MPFR_FLAGS_ALL;
                  inex2 = mpfr_pow_si (z2, x, yy, (mp_rnd_t) rnd);
                  cmpres (i, j, rnd, z1, inex1, z2, inex2,
                          "mpfr_pow_si, flags set");

                  mpz_init (yyy);
                  mpfr_get_z (yyy, y, GMP_RNDN);
                  mpfr_clear_flags ();
                  inex2 = mpfr_pow_z (z2, x, yyy, (mp_rnd_t) rnd);
                  cmpres (i, j, rnd, z1, inex1, z2, inex2,
                          "mpfr_pow_z, flags cleared");
                  __gmpfr_flags = MPFR_FLAGS_ALL;
                  inex2 = mpfr_pow_z (z2, x, yyy, (mp_rnd_t) rnd);
                  cmpres (i, j, rnd, z1, inex1, z2, inex2,
                          "mpfr_pow_z, flags set");
                  mpz_clear (yyy);
                }
            }

          /* If x >= 0 and x is not -0, we can test mpfr_ui_pow. */
          if (mpfr_integer_p (x) && MPFR_IS_POS (x))
            {
              unsigned long xx = mpfr_get_ui (x, GMP_RNDN);

              mpfr_clear_flags ();
              inex2 = mpfr_ui_pow (z2, xx, y, (mp_rnd_t) rnd);
              cmpres (i, j, rnd, z1, inex1, z2, inex2,
                      "mpfr_ui_pow, flags cleared");
              __gmpfr_flags = MPFR_FLAGS_ALL;
              inex2 = mpfr_ui_pow (z2, xx, y, (mp_rnd_t) rnd);
              cmpres (i, j, rnd, z1, inex1, z2, inex2,
                      "mpfr_ui_pow, flags set");
            }
        }
  mpfr_clears (x, y, z1, z2, (mpfr_ptr) 0);
}

int
main (void)
{
  tests_start_mpfr ();
  tst ();
  tests_end_mpfr ();
  return 0;
}
