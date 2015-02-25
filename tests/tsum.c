/* tsum -- test file for the list summation function

Copyright 2004-2015 Free Software Foundation, Inc.
Contributed by the AriC and Caramel projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

/* TODO: add some generic random test with cancellations. Something like:
   1. Generate random numbers with random precisions.
   2. Compute the sum s at some random precision and some rounding mode.
   3. While s != 0:
   4. Include -s in the array.
   5. Reorder the terms randomly.
   6. Recompute a new sum s' at some random precision and some rounding mode.
   7. Check that |s'| < ulp(s), with a factor 1/2 for MPFR_RNDN.
   8. Reiterate at (3) with s = s'.
   Also add tests with intermediate overflows and tests with underflows
   (this matters here as we don't have subnormals).
   Note: partly done in the cancel() test (r8857); the result is not tested
   yet, but this test already shows an efficiency problem.
*/

#include "mpfr-test.h"

#ifndef MPFR_NCANCEL
#define MPFR_NCANCEL 10
#endif

static mpfr_prec_t
get_prec_max (mpfr_t *t, int n)
{
  mpfr_exp_t e, min, max;
  int i;

  min = MPFR_EMAX_MAX;
  max = MPFR_EMIN_MIN;
  for (i = 0; i < n; i++)
    if (MPFR_IS_PURE_FP (t[i]))
      {
        e = MPFR_GET_EXP (t[i]);
        if (e > max)
          max = e;
        e -= MPFR_GET_PREC (t[i]);
        if (e < min)
          min = e;
      }

  return min > max ? MPFR_PREC_MIN /* no pure FP values */
    : max - min + __gmpfr_ceil_log2 (n);
}

static void
get_exact_sum (mpfr_t sum, mpfr_t *tab, int n)
{
  int i;

  mpfr_set_prec (sum, get_prec_max (tab, n));
  mpfr_set_ui (sum, 0, MPFR_RNDN);
  for (i = 0; i < n; i++)
    if (mpfr_add (sum, sum, tab[i], MPFR_RNDN))
      {
        printf ("FIXME: get_exact_sum is buggy.\n");
        exit (1);
      }
}

static void
generic_tests (void)
{
  mpfr_t exact_sum, sum1, sum2;
  mpfr_t *t;
  mpfr_ptr *p;
  mpfr_prec_t precmax = 444;
  int i, m, nmax = 500;
  int rnd_mode;

  t = (mpfr_t *) (*__gmp_allocate_func) (nmax * sizeof(mpfr_t));
  p = (mpfr_ptr *) (*__gmp_allocate_func) (nmax * sizeof(mpfr_ptr));
  for (i = 0; i < nmax; i++)
    {
      mpfr_init2 (t[i], precmax);
      p[i] = t[i];
    }
  mpfr_inits2 (precmax, exact_sum, sum1, sum2, (mpfr_ptr) 0);

  for (m = 1; m < 4000; m++)
    {
      int non_uniform, n;
      mpfr_prec_t prec;

      non_uniform = randlimb () % 10;
      n = (randlimb () % nmax) + 1;
      prec = MPFR_PREC_MIN + (randlimb () % (precmax - MPFR_PREC_MIN + 1));
      mpfr_set_prec (sum1, prec);
      mpfr_set_prec (sum2, prec);

      for (i = 0; i < n; i++)
        {
          mpfr_set_prec (t[i], MPFR_PREC_MIN +
                         (randlimb () % (precmax - MPFR_PREC_MIN + 1)));
          mpfr_urandomb (t[i], RANDS);
          if (non_uniform && MPFR_NOTZERO (t[i]))
            mpfr_set_exp (t[i], randlimb () % 1000);
        }
      get_exact_sum (exact_sum, t, n);
      RND_LOOP (rnd_mode)
        {
          int inex1, inex2;

          inex1 = mpfr_set (sum1, exact_sum, (mpfr_rnd_t) rnd_mode);
          inex2 = mpfr_sum (sum2, p, n, (mpfr_rnd_t) rnd_mode);
          if (! mpfr_equal_p (sum1, sum2) || inex1 != inex2)
            {
              printf ("generic_tests failed on m = %d, %s\n", m,
                      mpfr_print_rnd_mode ((mpfr_rnd_t) rnd_mode));
              printf ("Expected ");
              mpfr_dump (sum1);
              printf ("with inex = %d\n", inex1);
              printf ("Got      ");
              mpfr_dump (sum2);
              printf ("with inex = %d\n", inex2);
              exit (1);
            }
        }
    }

  for (i = 0; i < nmax; i++)
    mpfr_clear (t[i]);
  mpfr_clears (exact_sum, sum1, sum2, (mpfr_ptr) 0);
  (*__gmp_free_func) (t, nmax * sizeof(mpfr_t));
  (*__gmp_free_func) (p, nmax * sizeof(mpfr_ptr));
}

static
void check_special (void)
{
  mpfr_t tab[3], r;
  mpfr_ptr tabp[3];
  int i;

  mpfr_inits (tab[0], tab[1], tab[2], r, (mpfr_ptr) 0);
  tabp[0] = tab[0];
  tabp[1] = tab[1];
  tabp[2] = tab[2];

  i = mpfr_sum (r, tabp, 0, MPFR_RNDN);
  if (!MPFR_IS_ZERO (r) || !MPFR_IS_POS (r) || i != 0)
    {
      printf ("Special case n==0 failed!\n");
      exit (1);
    }

  mpfr_set_ui (tab[0], 42, MPFR_RNDN);
  i = mpfr_sum (r, tabp, 1, MPFR_RNDN);
  if (mpfr_cmp_ui (r, 42) || i != 0)
    {
      printf ("Special case n==1 failed!\n");
      exit (1);
    }

  mpfr_set_ui (tab[1], 17, MPFR_RNDN);
  MPFR_SET_NAN (tab[2]);
  i = mpfr_sum (r, tabp, 3, MPFR_RNDN);
  if (!MPFR_IS_NAN (r) || i != 0)
    {
      printf ("Special case NAN failed!\n");
      exit (1);
    }

  MPFR_SET_INF (tab[2]);
  MPFR_SET_POS (tab[2]);
  i = mpfr_sum (r, tabp, 3, MPFR_RNDN);
  if (!MPFR_IS_INF (r) || !MPFR_IS_POS (r) || i != 0)
    {
      printf ("Special case +INF failed!\n");
      exit (1);
    }

  MPFR_SET_INF (tab[2]);
  MPFR_SET_NEG (tab[2]);
  i = mpfr_sum (r, tabp, 3, MPFR_RNDN);
  if (!MPFR_IS_INF (r) || !MPFR_IS_NEG (r) || i != 0)
    {
      printf ("Special case -INF failed!\n");
      exit (1);
    }

  MPFR_SET_ZERO (tab[1]);
  i = mpfr_sum (r, tabp, 2, MPFR_RNDN);
  if (mpfr_cmp_ui (r, 42) || i != 0)
    {
      printf ("Special case 42+0 failed!\n");
      exit (1);
    }

  MPFR_SET_NAN (tab[0]);
  i = mpfr_sum (r, tabp, 3, MPFR_RNDN);
  if (!MPFR_IS_NAN (r) || i != 0)
    {
      printf ("Special case NAN+0+-INF failed!\n");
      exit (1);
    }

  mpfr_set_inf (tab[0], 1);
  mpfr_set_ui  (tab[1], 59, MPFR_RNDN);
  mpfr_set_inf (tab[2], -1);
  i = mpfr_sum (r, tabp, 3, MPFR_RNDN);
  if (!MPFR_IS_NAN (r) || i != 0)
    {
      printf ("Special case +INF + 59 +-INF failed!\n");
      exit (1);
    }

  mpfr_clears (tab[0], tab[1], tab[2], r, (mpfr_ptr) 0);
}

#define NC 7
#define NS 6

static void
check_more_special (void)
{
  char *str[NC] = { "NaN", "+Inf", "-Inf", "+0", "-0", "+1", "-1" };
  int i, r, k[NS];
  mpfr_t c[NC], s[NS], sum;
  mpfr_ptr p[NS];
  int inex;

  for (i = 0; i < NC; i++)
    {
      int ret;
      mpfr_init2 (c[i], 8);
      ret = mpfr_set_str (c[i], str[i], 0, MPFR_RNDN);
      MPFR_ASSERTN (ret == 0);
    }
  for (i = 0; i < NS; i++)
    mpfr_init2 (s[i], 8);
  mpfr_init2 (sum, 8);

  RND_LOOP(r)
    {
      i = 0;
      while (1)
        {
          while (i < NS)
            {
              p[i] = c[0];
              mpfr_set_nan (s[i]);
              k[i++] = 0;
            }
          inex = mpfr_sum (sum, p, NS, (mpfr_rnd_t) r);
          if (! ((MPFR_IS_NAN (sum) && MPFR_IS_NAN (s[NS-1])) ||
                 (mpfr_equal_p (sum, s[NS-1]) &&
                  MPFR_SIGN (sum) == MPFR_SIGN (s[NS-1]))) || inex != 0)
            {
              printf ("Error in check_more_special on %s",
                      mpfr_print_rnd_mode ((mpfr_rnd_t) r));
              for (i = 0; i < NS; i++)
                printf (" %d", k[i]);
              printf (" with\n");
              for (i = 0; i < NS; i++)
                {
                  printf ("  p[%d] = %s = ", i, str[k[i]]);
                  mpfr_dump (p[i]);
                }
              printf ("Expected ");
              mpfr_dump (s[NS-1]);
              printf ("with inex = 0\n");
              printf ("Got      ");
              mpfr_dump (sum);
              printf ("with inex = %d\n", inex);
              exit (1);
            }
          while (k[--i] == NC-1)
            if (i == 0)
              goto next_rnd;
          p[i] = c[++k[i]];
          if (i == 0)
            mpfr_set (s[i], p[i], (mpfr_rnd_t) r);
          else
            mpfr_add (s[i], s[i-1], p[i], (mpfr_rnd_t) r);
          i++;
        }
    next_rnd: ;
    }

  for (i = 0; i < NC; i++)
    mpfr_clear (c[i]);
  for (i = 0; i < NS; i++)
    mpfr_clear (s[i]);
  mpfr_clear (sum);
}

/* bug reported by Joseph S. Myers on 2013-10-27
   https://sympa.inria.fr/sympa/arc/mpfr/2013-10/msg00015.html */
static void
bug20131027 (void)
{
  mpfr_t sum, t[4];
  mpfr_ptr p[4];
  char *s[4] = {
    "0x1p1000",
    "-0x0.fffffffffffff80000000000000001p1000",
    "-0x1p947",
    "0x1p880"
  };
  int i, r;

  mpfr_init2 (sum, 53);

  for (i = 0; i < 4; i++)
    {
      mpfr_init2 (t[i], i == 0 ? 53 : 1000);
      mpfr_set_str (t[i], s[i], 0, MPFR_RNDN);
      p[i] = t[i];
    }

  RND_LOOP(r)
    {
      int expected_sign = (mpfr_rnd_t) r == MPFR_RNDD ? -1 : 1;
      int inex;

      inex = mpfr_sum (sum, p, 4, (mpfr_rnd_t) r);

      if (MPFR_NOTZERO (sum) || MPFR_SIGN (sum) != expected_sign || inex != 0)
        {
          printf ("mpfr_sum incorrect in bug20131027 for %s:\n"
                  "expected %c0 with inex = 0, got ",
                  mpfr_print_rnd_mode ((mpfr_rnd_t) r),
                  expected_sign > 0 ? '+' : '-');
          mpfr_dump (sum);
          printf ("with inex = %d\n", inex);
          exit (1);
        }
    }

  for (i = 0; i < 4; i++)
    mpfr_clear (t[i]);
  mpfr_clear (sum);
}

/* TODO: A test with more inputs (but can't be compared to mpfr_add). */
static void
check_extreme (void)
{
  mpfr_t u, v, w, x, y;
  mpfr_ptr t[2];
  int i, inex1, inex2, r;

  t[0] = u;
  t[1] = v;

  mpfr_inits2 (32, u, v, w, x, y, (mpfr_ptr) 0);
  mpfr_setmin (u, mpfr_get_emax ());
  mpfr_setmax (v, mpfr_get_emin ());
  mpfr_setmin (w, mpfr_get_emax () - 40);
  RND_LOOP (r)
    for (i = 0; i < 2; i++)
      {
        mpfr_set_prec (x, 64);
        inex1 = mpfr_add (x, u, w, MPFR_RNDN);
        MPFR_ASSERTN (inex1 == 0);
        inex1 = mpfr_prec_round (x, 32, (mpfr_rnd_t) r);
        inex2 = mpfr_sum (y, t, 2, (mpfr_rnd_t) r);
        if (! mpfr_equal_p (x, y) || inex1 != inex2)
          {
            printf ("Error in check_extreme (%s, i = %d)\n",
                    mpfr_print_rnd_mode ((mpfr_rnd_t) r), i);
            printf ("Expected ");
            mpfr_dump (x);
            printf ("with inex = %d\n", inex1);
            printf ("Got      ");
            mpfr_dump (y);
            printf ("with inex = %d\n", inex2);
            exit (1);
          }
        mpfr_neg (v, v, MPFR_RNDN);
        mpfr_neg (w, w, MPFR_RNDN);
      }
  mpfr_clears (u, v, w, x, y, (mpfr_ptr) 0);
}

static void
cancel (void)
{
  mpfr_t x[2 * MPFR_NCANCEL];
  mpfr_ptr px[2 * MPFR_NCANCEL];
  int i, j, n;

  for (i = 0; i < 8; i++)
    {
      for (n = 0; n < numberof (x); n++)
        {
          mpfr_prec_t p;
          mpfr_rnd_t rnd;

          px[n] = x[n];
          p = MPFR_PREC_MIN + (randlimb () % 256);
          mpfr_init2 (x[n], p);
          if (n < MPFR_NCANCEL)
            {
              mpfr_exp_t e;

              e = (i & 1) ? 0 : mpfr_get_emin ();
              tests_default_random (x[n], 256, e,
                                    ((i & 2) ? e + 2000 : mpfr_get_emax ()));
            }
          else
            {
              /* random permutation with n random transpositions */
              for (j = 0; j < n; j++)
                {
                  int k1, k2;

                  k1 = randlimb () % (n-1);
                  k2 = randlimb () % (n-1);
                  mpfr_swap (x[k1], x[k2]);
                }
              rnd = RND_RAND ();
#if DEBUG
              printf ("mpfr_sum cancellation test\n");
              for (j = 0; j < n; j++)
                {
                  printf ("  x%d[%3ld] = ", j, mpfr_get_prec(x[j]));
                  mpfr_out_str (stdout, 16, 0, x[j], MPFR_RNDN);
                  printf ("\n");
                }
              printf ("  rnd = %s, output prec = %ld\n",
                      mpfr_print_rnd_mode (rnd), mpfr_get_prec (x[n]));
#endif
              mpfr_sum (x[n], px, n, rnd);
              mpfr_neg (x[n], x[n], MPFR_RNDN);
            }
        }

      while (--n >= 0)
        mpfr_clear (x[n]);
    }
}

int
main (void)
{
  tests_start_mpfr ();

  check_special ();
  check_more_special ();
  bug20131027 ();
  generic_tests ();
  check_extreme ();
  cancel ();

  tests_end_mpfr ();
  return 0;
}
