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

  for (m = 0; m < 4000; m++)
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
          if (m % 8 != 0 && (m % 8 == 1 || (randlimb () & 1)))
            mpfr_neg (t[i], t[i], MPFR_RNDN);
          if (non_uniform && MPFR_NOTZERO (t[i]))
            mpfr_set_exp (t[i], randlimb () % 1000);
          /* putchar ("-0+"[SIGN (mpfr_sgn (t[i])) + 1]); */
        }
      /* putchar ('\n'); */
      get_exact_sum (exact_sum, t, n);
      RND_LOOP (rnd_mode)
        {
          int inex1, inex2;

          inex1 = mpfr_set (sum1, exact_sum, (mpfr_rnd_t) rnd_mode);
          inex2 = mpfr_sum (sum2, p, n, (mpfr_rnd_t) rnd_mode);
          if (!(mpfr_equal_p (sum1, sum2) && SAME_SIGN (inex1, inex2)))
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

/* i * 2^46 + j * 2^45 + k * 2^44 + f * 2^(-2),
   with -1 <= i, j, k <= 1, i != 0, -3 <= f <= 3,
   ulp(exact sum) = 2^0 and ulp(exact sum) = 2^44. */
static void
check1 (void)
{
  mpfr_t sum1, sum2, s1, s2, s3, t[4];
  mpfr_ptr p[4];
  int i, j, k, f, prec, r, inex1, inex2;

  mpfr_init2 (sum1, 47);
  mpfr_init2 (sum2, 47);
  mpfr_init2 (s1, 3);
  mpfr_init2 (s2, 3);
  mpfr_init2 (s3, 49);
  for (i = 0; i < 4; i++)
    {
      mpfr_init2 (t[i], 2);
      p[i] = t[i];
    }

  for (i = -1; i <= 1; i += 2)
    {
      mpfr_set_si_2exp (t[0], i, 46, MPFR_RNDN);
      for (j = -1; j <= 1; j++)
        {
          mpfr_set_si_2exp (t[1], j, 45, MPFR_RNDN);
          inex1 = mpfr_add (s1, t[0], t[1], MPFR_RNDN);
          MPFR_ASSERTN (inex1 == 0);
          for (k = -1; k <= 1; k++)
            {
              mpfr_set_si_2exp (t[2], k, 44, MPFR_RNDN);
              inex1 = mpfr_add (s2, s1, t[2], MPFR_RNDN);
              MPFR_ASSERTN (inex1 == 0);
              for (f = -3; f <= 3; f++)
                {
                  mpfr_set_si_2exp (t[3], f, -2, MPFR_RNDN);
                  inex1 = mpfr_add (s3, s2, t[3], MPFR_RNDN);
                  MPFR_ASSERTN (inex1 == 0);
                  for (prec = mpfr_get_exp (s3);
                       prec >= MPFR_PREC_MIN;
                       prec -= 44)
                    {
                      mpfr_set_prec (sum1, prec);
                      mpfr_set_prec (sum2, prec);
                      RND_LOOP (r)
                        {
                          inex1 = mpfr_set (sum1, s3, (mpfr_rnd_t) r);
                          inex2 = mpfr_sum (sum2, p, 4, (mpfr_rnd_t) r);
                          MPFR_ASSERTN (mpfr_check (sum1));
                          MPFR_ASSERTN (mpfr_check (sum2));
                          if (!(mpfr_equal_p (sum1, sum2) &&
                                SAME_SIGN (inex1, inex2)))
                            {
                              printf ("Error in check1 on %s, prec = %d, "
                                      "i = %d, j = %d, k = %d, f = %d\n",
                                      mpfr_print_rnd_mode ((mpfr_rnd_t) r),
                                      prec, i, j, k, f);
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
                }
            }
        }
    }

  for (i = 0; i < 4; i++)
    mpfr_clear (t[i]);
  mpfr_clears (sum1, sum2, s1, s2, s3, (mpfr_ptr) 0);
}

/* With N = 2 * GMP_NUMB_BITS:
   i * 2^N + j + k * 2^(-1) + f1 * 2^(-N) + f2 * 2^(-N),
   with i = -1 or 1, j = 0 or i, -1 <= k <= 1, -1 <= f1 <= 1, -1 <= f2 <= 1
   ulp(exact sum) = 2^0. */
static void
check2 (void)
{
  mpfr_t sum1, sum2, s1, s2, s3, s4, t[5];
  mpfr_ptr p[5];
  int i, j, k, f1, f2, prec, r, inex1, inex2;

#define N (2 * GMP_NUMB_BITS)

  mpfr_init2 (sum1, N+1);
  mpfr_init2 (sum2, N+1);
  mpfr_init2 (s1, N+1);
  mpfr_init2 (s2, N+2);
  mpfr_init2 (s3, 2*N+1);
  mpfr_init2 (s4, 2*N+1);
  for (i = 0; i < 5; i++)
    {
      mpfr_init2 (t[i], 2);
      p[i] = t[i];
    }

  for (i = -1; i <= 1; i += 2)
    {
      mpfr_set_si_2exp (t[0], i, N, MPFR_RNDN);
      for (j = 0; j != 2*i; j += i)
        {
          mpfr_set_si (t[1], j, MPFR_RNDN);
          inex1 = mpfr_add (s1, t[0], t[1], MPFR_RNDN);
          MPFR_ASSERTN (inex1 == 0);
          for (k = -1; k <= 1; k++)
            {
              mpfr_set_si_2exp (t[2], k, -1, MPFR_RNDN);
              inex1 = mpfr_add (s2, s1, t[2], MPFR_RNDN);
              MPFR_ASSERTN (inex1 == 0);
              for (f1 = -1; f1 <= 1; f1++)
                {
                  mpfr_set_si_2exp (t[3], f1, -N, MPFR_RNDN);
                  inex1 = mpfr_add (s3, s2, t[3], MPFR_RNDN);
                  MPFR_ASSERTN (inex1 == 0);
                  for (f2 = -1; f2 <= 1; f2++)
                    {
                      mpfr_set_si_2exp (t[4], f2, -N, MPFR_RNDN);
                      inex1 = mpfr_add (s4, s3, t[4], MPFR_RNDN);
                      MPFR_ASSERTN (inex1 == 0);
                      prec = mpfr_get_exp (s4);
                      mpfr_set_prec (sum1, prec);
                      mpfr_set_prec (sum2, prec);
                      RND_LOOP (r)
                        {
                          inex1 = mpfr_set (sum1, s4, (mpfr_rnd_t) r);
                          inex2 = mpfr_sum (sum2, p, 5, (mpfr_rnd_t) r);
                          MPFR_ASSERTN (mpfr_check (sum1));
                          MPFR_ASSERTN (mpfr_check (sum2));
                          if (!(mpfr_equal_p (sum1, sum2) &&
                                SAME_SIGN (inex1, inex2)))
                            {
                              printf ("Error in check2 on %s, prec = %d, "
                                      "i = %d, j = %d, k = %d, f1 = %d, "
                                      "f2 = %d\n",
                                      mpfr_print_rnd_mode ((mpfr_rnd_t) r),
                                      prec, i, j, k, f1, f2);
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
                }
            }
        }
    }

  for (i = 0; i < 5; i++)
    mpfr_clear (t[i]);
  mpfr_clears (sum1, sum2, s1, s2, s3, s4, (mpfr_ptr) 0);
}

/* t[i] = (2^17 - 1) * 2^(17*(i-8)) for 0 <= i <= 16.
 * t[17] = 2^(17*9+1) * j for -4 <= j <= 4.
 * t[18] = 2^(-1) * k for -1 <= k <= 1.
 * t[19] = 2^(-17*8) * m for -3 <= m <= 3.
 * prec = 17*9+4
 */
static void
check3 (void)
{
  mpfr_t sum1, sum2, s1, s2, s3, s4, t[20];
  mpfr_ptr p[20];
  int i, s, j, k, m, r, inex1, inex2;
  int prec = 17*9+4;

  mpfr_init2 (sum1, prec);
  mpfr_init2 (sum2, prec);
  mpfr_init2 (s1, 17*17);
  mpfr_init2 (s2, 17*17+4);
  mpfr_init2 (s3, 17*17+4);
  mpfr_init2 (s4, 17*17+5);
  mpfr_set_ui (s1, 0, MPFR_RNDN);
  for (i = 0; i < 20; i++)
    {
      mpfr_init2 (t[i], 20);
      p[i] = t[i];
      if (i < 17)
        {
          mpfr_set_ui_2exp (t[i], 0x1ffff, 17*(i-8), MPFR_RNDN);
          inex1 = mpfr_add (s1, s1, t[i], MPFR_RNDN);
          MPFR_ASSERTN (inex1 == 0);
        }
    }

  for (s = 1; s >= -1; s -= 2)
    {
      for (j = -4; j <= 4; j++)
        {
          mpfr_set_si_2exp (t[17], j, 17*9+1, MPFR_RNDN);
          inex1 = mpfr_add (s2, s1, t[17], MPFR_RNDN);
          MPFR_ASSERTN (inex1 == 0);
          for (k = -1; k <= 1; k++)
            {
              mpfr_set_si_2exp (t[18], k, -1, MPFR_RNDN);
              inex1 = mpfr_add (s3, s2, t[18], MPFR_RNDN);
              MPFR_ASSERTN (inex1 == 0);
              for (m = -3; m <= 3; m++)
                {
                  mpfr_set_si_2exp (t[19], m, -17*8, MPFR_RNDN);
                  inex1 = mpfr_add (s4, s3, t[19], MPFR_RNDN);
                  MPFR_ASSERTN (inex1 == 0);
                  RND_LOOP (r)
                    {
                      inex1 = mpfr_set (sum1, s4, (mpfr_rnd_t) r);
                      inex2 = mpfr_sum (sum2, p, 20, (mpfr_rnd_t) r);
                      MPFR_ASSERTN (mpfr_check (sum1));
                      MPFR_ASSERTN (mpfr_check (sum2));
                      if (!(mpfr_equal_p (sum1, sum2) &&
                            SAME_SIGN (inex1, inex2)))
                        {
                          printf ("Error in check3 on %s, "
                                  "s = %d, j = %d, k = %d, m = %d\n",
                                  mpfr_print_rnd_mode ((mpfr_rnd_t) r),
                                  s, j, k, m);
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
            }
        }
      for (i = 0; i < 17; i++)
        mpfr_neg (t[i], t[i], MPFR_RNDN);
      mpfr_neg (s1, s1, MPFR_RNDN);
    }

  for (i = 0; i < 20; i++)
    mpfr_clear (t[i]);
  mpfr_clears (sum1, sum2, s1, s2, s3, s4, (mpfr_ptr) 0);
}

/* Test of s * (q * 2^(n-1) - 2^k) + h + i * 2^(-2) + j * 2^(-2)
 * with h = -1 or 1, -1 <= i odd <= j <= 3, 2 <= q <= 3, s = -1 or 1,
 * prec n-k.
 * On a 64-bit machine:
 * MPFR_RNDN, tmd=2, rbit=0, sst=0, negative is checked with the inputs
 *   -3*2^58, 2^5, -1, 2^(-2), 3*2^(-2)
 * MPFR_RNDN, tmd=2, rbit=0, sst=1, negative is checked with the inputs
 *   -3*2^58, 2^5, -1, 3*2^(-2), 3*2^(-2)
 */
static void
check4 (void)
{
  mpfr_t sum1, sum2, s1, s2, s3, s4, t[5];
  mpfr_ptr p[5];
  int h, i, j, k, n, q, r, s, prec, inex1, inex2;

  mpfr_inits2 (257, sum1, sum2, s1, s2, s3, s4, (mpfr_ptr) 0);
  for (i = 0; i < 5; i++)
    {
      mpfr_init2 (t[i], 2);
      p[i] = t[i];
    }

  /* No GNU style for the many nested loops... */
  for (k = 1; k <= 64; k++) {
    mpfr_set_si_2exp (t[0], -1, k, MPFR_RNDN);
    for (n = k + 2; n <= k + 65; n++) {
      prec = n - k;
      mpfr_set_prec (sum1, prec);
      mpfr_set_prec (sum2, prec);
      for (q = 2; q <= 3; q++) {
        mpfr_set_si_2exp (t[1], q, n - 1, MPFR_RNDN);
        inex1 = mpfr_add (s1, t[0], t[1], MPFR_RNDN);
        MPFR_ASSERTN (inex1 == 0);
        for (s = -1; s <= 1; s += 2) {
          mpfr_neg (t[0], t[0], MPFR_RNDN);
          mpfr_neg (t[1], t[1], MPFR_RNDN);
          mpfr_neg (s1, s1, MPFR_RNDN);
          for (h = -1; h <= 1; h += 2) {
            mpfr_set_si (t[2], h, MPFR_RNDN);
            inex1 = mpfr_add (s2, s1, t[2], MPFR_RNDN);
            MPFR_ASSERTN (inex1 == 0);
            for (i = -1; i <= 3; i += 2) {
              mpfr_set_si_2exp (t[3], i, -2, MPFR_RNDN);
              inex1 = mpfr_add (s3, s2, t[3], MPFR_RNDN);
              MPFR_ASSERTN (inex1 == 0);
              for (j = i; j <= 3; j++) {
                mpfr_set_si_2exp (t[4], j, -2, MPFR_RNDN);
                inex1 = mpfr_add (s4, s3, t[4], MPFR_RNDN);
                MPFR_ASSERTN (inex1 == 0);
                RND_LOOP (r) {
                  inex1 = mpfr_set (sum1, s4, (mpfr_rnd_t) r);
                  inex2 = mpfr_sum (sum2, p, 5, (mpfr_rnd_t) r);
                  MPFR_ASSERTN (mpfr_check (sum1));
                  MPFR_ASSERTN (mpfr_check (sum2));
                  if (!(mpfr_equal_p (sum1, sum2) &&
                        SAME_SIGN (inex1, inex2)))
                    {
                      printf ("Error in check4 on %s, "
                              "k = %d, n = %d (prec %d), "
                              "q = %d, s = %d, h = %d, i = %d, j = %d\n",
                              mpfr_print_rnd_mode ((mpfr_rnd_t) r),
                              k, n, prec, q, s, h, i, j);
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
            }
          }
        }
      }
    }
  }

  for (i = 0; i < 5; i++)
    mpfr_clear (t[i]);
  mpfr_clears (sum1, sum2, s1, s2, s3, s4, (mpfr_ptr) 0);
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

/* Occurs in branches/new-sum/src/sum.c@9344 on a 64-bit machine. */
static void
bug20150327 (void)
{
  mpfr_t sum1, sum2, t[3];
  mpfr_ptr p[3];
  char *s[3] = { "0.10000111110101000010101011100001", "1E-100", "0.1E95" };
  int i, r;

  mpfr_inits2 (58, sum1, sum2, (mpfr_ptr) 0);

  for (i = 0; i < 3; i++)
    {
      mpfr_init2 (t[i], 64);
      mpfr_set_str (t[i], s[i], 2, MPFR_RNDN);
      p[i] = t[i];
    }

  RND_LOOP(r)
    {
      int inex1, inex2;

      mpfr_set (sum1, t[2], MPFR_RNDN);
      inex1 = -1;
      if (MPFR_IS_LIKE_RNDU ((mpfr_rnd_t) r, 1))
        {
          mpfr_nextabove (sum1);
          inex1 = 1;
        }

      inex2 = mpfr_sum (sum2, p, 3, (mpfr_rnd_t) r);

      if (!(mpfr_equal_p (sum1, sum2) && SAME_SIGN (inex1, inex2)))
        {
          printf ("mpfr_sum incorrect in bug20150327 for %s:\n",
                  mpfr_print_rnd_mode ((mpfr_rnd_t) r));
          printf ("Expected ");
          mpfr_dump (sum1);
          printf ("with inex = %d\n", inex1);
          printf ("Got      ");
          mpfr_dump (sum2);
          printf ("with inex = %d\n", inex2);
          exit (1);
        }
    }

  for (i = 0; i < 3; i++)
    mpfr_clear (t[i]);
  mpfr_clears (sum1, sum2, (mpfr_ptr) 0);
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
        if (!(mpfr_equal_p (x, y) && SAME_SIGN (inex1, inex2)))
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

/* Generic random tests with cancellations */
static void
cancel (void)
{
  mpfr_t x[2 * MPFR_NCANCEL], bound;
  mpfr_ptr px[2 * MPFR_NCANCEL];
  int i, j, n;

  mpfr_init2 (bound, 2);

  for (i = 0; i < 8; i++)
    {
      mpfr_set_inf (bound, 1);
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
              if (mpfr_zero_p (x[n]))
                {
                  n++;
                  break;
                }
              mpfr_neg (x[n], x[n], MPFR_RNDN);
              if (mpfr_greater_p (x[n], bound))
                {
                  printf ("Error in cancel on i = %d, n = %d\n", i, n);
                  exit (1);
                }
              mpfr_set_ui_2exp (bound, 1,
                                mpfr_get_exp (x[n]) - p - (rnd == MPFR_RNDN),
                                MPFR_RNDN);
              /* The next sum will be <= bound in absolute value
                 (the equality can be obtained in all rounding modes
                 since the sum will be rounded). */
            }
        }

      while (--n >= 0)
        mpfr_clear (x[n]);
    }

  mpfr_clear (bound);
}

#ifndef NOVFL
# define NOVFL 30
#endif

static void
check_overflow (void)
{
  mpfr_t sum1, sum2, x, y;
  mpfr_ptr t[2 * NOVFL];
  mpfr_exp_t emin, emax;
  int i, r;

  emin = mpfr_get_emin ();
  emax = mpfr_get_emax ();
  set_emin (MPFR_EMIN_MIN);
  set_emax (MPFR_EMAX_MAX);

  mpfr_inits2 (32, sum1, sum2, x, y, (mpfr_ptr) 0);
  mpfr_setmax (x, mpfr_get_emax ());
  mpfr_neg (y, x, MPFR_RNDN);

  for (i = 0; i < 2 * NOVFL; i++)
    t[i] = i < NOVFL ? x : y;

  /* Two kinds of test:
   *   i = 1: overflow.
   *   i = 2: intermediate overflow (exact sum is 0).
   */
  for (i = 1; i <= 2; i++)
    RND_LOOP(r)
      {
        int inex1, inex2;

        inex1 = mpfr_add (sum1, x, i == 1 ? x : y, (mpfr_rnd_t) r);
        inex2 = mpfr_sum (sum2, t, i * NOVFL, (mpfr_rnd_t) r);
        MPFR_ASSERTN (mpfr_check (sum1));
        MPFR_ASSERTN (mpfr_check (sum2));
        if (!(mpfr_equal_p (sum1, sum2) && SAME_SIGN (inex1, inex2)))
          {
            printf ("Error in check_overflow on %s, i = %d\n",
                    mpfr_print_rnd_mode ((mpfr_rnd_t) r), i);
            printf ("Expected ");
            mpfr_dump (sum1);
            printf ("with inex = %d\n", inex1);
            printf ("Got      ");
            mpfr_dump (sum2);
            printf ("with inex = %d\n", inex2);
            exit (1);
          }
      }

  mpfr_clears (sum1, sum2, x, y, (mpfr_ptr) 0);

  set_emin (emin);
  set_emax (emax);
}

#ifndef NUNFL
# define NUNFL 9
#endif

/* t[0] = 2^(-k) - sum(t[i],i=1..n)
 */
static void
check_underflow (void)
{
  mpfr_t sum1, sum2, t[NUNFL];
  mpfr_ptr p[NUNFL];
  mpfr_prec_t precmax = 444;
  mpfr_exp_t emin, emax;
  unsigned int ex_flags, flags;
  int c, i;

  emin = mpfr_get_emin ();
  emax = mpfr_get_emax ();
  set_emin (MPFR_EMIN_MIN);
  set_emax (MPFR_EMAX_MAX);

  ex_flags = MPFR_FLAGS_UNDERFLOW | MPFR_FLAGS_INEXACT;

  mpfr_init2 (sum1, MPFR_PREC_MIN);
  mpfr_init2 (sum2, precmax);

  for (i = 0; i < NUNFL; i++)
    {
      mpfr_init2 (t[i], precmax);
      p[i] = t[i];
    }

  for (c = 0; c < 8; c++)
    {
      mpfr_prec_t fprec;
      int n, neg, r;

      fprec = MPFR_PREC_MIN + (randlimb () % (precmax - MPFR_PREC_MIN + 1));
      n = 3 + (randlimb () % (NUNFL - 2));
      MPFR_ASSERTN (n <= NUNFL);

      mpfr_set_prec (sum2, (randlimb () & 1) ? MPFR_PREC_MIN : precmax);
      mpfr_set_prec (t[0], fprec + 64);
      mpfr_set_zero (t[0], 1);

      for (i = 1; i < n; i++)
        {
          int inex;

          mpfr_set_prec (t[i], MPFR_PREC_MIN +
                         (randlimb () % (fprec - MPFR_PREC_MIN + 1)));
          do
            mpfr_urandomb (t[i], RANDS);
          while (MPFR_IS_ZERO (t[i]));
          mpfr_set_exp (t[i], MPFR_EMIN_MIN);
          inex = mpfr_sub (t[0], t[0], t[i], MPFR_RNDN);
          MPFR_ASSERTN (inex == 0);
        }

      neg = randlimb () & 1;
      if (neg)
        mpfr_nextbelow (t[0]);
      else
        mpfr_nextabove (t[0]);

      RND_LOOP(r)
        {
          int inex1, inex2;

          mpfr_set_zero (sum1, 1);
          if (neg)
            mpfr_nextbelow (sum1);
          else
            mpfr_nextabove (sum1);
          inex1 = mpfr_div_2ui (sum1, sum1, 2, (mpfr_rnd_t) r);

          mpfr_clear_flags ();
          inex2 = mpfr_sum (sum2, p, n, (mpfr_rnd_t) r);
          flags = __gmpfr_flags;

          MPFR_ASSERTN (mpfr_check (sum1));
          MPFR_ASSERTN (mpfr_check (sum2));

          if (flags != ex_flags)
            {
              printf ("Bad flags in check_underflow on %s, c = %d\n",
                      mpfr_print_rnd_mode ((mpfr_rnd_t) r), c);
              printf ("Expected flags:");
              flags_out (ex_flags);
              printf ("Got flags:     ");
              flags_out (flags);
              printf ("sum = ");
              mpfr_dump (sum2);
              exit (1);
            }

          if (!(mpfr_equal_p (sum1, sum2) && SAME_SIGN (inex1, inex2)))
            {
              printf ("Error in check_underflow on %s, c = %d\n",
                      mpfr_print_rnd_mode ((mpfr_rnd_t) r), c);
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

  for (i = 0; i < NUNFL; i++)
    mpfr_clear (t[i]);
  mpfr_clears (sum1, sum2, (mpfr_ptr) 0);

  set_emin (emin);
  set_emax (emax);
}

static void
check_coverage (void)
{
#ifdef MPFR_COV_CHECK
  int r, i, j, k, p;
  int err = 0;

  for (r = 0; r < MPFR_RND_MAX; r++)
    for (i = 0; i < 1 + ((mpfr_rnd_t) r == MPFR_RNDN); i++)
      for (j = 0; j < 2; j++)
        for (k = 0; k < 3; k++)
          for (p = 0; p < 2; p++)
            if (!__gmpfr_cov_sum_tmd[r][i][j][k][p])
              {
                printf ("TMD not tested on %s, tmd=%d, rbit=%d, sst=%d, %s\n",
                        mpfr_print_rnd_mode ((mpfr_rnd_t) r), i+1, j, k-1,
                        p ? "positive" : "negative");
                err = 1;
              }

  if (err)
    exit (1);
#endif
}

int
main (void)
{
  tests_start_mpfr ();

  check_special ();
  check_more_special ();
  check1 ();
  check2 ();
  check3 ();
  check4 ();
  bug20131027 ();
  bug20150327 ();
  generic_tests ();
  check_extreme ();
  cancel ();
  check_overflow ();
  check_underflow ();

  check_coverage ();
  tests_end_mpfr ();
  return 0;
}
