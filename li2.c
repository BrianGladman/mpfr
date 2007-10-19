/* mpfr_li2 -- Dilogarithm.

Copyright 2007 Free Software Foundation, Inc.
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

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* assuming B[0]...B[2(n-1)] are computed, computes and stores B[2n]*(2n+1)!

   t/(exp(t)-1) = sum(B[j]*t^j/j!, j=0..infinity)
   thus t = (exp(t)-1) * sum(B[j]*t^j/j!, n=0..infinity).
   Taking the coefficient of degree n+1 > 1, we get:
   0 = sum(1/(n+1-k)!*B[k]/k!, k=0..n)
   which gives:
   B[n] = -sum(binomial(n+1,k)*B[k], k=0..n-1)/(n+1).

   Let C[n] = B[n]*(n+1)!.
   Then C[n] = -sum(binomial(n+1,k)*C[k]*n!/(k+1)!,  k=0..n-1),
   which proves that the C[n] are integers.
*/
static mpz_t *
bernoulli (mpz_t * b, unsigned long n)
{
  if (n == 0)
    {
      b = (mpz_t *) (*__gmp_allocate_func) (sizeof (mpz_t));
      mpz_init_set_ui (b[0], 1);
    }
  else
    {
      mpz_t t;
      unsigned long k;

      b = (mpz_t *) (*__gmp_reallocate_func)
        (b, n * sizeof (mpz_t), (n + 1) * sizeof (mpz_t));
      mpz_init (b[n]);
      /* b[n] = -sum(binomial(2n+1,2k)*C[k]*(2n)!/(2k+1)!,  k=0..n-1) */
      mpz_init_set_ui (t, 2 * n + 1);
      mpz_mul_ui (t, t, 2 * n - 1);
      mpz_mul_ui (t, t, 2 * n);
      mpz_mul_ui (t, t, n);
      mpz_div_ui (t, t, 3);     /* exact: t=binomial(2*n+1,2*k)*(2*n)!/(2*k+1)!
                                   for k=n-1 */
      mpz_mul (b[n], t, b[n - 1]);
      for (k = n - 1; k-- > 0;)
        {
          mpz_mul_ui (t, t, 2 * k + 1);
          mpz_mul_ui (t, t, 2 * k + 2);
          mpz_mul_ui (t, t, 2 * k + 2);
          mpz_mul_ui (t, t, 2 * k + 3);
          mpz_div_ui (t, t, 2 * (n - k) + 1);
          mpz_div_ui (t, t, 2 * (n - k));
          mpz_addmul (b[n], t, b[k]);
        }
      /* take into account C[1] */
      mpz_mul_ui (t, t, 2 * n + 1);
      mpz_div_2exp (t, t, 1);
      mpz_sub (b[n], b[n], t);
      mpz_neg (b[n], b[n]);
      mpz_clear (t);
    }
  return b;
}

/* Compute the alternating series
   s = S(z) = \sum_{k=0}^infty B_{2k} (z))^{2k+1} / (2k+1)! 
   with 0 < z <= log(2) to the precision of s rounded in the direction
   rnd_mode. 
   Return the maximum index of the truncature which is useful 
   for determinating the relative error.
*/
static int
li2_series (mpfr_t sum, mpfr_srcptr z, mpfr_rnd_t rnd_mode)
{
  int inexact;
  unsigned int i, Bm, Bmax;
  mpfr_t s, u, v, w;
  mpfr_prec_t sump, p;
  mp_exp_t se, err;
  mpz_t *B;
  MPFR_ZIV_DECL (loop);

  /* The series converges for |z| < 2 pi, but in mpfr_li2 the argument is 
     reduced so that 0 < z <= log(2). Here is additionnal check that z is
     (nearly) correct */
  MPFR_ASSERTD (MPFR_IS_STRICTPOS (z));
  MPFR_ASSERTD (mpfr_cmp_d (z, 0.6953125) <= 0);

  sump = MPFR_PREC (sum);  /* target precision */
  p = sump + MPFR_INT_CEIL_LOG2 (sump) + 4; /* the working precision */
  mpfr_init2 (s, p);
  mpfr_init2 (u, p);
  mpfr_init2 (v, p);
  mpfr_init2 (w, p);

  B = bernoulli ((mpz_t *) 0, 0);
  Bm = Bmax = 1;

  MPFR_ZIV_INIT (loop, p);
  for (;;)
    {
      mpfr_sqr (u, z, GMP_RNDU);
      mpfr_set (v, z, GMP_RNDU);
      mpfr_set (s, z, GMP_RNDU);
      se = MPFR_GET_EXP (s);
      err = 0;

      for (i = 1;; i++)
	{
	  if (i >= Bmax)
	    B = bernoulli (B, Bmax++);  /* B_2i * (2i+1)!, exact */

	  mpfr_mul (v, u, v, GMP_RNDU);
	  mpfr_div_ui (v, v, 2 * i, GMP_RNDU);
	  mpfr_div_ui (v, v, 2 * i, GMP_RNDU);
	  mpfr_div_ui (v, v, 2 * i + 1, GMP_RNDU);
	  mpfr_div_ui (v, v, 2 * i + 1, GMP_RNDU);  
	  /* here, v_2i = v_{2i-2} / (2i * (2i+1))^2 */

	  mpfr_mul_z (w, v, B[i], GMP_RNDN);  
	  /* here, w_2i = v_2i * B_2i * (2i+1)! with
	     error(w_2i) < 2^(5 * i + 8) ulp(w_2i) (see algorithm.tex) */ 

	  mpfr_add (s, s, w, GMP_RNDN);
	  
	  err = MAX (err + se, 5 * i + 8 + MPFR_GET_EXP (w)) 
	    - MPFR_GET_EXP (s) + 1;
	  se = MPFR_GET_EXP (s);
	  if (MPFR_GET_EXP (w) <= se - (mp_exp_t)p)
	    break;
	}
      
      /* the previous value of err is the rounding error, 
	 the truncation error is less than EXP(z) - 4 * i - 4 
	 (see algorithm.tex)*/
      err = MAX (err, MPFR_GET_EXP (z) - 4 * i - 4) + 1;
      if (MPFR_CAN_ROUND (s, p - err, sump, rnd_mode))
	break;

      MPFR_ZIV_NEXT (loop, p);
      mpfr_set_prec (s, p);
      mpfr_set_prec (u, p);
      mpfr_set_prec (v, p);
      mpfr_set_prec (w, p);
    }
  MPFR_ZIV_FREE (loop);
  inexact = mpfr_set (sum, s, rnd_mode);

  Bm = Bmax;
  while (Bm--)
    mpz_clear (B[Bm]);
  (*__gmp_free_func) (B, Bmax * sizeof (mpz_t));
  mpfr_clears (s, u, v, w, (void *) 0);

  /* Let K be the returned value.
     As we compute an alternating series, the truncation error has the same
     sign as the next term w_{K+2} which is positive iff K%4 == 0.
     Assume that error(z) <= (1+t) z', where z' is the actual value, then
     error(s) <= 2 * (K+1) * t (see algorithm.tex).
  */
  return 2*i;
}

/* Compute the real part of the dilogarithm defined by
   Li2(x) = -\Int_{t=0}^x log(1-t)/t dt */
int
mpfr_li2 (mpfr_ptr y, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  int inexact;
  mp_exp_t err;
  mpfr_prec_t yp, m;
  MPFR_ZIV_DECL (loop);
  MPFR_SAVE_EXPO_DECL (expo);

  MPFR_LOG_FUNC (("x[%#R]=%R rnd=%d", x, x, rnd_mode), ("y[%#R]=%R", y));

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    {
      if (MPFR_IS_NAN (x))
        {
          MPFR_SET_NAN (y);
          MPFR_RET_NAN;
        }
      else if (MPFR_IS_INF (x))
        {
          MPFR_SET_NEG (y);
          MPFR_SET_INF (y);
          MPFR_RET (0);
        }
      else                      /* x is zero */
        {
          MPFR_ASSERTD (MPFR_IS_ZERO (x));
          MPFR_SET_SAME_SIGN (y, x);
          MPFR_SET_ZERO (y);
          MPFR_RET (0);
        }
    }

  /* Li2(x) = x + x^2/4 + x^3/9 + ... */ 
  if (MPFR_IS_POS (x))
    MPFR_FAST_COMPUTE_IF_SMALL_INPUT (y, x,  -MPFR_GET_EXP (x), 0, 1,
				      rnd_mode, {});
  else
    MPFR_FAST_COMPUTE_IF_SMALL_INPUT (y, x,  -MPFR_GET_EXP (x), 2, 0,
				      rnd_mode, {});

  MPFR_SAVE_EXPO_MARK (expo);
  yp = MPFR_PREC (y);
  m = yp + MPFR_INT_CEIL_LOG2 (yp) + 13;

  if (MPFR_LIKELY ((mpfr_cmp_ui (x, 0) > 0) && (mpfr_cmp_d (x, 0.5) <= 0)))
    /* 0 < x <= 1/2: Li2(x) = S(-log(1-x))-log^2(1-x)/4 */
    {
      mpfr_t s, u;
      mp_exp_t d, expo_l;
      int k;

      mpfr_init2 (u, m);
      mpfr_init2 (s, m);

      MPFR_ZIV_INIT (loop, m);
      for (;;)
        {
          mpfr_ui_sub (u, 1, x, GMP_RNDN);
          mpfr_log (u, u, GMP_RNDU);
          mpfr_neg (u, u, GMP_RNDN);    /* u = -log(1-x) */
          expo_l = MPFR_GET_EXP (u);
          k = li2_series (s, u, GMP_RNDU);
          d = 2 - expo_l + MPFR_INT_CEIL_LOG2 (k + 1) + MPFR_GET_EXP (s);

          mpfr_sqr (u, u, GMP_RNDU);
          mpfr_div_2ui (u, u, 2, GMP_RNDU);     /* u = log^2(1-x) / 4 */
          mpfr_sub (s, s, u, GMP_RNDN);

          /* error(s) <= (0.5 + 2^(d-EXP(s)) + 2^(6+expo_l-EXP(s))) ulp(s) */
          err = (mp_exp_t) MAX (d, 6 + expo_l) - MPFR_GET_EXP (s) + 1;
          if (MPFR_CAN_ROUND (s, m - err, yp, rnd_mode))
            break;

          MPFR_ZIV_NEXT (loop, m);
          mpfr_set_prec (u, m);
          mpfr_set_prec (s, m);
        }
      MPFR_ZIV_FREE (loop);
      inexact = mpfr_set (y, s, rnd_mode);

      mpfr_clear (u);
      mpfr_clear (s);
      MPFR_SAVE_EXPO_FREE (expo);
      return mpfr_check_range (y, inexact, rnd_mode);
    }
  else if (!mpfr_cmp_ui (x, 1))
    /* Li2(1)= pi^2 / 6 */
    {
      mpfr_t u;
      mpfr_init2 (u, m);

      MPFR_ZIV_INIT (loop, m);
      for (;;)
        {
          mpfr_const_pi (u, GMP_RNDU);
          mpfr_sqr (u, u, GMP_RNDN);
          mpfr_div_ui (u, u, 6, GMP_RNDN);

          err = m - 4;          /* error(u) <= 19/2 ulp(u) */
          if (MPFR_CAN_ROUND (u, err, yp, rnd_mode))
            break;

          MPFR_ZIV_NEXT (loop, m);
          mpfr_set_prec (u, m);
        }
      MPFR_ZIV_FREE (loop);
      inexact = mpfr_set (y, u, rnd_mode);

      mpfr_clear (u);
      MPFR_SAVE_EXPO_FREE (expo);
      return mpfr_check_range (y, inexact, rnd_mode);
    }
  else if (mpfr_cmp_ui (x, 2) >= 0)
    /* x >= 2: Li2(x) = -S(-log(1-1/x))-log^2(x)/2+log^2(1-1/x)/4-pi^2/3 */
    {
      int k;
      mp_exp_t d, expo_l;
      mpfr_t s, u, xx;
      mpfr_init2 (u, m);
      mpfr_init2 (s, m);
      mpfr_init2 (xx, m);

      MPFR_ZIV_INIT (loop, m);
      for (;;)
        {
          mpfr_ui_div (xx, 1, x, GMP_RNDN);
          mpfr_neg (xx, xx, GMP_RNDN);
          mpfr_log1p (u, xx, GMP_RNDN);
          mpfr_neg (u, u, GMP_RNDU);    /* u = -log(1-1/x) */
          expo_l = MPFR_GET_EXP (u);
          k = li2_series (s, u, GMP_RNDN);
          d = 2 - expo_l + MPFR_INT_CEIL_LOG2 (k + 1) + MPFR_GET_EXP (s);

          mpfr_neg (s, s, GMP_RNDN);
          mpfr_sqr (u, u, GMP_RNDN);
          mpfr_div_2ui (u, u, 2, GMP_RNDN);     /* u= log^2(1-1/x)/4 */
          mpfr_add (s, s, u, GMP_RNDN);
          d = MAX (d, 4 + expo_l) - MPFR_GET_EXP (s) + 1;

          mpfr_log (u, x, GMP_RNDU);
          mpfr_sqr (u, u, GMP_RNDN);
          mpfr_div_2ui (u, u, 1, GMP_RNDN);     /* u = log^2(x)/2 */
          mpfr_sub (s, s, u, GMP_RNDN);
          d = MAX (d, 3 + MPFR_GET_EXP (u)) - MPFR_GET_EXP (s) + 1;

          mpfr_const_pi (u, GMP_RNDU);
          mpfr_sqr (u, u, GMP_RNDN);
          mpfr_div_ui (u, u, 3, GMP_RNDN);      /* u = pi^2/3 */
          mpfr_add (s, s, u, GMP_RNDN);

          err =
            (mp_exp_t) MAX (d, 3 + MPFR_GET_EXP (u)) - MPFR_GET_EXP (s) + 1;
          if (MPFR_CAN_ROUND (s, m - err, yp, rnd_mode))
            break;

          MPFR_ZIV_NEXT (loop, m);
          mpfr_set_prec (u, m);
          mpfr_set_prec (s, m);
          mpfr_set_prec (xx, m);
        }
      MPFR_ZIV_FREE (loop);
      inexact = mpfr_set (y, s, rnd_mode);

      mpfr_clears (s, u, xx, (void *) 0);
      MPFR_SAVE_EXPO_FREE (expo);
      return mpfr_check_range (y, inexact, rnd_mode);
    }
  else if (mpfr_cmp_ui (x, 1) > 0)
    /* 2 > x > 1: Li2(x) = S(log(x))+log^2(x)/4-log(x)log(x-1)+pi^2/6 */
    {
      int k;
      long d;
      mpfr_t s, u, v, xx;
      mpfr_init2 (s, m);
      mpfr_init2 (u, m);
      mpfr_init2 (v, m);
      mpfr_init2 (xx, m);

      MPFR_ZIV_INIT (loop, m);
      for (;;)
        {
          mpfr_log (v, x, GMP_RNDU);
          k = li2_series (s, v, rnd_mode);
          d = 4 * (k + 1);      /* s > 0, error(s)<= d ulp(s) */

          mpfr_sqr (u, v, GMP_RNDN);
          mpfr_div_2ui (u, u, 2, GMP_RNDN);     /* u = log^2(x)/4 */
          mpfr_add (s, s, u, GMP_RNDN);
          d += 9;               /* s*u > 0 */

          mpfr_sub_ui (xx, x, 1, GMP_RNDN);
          mpfr_log (u, xx, GMP_RNDU);
          mpfr_mul (u, v, u, GMP_RNDN); /* u = log(x) * log(x-1) */
          mpfr_sub (s, s, u, GMP_RNDN);
          d += 2;               /* 0.5 > (-u), s*(-u) > 0 */

          mpfr_const_pi (u, GMP_RNDN);
          mpfr_sqr (u, u, GMP_RNDN);
          mpfr_div_ui (u, u, 6, GMP_RNDN);      /* u = pi^2/6 */
          mpfr_add (s, s, u, GMP_RNDN);
          err = (mp_exp_t) MPFR_INT_CEIL_LOG2 (d + 9);
          if (MPFR_CAN_ROUND (s, m - err, yp, rnd_mode))
            break;

          MPFR_ZIV_NEXT (loop, m);
          mpfr_set_prec (s, m);
          mpfr_set_prec (u, m);
          mpfr_set_prec (v, m);
          mpfr_set_prec (xx, m);
        }
      MPFR_ZIV_FREE (loop);
      inexact = mpfr_set (y, s, rnd_mode);

      mpfr_clears (s, u, v, xx, (void *) 0);
      MPFR_SAVE_EXPO_FREE (expo);
      return mpfr_check_range (y, inexact, rnd_mode);
    }
  else if (mpfr_cmp_d (x, .5) > 0)
    /* 1 > x > 1/2: Li2(x) = -S(-log(x))+log^2(x)/4-log(x)log(1-x)+pi^2/6 */
    {
      int k;
      mp_exp_t d, expo_l;
      mpfr_t s, u, v, xx;
      mpfr_init2 (s, m);
      mpfr_init2 (u, m);
      mpfr_init2 (v, m);
      mpfr_init2 (xx, m);


      MPFR_ZIV_INIT (loop, m);
      for (;;)
        {
          mpfr_log (u, x, GMP_RNDN);
          mpfr_neg (u, u, GMP_RNDN);
          k = li2_series (s, u, GMP_RNDN);
          mpfr_neg (s, s, GMP_RNDN);
          expo_l = MPFR_GET_EXP (u);
          d = 2 - expo_l + MPFR_INT_CEIL_LOG2 (k + 1) + MPFR_GET_EXP (s);

          mpfr_ui_sub (xx, 1, x, GMP_RNDN);
          mpfr_log (v, xx, GMP_RNDN);
          mpfr_mul (v, v, u, GMP_RNDN); /* v = - log(x) * log(1-x) */
          mpfr_add (s, s, v, GMP_RNDN);
          d = MAX (d, MPFR_GET_EXP (v)) - MPFR_GET_EXP (s) + 1;

          mpfr_sqr (u, u, GMP_RNDN);
          mpfr_div_2ui (u, u, 2, GMP_RNDN);     /* u = log^2(x)/4 */
          mpfr_add (s, s, u, GMP_RNDN);
          d = MAX (d - MPFR_GET_EXP (s), 2 + expo_l - MPFR_GET_EXP (s)) + 1;

          mpfr_const_pi (u, GMP_RNDN);
          mpfr_sqr (u, u, GMP_RNDN);
          mpfr_div_ui (u, u, 6, GMP_RNDN);      /* u = pi^2/6 */
          mpfr_add (s, s, u, GMP_RNDN);
          err =
            (mp_exp_t) MAX (d, 3 + MPFR_GET_EXP (u)) - MPFR_GET_EXP (s) + 1;
          if (MPFR_CAN_ROUND (s, m - err, yp, rnd_mode))
            break;

          MPFR_ZIV_NEXT (loop, m);
          mpfr_set_prec (s, m);
          mpfr_set_prec (u, m);
          mpfr_set_prec (v, m);
          mpfr_set_prec (xx, m);
        }
      MPFR_ZIV_FREE (loop);
      inexact = mpfr_set (y, s, rnd_mode);

      mpfr_clears (s, u, v, xx, (void *) 0);
      MPFR_SAVE_EXPO_FREE (expo);
      return mpfr_check_range (y, inexact, rnd_mode);
    }
  else if (mpfr_cmp_si (x, -1) >= 0)
    /* 0> x >= -1: Li2(x) = -S(log(1-x))-log^2(1-x)/4 */
    {
      int k;
      mp_exp_t d, expo_l;
      mpfr_t s, u, xx;
      mpfr_init2 (s, m);
      mpfr_init2 (u, m);
      mpfr_init2 (xx, m);

      MPFR_ZIV_INIT (loop, m);
      for (;;)
        {
	  mpfr_neg (xx, x, GMP_RNDN);
	  mpfr_log1p (u, xx, GMP_RNDN);
	  k = li2_series (s, u, GMP_RNDN);
	  mpfr_neg (s, s, GMP_RNDN);
          expo_l = MPFR_GET_EXP (u);
          d = 2 - expo_l + MPFR_INT_CEIL_LOG2 (k + 1) + MPFR_GET_EXP (s);

	  mpfr_sqr (u, u, GMP_RNDN);
	  mpfr_div_2ui (u, u, 2, GMP_RNDN); /* u = log^2(1-x)/4 */
	  mpfr_sub (s, s, u, GMP_RNDN);
          err =
            (mp_exp_t) MAX (d, 2 + MPFR_GET_EXP (u)) - MPFR_GET_EXP (s) + 1;
          if (MPFR_CAN_ROUND (s, m - err, yp, rnd_mode))
            break;

          MPFR_ZIV_NEXT (loop, m);
          mpfr_set_prec (s, m);
          mpfr_set_prec (u, m);
          mpfr_set_prec (xx, m);
        }
      MPFR_ZIV_FREE (loop);
      inexact = mpfr_set (y, s, rnd_mode);

      mpfr_clears (s, u, xx, (void *) 0);
      MPFR_SAVE_EXPO_FREE (expo);
      return mpfr_check_range (y, inexact, rnd_mode);
    }
  else
    /* x < -1: Li2(x) 
       = S(log(1-1/x))-ln^2(-x)/4-log(1-x)log(-x)/2+log^2(1-x)/4-pi^2/6 */
    {
      int k;
      mp_exp_t d, expo_l;
      mpfr_t s, u, v, w, xx;
      mpfr_init2 (s, m);
      mpfr_init2 (u, m);
      mpfr_init2 (v, m);
      mpfr_init2 (w, m);
      mpfr_init2 (xx, m);

      MPFR_ZIV_INIT (loop, m);
      for (;;)
        {
	  mpfr_ui_div (xx, 1, x, GMP_RNDN);
	  mpfr_neg (xx, xx, GMP_RNDN);
	  mpfr_log1p (u, xx, GMP_RNDU);
	  expo_l = MPFR_GET_EXP (u);
	  k = li2_series (s, u, GMP_RNDN);
	  d = 2 - expo_l + MPFR_INT_CEIL_LOG2 (k + 1) + MPFR_GET_EXP (s);

	  mpfr_ui_sub (xx, 1, x, GMP_RNDN);
	  mpfr_log (u, xx, GMP_RNDU);
	  mpfr_neg (xx, x, GMP_RNDN);
	  mpfr_log (v, xx, GMP_RNDU);
	  mpfr_mul (w, v, u, GMP_RNDN);
	  mpfr_div_2ui (w, w, 1, GMP_RNDN); /* w = log(-x) * log(1-x) / 2 */
	  mpfr_sub (s, s, w, GMP_RNDN);
	  d = MAX (d, 3 + MPFR_GET_EXP (w)) - MPFR_GET_EXP (s) + 1;

	  mpfr_sqr (w, v, GMP_RNDN);
	  mpfr_div_2ui (w, w, 2, GMP_RNDN); /* w = log^2(-x) / 4 */
	  mpfr_sub (s, s, w, GMP_RNDN);
	  d = MAX (d, 3) - MPFR_GET_EXP (s) + 1;

	  mpfr_sqr (w, u, GMP_RNDN);
	  mpfr_div_2ui (w, w, 2, GMP_RNDN); /* w = log^2(1-x) / 4 */
	  mpfr_add (s, s, w, GMP_RNDN);
	  d = MAX (d, 4) - MPFR_GET_EXP (s) + 1;

	  mpfr_const_pi (w, GMP_RNDN);
	  mpfr_sqr (w, w, GMP_RNDN);
	  mpfr_div_ui (w, w, 6, GMP_RNDN);  /* w = pi^2 / 6 */
	  mpfr_sub (s, s, w, GMP_RNDN);
	  err = 
	    (mp_exp_t)MAX (d, 3 + MPFR_GET_EXP (u)) - MPFR_GET_EXP (s) + 1;
	  if (MPFR_CAN_ROUND (s, m - err, yp, rnd_mode))
	    break;

	  MPFR_ZIV_NEXT (loop, m);
	  mpfr_set_prec (s, m);
	  mpfr_set_prec (u, m);
	  mpfr_set_prec (v, m);
	  mpfr_set_prec (w, m);
	  mpfr_set_prec (xx, m);
	}
      MPFR_ZIV_FREE (loop);
      inexact = mpfr_set (y, s, rnd_mode);

      mpfr_clears (s, u, v, w, xx, (void *) 0);
      MPFR_SAVE_EXPO_FREE (expo);
      return mpfr_check_range (y, inexact, rnd_mode);
    }

  MPFR_ASSERTN (0);             /* should never reach this point */
}
