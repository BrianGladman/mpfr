/* mpfr_trigamma -- trigamma function of a floating-point number

Copyright 2024 Free Software Foundation, Inc.
Contributed by the AriC and Caramba projects, INRIA.

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
along with the GNU MPFR Library; see the file COPYING.LESSER.
If not, see <https://www.gnu.org/licenses/>. */

#include "mpfr-impl.h"

/* References:
   [1] Algorithm AS 121: Trigamma Function. B. E. Schneider,
       Journal of the Royal Statistical Society. Series C (Applied Statistics),
       Vol. 27, No. 1 (1978), pp. 97-99 (3 pages)
*/

/* compute trigamma(x) for x < 1/2 using the reflection formula:
   trigamma(1-x) + trigamma(x) = pi^2*(1+cot(pi*x)^2)
   thus for x < 1/2:
   (a) evaluate z = trigamma(y), where y = 1-x >= 1/2
   (b) return pi^2*(1+cot(pi*x)^2) - z
*/
static int
mpfr_trigamma_reflection (mpfr_ptr y, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  /* trigamma(n) = +Inf for n integer, n <= 0 */
  if (mpfr_integer_p (x))
    {
      MPFR_SET_INF(y);
      MPFR_SET_POS(y);
      return 0;
    }

  return mpfr_set (y, x, rnd_mode); // not yet implemented
}

/* special case for x=1/2, where trigamma(1/2) = pi^2/2 */
static int
half (mpfr_ptr y, mpfr_rnd_t rnd_mode)
{
  mpfr_t t;
  int inex;
  mpfr_prec_t p = MPFR_PREC(y), extra;
  /* We can check by exhaustive search that the algorithm below returns the
     correct rounding for extra = 7 for p <= 152,
     then for extra = 14 for p <= 9629 bits.
     Then extra = 64 is probably overkill for all possible values of p,
     but it would be more rigorous to have a Ziv loop (FIXME). */
  extra = (p <= 152) ? 7 : (p <= 9629) ? 14 : 64;
  mpfr_init2 (t, MPFR_PREC(y) + extra);
  mpfr_const_pi (t, MPFR_RNDN);
  mpfr_sqr (t, t, MPFR_RNDN);
  inex = mpfr_div_2ui (y, t, 1, rnd_mode);
  mpfr_clear (t);
  return inex;
}

/* Put in s an approximation of trigamma(x).
   Assumes x >= 2.
   Assumes s does not overlap with x.
   Returns an integer e such that the error is bounded by 2^e ulps
   of the result s.
   Use the formula (6.4.11) with n=1 from Abramowitz & Stegun:
   trigamma(x) = 1/x + 1/(2x^2) + sum(B[2j]/x^(2j+1), j=1..infinity)
   where B[2j] are Bernoulli numbers, which we rewrite as:
   trigamma(x) = 1/x * (1 + 1/(2x) + sum(B[2j]/x^(2j), j=1..infinity))
*/
static mpfr_exp_t
mpfr_trigamma_approx (mpfr_ptr s, mpfr_srcptr x)
{
  mpfr_prec_t p = MPFR_PREC (s);
  mpfr_t t, u, invxx;
  mpfr_exp_t e, exps, f, expu;
  unsigned long n;

  MPFR_ASSERTN (MPFR_IS_POS (x) && MPFR_GET_EXP (x) >= 2);

  mpfr_init2 (t, p);
  mpfr_init2 (u, p);
  mpfr_init2 (invxx, p);

  mpfr_set_ui (s, 1, MPFR_RNDN);     /* exact */
  mpfr_ui_div (t, 1, x, MPFR_RNDN);  /* error <= 1/2 ulp */
  mpfr_div_2ui (t, t, 1, MPFR_RNDN); /* exact */
  mpfr_add (s, s, t, MPFR_RNDN);     /* error <= 1/2 ulp */
  /* since x >= 2, we have t <= 1/2 thus the 1/2 ulp error on t = 1/x
     translates to 1/4 ulp(1), and to 1/8 ulp(1) after t = t/2,
     which stays <= 1/8 ulp(1) after the addition of s and t, thus the error
     so far is bounded by 1/8 + 1/2 < 1 ulp(s) */
  e = 1; /* initial error in ulp(s) */
  /* Note: the values 'theta' below may represent different values,
     all with |theta| <= 2^-p, following Higham's method */
  mpfr_sqr (invxx, x, MPFR_RNDZ);     /* invxx = x^2 * (1 + theta) */
  mpfr_ui_div (invxx, 1, invxx, MPFR_RNDU); /* invxx = 1/x^2 * (1 + theta)^2 */

  /* in the following we note err=xxx when the ratio between the approximation
     and the exact result can be written (1 + theta)^xxx for |theta| <= 2^-p */
  mpfr_set_ui (t, 1, MPFR_RNDN); /* err = 0 */
  for (n = 1;; n++)
    {
      /* The main term is Bernoulli[2n]/x^(2n) = b[n]/(2n+1)!/x^(2n)
         = b[n]*t[n] where t[n]/t[n-1] = 1/(2n+1)/x^2. */
      mpfr_mul (t, t, invxx, MPFR_RNDU);        /* err = err + 3 */
      mpfr_div_ui (t, t, 2 * n + 1, MPFR_RNDU); /* err = err + 1 */
      /* we thus have err = 4n here */
      mpfr_mul_z (u, u, mpfr_bernoulli_cache(n), MPFR_RNDU);
        /* err = 4n+1, and the absolute error is bounded by 4n+1 ulp(u)
           [Rule 11] from algorithms.pdf */
      /* if the terms 'u' are decreasing by a factor two at least,
         then the error coming from those is bounded by
         sum((4n+1)/2^n, n=1..infinity) = 9 */
      exps = MPFR_GET_EXP (s);
      expu = MPFR_GET_EXP (u);
      if (expu < exps - (mpfr_exp_t) p)
        break;
      mpfr_add (s, s, u, MPFR_RNDN);
      if (MPFR_GET_EXP (s) < exps)
        e <<= exps - MPFR_GET_EXP (s);
      e ++; /* error in mpfr_add */
      f = 4 * n + 1;
      /* convert the 4n+1 ulp(u) error into ulp(s) */
      while (expu < exps)
        {
          f = (1 + f) / 2;
          expu ++;
        }
      e += f; /* total rounding error coming from 'u' term */
    }

  mpfr_clear (t);
  mpfr_clear (u);
  mpfr_clear (invxx);

  /* multiply the sum s by 1/x: if the error is bounded by e ulp(s),
     then it is bounded relatively by 2*e*2^-p, thus after the division
     below it is bounded relatively by (1+2*e*2^-p)*(1 + 2^-p) - 1 <
     (2e+2)*2^-p thus by 2e+2 ulps (again by Rule 1). */
  mpfr_div (s, s, x, MPFR_RNDN);

  e = 2 * e + 2;
  f = 0;
  while (e > 1)
    {
      f++;
      e = (e + 1) / 2;
      /* Invariant: 2^f * e does not decrease */
    }
  return f;
}

/* case x >= 1/2 */
static int
mpfr_trigamma_positive (mpfr_ptr y, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  mpfr_prec_t p = MPFR_PREC(y) + 10, q;
  mpfr_t t, u, x_plus_j;
  int inex;
  mpfr_exp_t errt, erru, expt;
  mpfr_prec_t j = 0, min;
  MPFR_ZIV_DECL (loop);

  if (mpfr_cmp_ui_2exp (x, 1, -1) == 0) /* x = 1/2 */
    return half (y, rnd_mode);

  /* now x > 1/2: we use the shift formula trigamma(x+1) = trigamma(x) - 1/x^2
     which yields
     trigamma(x) = 1/x^2 + 1/(x+1)^2 + ... + 1/(x+k-1)^2 + trigamma(x+k)
     until z = x+k is large enough such that we can use the formula:
     trigamma(z) = 1/z + 1/(2z^2) + sum(B[2j]/z^(2j+1), j=1..infinity) (2)
     where B[2j] are Bernoulli numbers.
  */

  /* Compute a precision q such that x+1 is exact. */
  if (MPFR_PREC(x) <= MPFR_GET_EXP(x))
    {
      /* The goal of the first assertion is to let the compiler ignore
         the second one when MPFR_EMAX_MAX <= MPFR_PREC_MAX. */
      MPFR_ASSERTD (MPFR_EXP(x) <= MPFR_EMAX_MAX);
      MPFR_ASSERTN (MPFR_EXP(x) <= MPFR_PREC_MAX);
      /* In that case, ulp(x) = 2^(EXP(x)-PREC(x)) >= 1,
         thus adding 1 will not change the precision (in case of binade
         change, we have x+1 = 2^EXP(x) which is exact). */
      q = MPFR_EXP(x);
    }
  else
    /* In that case, ulp(x) < 1, thus if we add 1 at bit of weight 0,
       we might get an overflow, and need PREC(x)+1 bits. */
    q = MPFR_PREC(x) + 1;

  mpfr_init2 (x_plus_j, q);

  mpfr_init2 (t, p);
  mpfr_init2 (u, p);
  MPFR_ZIV_INIT (loop, p);
  for(;;)
    {
      /* Since |B[2j]| ~ 4*sqrt(pi*j)*(j/(pi*e))^(2j), we have
         t[j] := B[2j]/z^(2j+1) ~ 4/z*sqrt(pi*j)*(j/(pi*e*z))^(2j).
         This yields t[j+1]/t[j] ~ ((j+1)/(pi*e*z))^(2j+2)/(j/(pi*e*z))^(2j)
                                 ~ ((j+1)/(pi*e*z))^2*(1+1/j)^(2j)
                                 ~ (j/(pi*e*z))^2*e^2 ~ (j/(pi*z))^2.
         The smaller term in the divergent series (2) is obtained approximately
         for j = pi*z, and this term is about exp(-2*pi*z).
         Since we want it to be less than 2^-p, this gives z > p*log(2)/(2*pi),
         i.e., x >= 0.1103 p. To be safe, we ensure x >= 0.25 * p.
      */
      min = (p + 3) / 4;
      if (min < 2)
        min = 2; /* ensures x_plus_j >= 2 at the end of the loop below */

      mpfr_set (x_plus_j, x, MPFR_RNDN);
      mpfr_set_ui (u, 0, MPFR_RNDN);
      j = 0;
      while (mpfr_cmp_ui (x_plus_j, min) < 0)
        {
          j ++;
          mpfr_ui_div (t, 1, x_plus_j, MPFR_RNDN);
          /* t = 1/(x+j) * (1 + theta1) with |theta1| < 2^-p */
          mpfr_sqr (t, t, MPFR_RNDN);
          /* t = 1/(x+j)^2 * (1 + theta1)^2 * (1 + theta2) with
             |theta1|, |theta2| < 2^-p, thus
             t = 1/(x+j)^2 * (1 + theta3) with |theta3| < 3.1 * 2^-p */
          mpfr_add (u, u, t, MPFR_RNDN);
          /* u = (u_old + t) * (1 + theta4) with |theta4| < 2^-p
               = (u_old + 1/(x+j)^2 * (1 + theta3)) * (1 + theta4)
               = (u_old + 1/(x+j)^2) * (1 + theta5)
               with |theta5| < (1 + theta3)*(1 + theta4) - 1 < 4.2 * 2^-p.
               The relative error for this step is thus bounded by
               4.2 * 2^-p at each step, thus at most by 4.2 ulps */
          inex = mpfr_add_ui (x_plus_j, x_plus_j, 1, MPFR_RNDZ);
          if (inex != 0) /* we lost one bit */
            {
              q ++;
              mpfr_prec_round (x_plus_j, q, MPFR_RNDZ);
              mpfr_nextabove (x_plus_j);
            }
          /* By induction, the total error is bounded by 4.2*j ulps.
             Indeed, assume the error on u_old was bounded by
             4.2*(j-1)*ulp(u_old), then the total error is bounded by:
             4.2*(j-1)*ulp(u_old) + 4.2*ulp(u) <= 4.2*j*ulp(u) since u_old < u.
          */
        }
    }
  j = 5 * u; /* upper bound for the error */
  for (erru = 0; j > 1; erru++, j = (j + 1) / 2);
  errt = mpfr_trigamma_approx (t, x_plus_j);
  MPFR_ZIV_FREE (loop);
  inex = mpfr_set (y, t, rnd_mode);
  mpfr_clear (t);
  mpfr_clear (u);
  mpfr_clear (x_plus_j);
  return inex;
}

/* trigamma is the 2nd derivative of log(gamma(x)) */
int
mpfr_trigamma (mpfr_ptr y, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  int inex;
  MPFR_SAVE_EXPO_DECL (expo);

  MPFR_LOG_FUNC
    (("x[%Pd]=%.*Rg rnd=%d", mpfr_get_prec(x), mpfr_log_prec, x, rnd_mode),
     ("y[%Pd]=%.*Rg inexact=%d", mpfr_get_prec(y), mpfr_log_prec, y, inex));

  if (MPFR_UNLIKELY(MPFR_IS_SINGULAR(x)))
    {
      if (MPFR_IS_NAN(x))
        {
          MPFR_SET_NAN(y);
          MPFR_RET_NAN;
        }
      else if (MPFR_IS_INF(x))
        {
          if (MPFR_IS_POS(x)) /* trigamma(+Inf) = +0 */
            {
              MPFR_SET_SAME_SIGN(y, x);
              MPFR_SET_ZERO(y);
              MPFR_RET(0);
            }
          else                /* trigamma(-Inf) = NaN */
            {
              MPFR_SET_NAN(y);
              MPFR_RET_NAN;
            }
        }
      else /* Zero case */
        {
          /* the following works also in case of overlap */
          MPFR_SET_INF(y);
          MPFR_SET_POS(y);
          MPFR_SET_DIVBY0 ();
          MPFR_RET(0);
        }
    }

  /* trigamma is undefined for negative integers */
  if (MPFR_IS_NEG(x) && mpfr_integer_p (x))
    {
      MPFR_SET_NAN(y);
      MPFR_RET_NAN;
    }

  /* now x is a normal number */

  MPFR_SAVE_EXPO_MARK (expo);
  /* For x very small, we have trigamma(x) = 1/x^2 + O(1),
     where the O(1) term is less than 2 for |x| < 2^-4.
     Let w = prec(y) + 20 be the working precision.
     If |x| < 2^e, then 1/x^2 > 2^(-2e), thus ulp_w(1/x^2) >= 2^(-2e+1-w).
     As long as -2e+1-w >= -1, we have ulp_w(1/x^2) >= 1/2,
     thus |trigamma(x) - 1/x^2| < 4 ulp_w(1/x^2).
  */
  mpfr_exp_t e = MPFR_GET_EXP (x);
  if (e <= -4) /* |x| < 2^-4 */
    {
      mpfr_prec_t w = MPFR_PREC(y) + 20;
      if (-2 * e + 1 - w >= -1)
        {
          mpfr_t t;
          mpfr_init2 (t, w);
          mpfr_mul (t, x, x, MPFR_RNDN);
          /* t = x^2 * (1 + theta1) with |theta1| < 2^-w */
          mpfr_si_div (t, 1, t, MPFR_RNDN);
          /* t = 1/x^2 / (1 + theta1)^2 * (1 + theta2)
             with |theta1|, |theta2} < 2^-w thus
             t = 1/x^2 * (1 + theta3) with |theta3| < 4*2^-w,
             and the rounding error is bounded by 4 ulps.
             Since the error from the O(1) term is also bounded by 4 ulps,
             the total error is bounded by 8 ulps. */
          if (MPFR_CAN_ROUND (t, p - 3, MPFR_PREC(y), rnd_mode))
            {
              inex = mpfr_set (y, t, rnd_mode);
              mpfr_clear (t);
              MPFR_SAVE_EXPO_UPDATE_FLAGS (expo, __gmpfr_flags);
              goto end;
            }
        }
    }

  if (MPFR_IS_NEG(x) || MPFR_EXP(x) < 0) /* x < 1/2 */
    inex = mpfr_trigamma_reflection (y, x, rnd_mode);
  else
    inex = mpfr_trigamma_positive (y, x, rnd_mode);

 end:
  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (y, inex, rnd_mode);
}
