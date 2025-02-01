/* mpfr_trigamma -- trigamma function of a floating-point number

Copyright 2024-2025 Free Software Foundation, Inc.
Contributed by the Pascaline and Caramba projects, INRIA.

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
   [2] Handbook of Mathematical Functions, Abramowitz & Stegun, 1964

   We assume there are no exact or midpoint cases for trigamma(x),
   which would make Ziv's algorithm loop forever.
*/

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
   trigamma(x) = 1/x * (1 + 1/(2x) + sum(B[2j]/x^(2j), j=1..infinity)).
   This is a divergent series, and as usual we assume the error when we sum up to
   term of index j is bounded by the absolute value of the (j+1)-th term.
*/
static mpfr_exp_t
mpfr_trigamma_approx (mpfr_ptr s, mpfr_srcptr x)
{
  mpfr_prec_t p;
  mpfr_t t, u, invxx;
  mpfr_exp_t e, f, expu;
  unsigned long n;

  MPFR_ASSERTN (MPFR_IS_POS (x) && MPFR_GET_EXP (x) >= 2);

  p = MPFR_GET_PREC (s);
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
         = b[n]*t[n] where t[n]/t[n-1] = 1/(2n)/(2n+1)/x^2 with t[0]=1. */
      mpfr_mul (t, t, invxx, MPFR_RNDU);        /* err = err + 3 */
      mpfr_div_ui (t, t, 2 * n, MPFR_RNDU);     /* err = err + 1 */
      mpfr_div_ui (t, t, 2 * n + 1, MPFR_RNDU); /* err = err + 1 */
      /* we thus have err = 5n here */
      mpfr_mul_z (u, t, mpfr_bernoulli_cache(n), MPFR_RNDU);
        /* err = 5n+1, and the absolute error is bounded by 5n+1 ulp(u)
           [Rule 11] from algorithms.pdf */
      /* if the terms 'u' are decreasing by a factor two at least,
         then the error coming from those is bounded by
         sum((5n+1)/2^n, n=1..infinity) = 11 */
      MPFR_ASSERTD(MPFR_GET_EXP(s) == 1);
      /* since x*trigamma(x) is decreasing on [2,+Inf) from about 1.29 to 1,
         s is always in the binade [1,2) here */
      expu = MPFR_GET_EXP (u);
      if (expu < 1 - (mpfr_exp_t) p) /* |u| < 1/2 ulp(s) since EXP(s) = 1 */
        break;
      mpfr_add (s, s, u, MPFR_RNDN);
      e ++; /* error in mpfr_add */
      f = 5 * n + 1;
      /* convert the 5n+1 ulp(u) error into ulp(s) */
      while (expu < 1)
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
  mpfr_prec_t px, py, p, q;
  mpfr_t t, u, x_plus_j;
  int inex;
  mpfr_exp_t expx, errt, erru, expt1, expt2;
  mpfr_prec_t guard = 10, j, min;
  MPFR_ZIV_DECL (loop);

  if (mpfr_cmp_ui_2exp (x, 1, -1) == 0) /* x = 1/2 */
    return half (y, rnd_mode);

  py = MPFR_GET_PREC (y);
  expx = MPFR_GET_EXP (x);

  /* for very large x, trigamma(x) = 1/x + 1/(2x^2) + O(1/x^3) according
     to formula (6.4.12) from Abramowitz & Stegun. Graphically we see that
     for x >= 1, 1/x + 1/(2x^2) < trigamma(x) < 1/x + 1/x^2. */
  while (py + guard < expx)
    {
      /* this ensures x >= 1, moreover with e := MPFR_PREC(y) + guard,
         e < EXP(x) ensures 2^e <= x since 2^(EXP(x)-1) <= x < 2^EXP(x),
         thus 1/x^2 <= 2^(-2e) */
      mpfr_init2 (t, py + guard);
      inex = mpfr_ui_div (t, 1, x, MPFR_RNDN);
      /* if x is a huge power of 2, then we can round as soon as 1/x^2
         <= 1/2 ulp_py(1/x), where py is the precision of y.
         If x = 2^k, then 1/2 ulp_py(1/x) = 2^(-k-py-1), which
         gives -2k <= -k-py-1, thus py <= k - 1. */
      if (inex == 0 && py <= expx - 2)
        {
          mpfr_set (y, t, rnd_mode);
          mpfr_clear (t);
          if (rnd_mode == MPFR_RNDA || rnd_mode == MPFR_RNDU)
            {
              mpfr_nextabove (y);
              return 1;
            }
          else
            return -1;
        }
      /* |t - trigamma(x)| <= 1/2*ulp(t) + |trigamma(x) - 1/x|
                           <= 1/2*ulp(t) + 1/(2x^2)
                           <= 1/2*ulp(t) + 2^(-2e-1)
                          <= ulp(t)
         since |t| >= 2^-e thus ulp(t) >= 2^(-e-PREC(y)-guard) = 2^(-2e) */
      if (MPFR_CAN_ROUND (t, py + guard, py, rnd_mode))
        {
          inex = mpfr_set (y, t, rnd_mode);
          mpfr_clear (t);
          return inex;
        }
      mpfr_clear (t);
      /* double the guard bits, as long as PREC(y) + guard < EXP(x).
         Note: similar to MPFR_ZIV_NEXT in a Ziv loop. */
      if (py + 2 * guard < expx)
        guard = 2 * guard;
      else if (guard < expx - py - 1)
        guard = expx - py - 1; /* largest possible value */
      else
        break;
    }

  /* now x > 1/2: we use the shift formula trigamma(x+1) = trigamma(x) - 1/x^2
     which yields
     trigamma(x) = 1/x^2 + 1/(x+1)^2 + ... + 1/(x+j)^2 + trigamma(x+j+1)
     until z = x+j+1 is large enough such that we can use the formula:
     trigamma(z) = 1/z + 1/(2z^2) + sum(B[2j]/z^(2j+1), j=1..infinity) (2)
     where B[2j] are Bernoulli numbers.
  */

  px = MPFR_GET_PREC (x);

  /* Compute a precision q such that x+1 is exact. */
  if (px <= expx)
    {
      /* The goal of the first assertion is to let the compiler ignore
         the second one when MPFR_EMAX_MAX <= MPFR_PREC_MAX. */
      MPFR_ASSERTD (expx <= MPFR_EMAX_MAX);
      MPFR_ASSERTN (expx <= MPFR_PREC_MAX);
      /* In that case, ulp(x) = 2^(EXP(x)-PREC(x)) >= 1,
         thus adding 1 will not change the precision (in case of binade
         change, we have x+1 = 2^EXP(x) which is exact). */
      q = expx;
    }
  else
    /* In that case, ulp(x) < 1, thus if we add 1 at bit of weight 0,
       we might get an overflow, and need PREC(x)+1 bits. */
    q = px + 1;

  mpfr_init2 (x_plus_j, q);

  p = py + 10;
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

      /* For the mpfr_cmp_ui below. Currently, mpfr_prec_t <= long,
         so that the compiler should remove this check. */
      MPFR_ASSERTN (min <= ULONG_MAX);

      mpfr_set (x_plus_j, x, MPFR_RNDN);
      mpfr_set_ui (u, 0, MPFR_RNDN);
      j = 0;
      while (mpfr_cmp_ui (x_plus_j, min) < 0)
        {
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
          j ++;
          /* By induction, the total error is bounded by 4.2*j ulps.
             Indeed, assume the error on u_old was bounded by
             4.2*(j-1)*ulp(u_old), then the total error is bounded by:
             4.2*(j-1)*ulp(u_old) + 4.2*ulp(u) <= 4.2*j*ulp(u) since u_old < u.
          */
        }
      /* u approximates 1/x^2 + 1/(x+1)^2 + ... + 1/(x+j-1)^2 */
      j = 5 * j; /* upper bound for the error */
      for (erru = 0; j > 1; erru++, j = (j + 1) / 2);
      errt = mpfr_trigamma_approx (t, x_plus_j);
      expt1 = MPFR_GET_EXP (t);
      /* now u approximates 1/x^2 + ... + 1/(x+j)^2 with error <= 2^erru ulp(u)
         and t approximates 1/z + 1/(2z^2) + sum(B[2j]/z^(2j+1), j=1..infinity)
         for z = x+j+1, with error <= 2^errt ulp(t) */
      mpfr_add (t, u, t, MPFR_RNDN); /* add both terms */
      MPFR_ASSERTD(MPFR_NOTZERO(t));
      expt2 = MPFR_GET_EXP (t);
      /* scale errt in case of cancellation */
      if (expt2 < expt1)
        errt += expt1 - expt2;
      /* scale erru in case of cancellation */
      if (MPFR_NOTZERO(u) && expt2 < MPFR_GET_EXP (u))
        erru += MPFR_EXP(u) - expt2;
      /* the error is now bounded by (2^errt + 2^erru) * ulp(t),
         for the new value of t */
      errt = (errt >= erru ? errt : erru) + 1;
      /* the error is bounded by 2^errt * ulp(t) */
      if (MPFR_CAN_ROUND (t, p - errt, MPFR_PREC(y), rnd_mode))
        break;
      MPFR_ZIV_NEXT (loop, p);
      mpfr_set_prec (t, p);
      mpfr_set_prec (u, p);
    }
  MPFR_ZIV_FREE (loop);
  inex = mpfr_set (y, t, rnd_mode);
  mpfr_clear (t);
  mpfr_clear (u);
  mpfr_clear (x_plus_j);
  return inex;
}

/* compute trigamma(x) for x < 1/2 (x <> 0) using the reflection formula
   (6.4.7) from Abramowitz & Stegun:
   trigamma(1-x) + trigamma(x) = pi^2*(1+cot(pi*x)^2)
   thus for x < 1/2:
   (a) evaluate z = trigamma(y), where y = 1-x >= 1/2
   (b) return pi^2*(1+cot(pi*x)^2) - z
*/
static int
mpfr_trigamma_reflection (mpfr_ptr y, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  mpfr_prec_t py, p;
  mpfr_t t, u, v;
  mpfr_exp_t e1, e2, expt, expv, expx, q;
  int inex;
  MPFR_ZIV_DECL (loop);

  /* trigamma(n) = +Inf for n integer, n <= 0 */
  if (mpfr_integer_p (x))
    {
      MPFR_SET_INF(y);
      MPFR_SET_POS(y);
      return 0;
    }

  py = MPFR_GET_PREC(y);
  p = py + 10;

  /* we want that 1-x is exact with precision q: if 0 < x < 1/2, then
     q = PREC(x)-EXP(x) is ok, otherwise if -1 <= x < 0, q = PREC(x)-EXP(x)+1
     is ok, otherwise for x < -1, PREC(x)+1 is ok if EXP(x) <= PREC(x),
     otherwise we need EXP(x) */
  expx = MPFR_GET_EXP (x);
  if (MPFR_IS_POS (x))           /* 0 < x < 1/2 */
    q = MPFR_GET_PREC(x) - expx;
  else if (expx <= 0)            /* -1/2 < x < 0 */
    q = MPFR_GET_PREC(x) - expx + 1;
  else if (expx <= MPFR_GET_PREC(x))
    q = MPFR_PREC(x) + 1;
  else
    q = expx;
  MPFR_ASSERTN (q <= MPFR_PREC_MAX);
  mpfr_init2 (u, q);
  MPFR_DBGRES(inex = mpfr_ui_sub (u, 1, x, MPFR_RNDN));
  MPFR_ASSERTN(inex == 0);

  mpfr_init2 (t, p);
  mpfr_init2 (v, p);

  MPFR_ZIV_INIT (loop, p);
  for (;;)
    {
      /* below we use theta for a variable with |theta|<=2^-p, where different
         instances of theta may represent different values */
      mpfr_const_pi (v, MPFR_RNDN);  /* v = pi*(1+theta) */
      mpfr_mul (t, v, x, MPFR_RNDN); /* t = pi*x*(1+theta)^2 */
      /* thus t = pi*x*(1+3*theta) (with a different value of theta)
         and the relative error is bounded by 3*2^-p, which by Rule 1 from
         algorithms.pdf translates into 3 ulps(t) */
      e1 = MPFR_GET_EXP(t) - (mpfr_exp_t) p + 2;
      /* bound for t: err(t) <= 2^e1 */

      /* compute cot(t) */
      mpfr_cot (t, t, MPFR_RNDN);
      expt = MPFR_GET_EXP(t);
      /* cot(t + h) = cot(t) + eps * (1 + cot(t)^2) with |eps| <= h <= 2^e1 */
      if (expt > 0) /* |cot(t)| > 1 */
        e1 = e1 + 2 * expt + 1; /* 1 + cot(t)^2 <= 2*cot(t)^2 */
      else
        e1 = e1 + 1; /* |cot(t)| <= 1 thus |1 + cot(t)^2| <= 2 */
      /* now |eps * (1 + cot(t)^2)| <= 2^e1 */
      /* add the rounding error from mpfr_cot, which is 1/2 ulp(t)
         = 2^(EXP(t)-p-1) */
      if (e1 >= expt - p - 1)
        e1 ++;
      else
        e1 = expt - p;
      /* now t = cot(pi*x) + eps with |eps| < 2^e1 */

      /* square t */
      mpfr_sqr (t, t, MPFR_RNDN);
      /* the induced error is 2*eps*cot(pi*x) + eps^2
         <= 2^(e1+1) * 2^expt + 2^(2*e1)
         <= 2^max(expt + 1, e1) + e1 + 1 */
      e1 = (e1 <= expt + 1) ? e1 + expt + 2 : 2 * e1 + 1;
      /* the induced error is bounded by 2^e1 */

      /* add the rounding error from mpfr_sqr */
      expt = MPFR_GET_EXP(t);
      if (e1 >= expt - p - 1)
        e1 ++;
      else
        e1 = expt - p;
      /* now t = cot(pi*x)^2 + eps with |eps| < 2^e1 */

      /* add 1 */
      mpfr_add_ui (t, t, 1, MPFR_RNDN);
      /* add the rounding error from the addition */
      expt = MPFR_GET_EXP(t);
      if (e1 >= expt - p - 1)
        e1 ++;
      else
        e1 = expt - p;
      /* now t = cot(pi*x)^2 + 1 + eps with |eps| < 2^e1 */

      /* multiply by pi */
      mpfr_mul (t, t, v, MPFR_RNDN);
      /* the induced error is |pi*eps| < 2^(e1+2) */
      e1 += 2;
      /* add the rounding error from mpfr_mul */
      expt = MPFR_GET_EXP(t);
      if (e1 >= expt - p - 1)
        e1 ++;
      else
        e1 = expt - p;
      /* now t = pi*(cot(pi*x)^2 + 1) + eps with |eps| < 2^e1 */

      /* multiply again by pi */
      mpfr_mul (t, t, v, MPFR_RNDN);
      /* the induced error is |pi*eps| < 2^(e1+2) */
      e1 += 2;
      /* add the rounding error from mpfr_mul */
      expt = MPFR_GET_EXP(t);
      if (e1 >= expt - p - 1)
        e1 ++;
      else
        e1 = expt - p;
      /* now t = pi^2*(cot(pi*x)^2 + 1) + eps with |eps| < 2^e1 */

      mpfr_trigamma_positive (v, u, MPFR_RNDN);   /* error <= 1/2 ulp */
      expv = MPFR_GET_EXP (v);
      mpfr_sub (v, t, v, MPFR_RNDN);
      if (MPFR_NOTZERO(v))
        {
          /* convert absolute error 2^e1 for t into ulp(v) */
          e1 -= MPFR_EXP(v) - p; /* ulp(v) = 2^(EXP(v) - PREC(v)) */
          /* the error on t is now bounded by 2^e1 * ulp(v) */
          /* now take into account the 1/2 ulp error on old v,
             and the 1/2 ulp error on (new) v */
          if (MPFR_EXP(v) >= expv) /* the new v is larger */
            e2 = 0; /* EXP(old_v) <= EXP(v) thus 1/2 ulp(old_v) <= 1/2 ulp(v)
                       thus 1/2 ulp(old_v) + 1/2 ulp(v) <= 2^0 * ulp(v) */
          else
            e2 = expv - MPFR_EXP(v); /* EXP(v) <= 2^k * EXP(old_v)
                                        thus 1/2 ulp(old_v) + 1/2 ulp(v)
                                        <= (2^(k-1) + 1) * ulp(v)
                                        <= 2^k * ulp(v) */
          /* add both errors */
          e1 = (e1 >= e2 ? e1 : e2) + 1;
          if (MPFR_CAN_ROUND (v, p - e1, py, rnd_mode))
            break;
        }
      MPFR_ZIV_NEXT (loop, p);
      mpfr_set_prec (t, p);
      mpfr_set_prec (v, p);
    }
  MPFR_ZIV_FREE (loop);

  inex = mpfr_set (y, v, rnd_mode);

  mpfr_clear (t);
  mpfr_clear (v);
  mpfr_clear (u);

  return inex;
}

/* trigamma is the 2nd derivative of log(gamma(x)) */
int
mpfr_trigamma (mpfr_ptr y, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  mpfr_exp_t e;
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
              MPFR_SET_ZERO(y);
              MPFR_SET_POS(y);
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

  /* trigamma is +Inf for negative integers */
  if (MPFR_IS_NEG(x) && mpfr_integer_p (x))
    {
      MPFR_SET_INF(y);
      MPFR_SET_POS(y);
      MPFR_SET_DIVBY0 ();
      MPFR_RET(0);
    }

  /* now x is a normal number */

  MPFR_SAVE_EXPO_MARK (expo);
  /* For x very small, we have trigamma(x) = 1/x^2 + O(1),
     where the O(1) term is positive and less than 2 for |x| < 2^-4.
     Let w = prec(y) + 20 be the working precision.
     If |x| < 2^e, then 1/x^2 > 2^(-2e), thus ulp_w(1/x^2) >= 2^(-2e+1-w).
     As long as -2e+1-w >= -1, we have ulp_w(1/x^2) >= 1/2,
     thus |trigamma(x) - 1/x^2| < 4 ulp_w(1/x^2).
  */
  e = MPFR_GET_EXP (x);
  if (e <= -4) /* |x| < 2^-4 */
    {
      int ok = 0;
      mpfr_prec_t w = MPFR_PREC(y) + 20;
      if (-2 * e + 1 - w >= -1)
        {
          mpfr_t t;
          mpfr_init2 (t, w);
          /* t = x^2 * (1 + theta1) with |theta1| < 2^-w */
          inex = mpfr_si_div (t, 1, x, MPFR_RNDN);
          /* if t = o(1/x) overflows with extended exponent range,
             then 1/x^2 will overflow too */
          if (MPFR_IS_INF(t))
            {
              MPFR_LOG_MSG (("t = o(1/x) overflows (rnd=%s)\n",
                             mpfr_print_rnd_mode (rnd_mode)));
              mpfr_clear (t);
              MPFR_SAVE_EXPO_FREE (expo);
              return mpfr_overflow (y, rnd_mode, MPFR_SIGN (x));
            }
          else if (inex == 0) /* x is a power of two: round to 1/x^2 + eps */
            {
              mpfr_sqr (y, t, MPFR_RNDN);
              if (rnd_mode == MPFR_RNDA || rnd_mode == MPFR_RNDU)
                {
                  mpfr_nextabove (y);
                  inex = 1;
                }
              else
                inex = -1;
              ok = 1;
            }
          /* t = 1/x^2 / (1 + theta1)^2 * (1 + theta2)
             with |theta1|, |theta2} < 2^-w thus
             t = 1/x^2 * (1 + theta3) with |theta3| < 4*2^-w,
             and the rounding error is bounded by 4 ulps.
             Since the error from the O(1) term is also bounded by 4 ulps,
             the total error is bounded by 8 ulps. */
          else if (MPFR_CAN_ROUND (t, w - 3, MPFR_PREC(y), rnd_mode))
            {
              inex = mpfr_sqr (y, t, rnd_mode);
              ok = 1;
            }
          mpfr_clear (t);
          if (ok)
            {
              MPFR_SAVE_EXPO_UPDATE_FLAGS (expo, __gmpfr_flags);
              goto end;
            }
        }
    }

  if (MPFR_IS_NEG(x) || e < 0) /* x < 1/2 */
    inex = mpfr_trigamma_reflection (y, x, rnd_mode);
  else
    inex = mpfr_trigamma_positive (y, x, rnd_mode);

 end:
  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (y, inex, rnd_mode);
}
