/* mpfr_gamma_inc -- incomplete gamma function

Copyright 2016 Free Software Foundation, Inc.
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

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* The incomplete gamma function is defined as:

   gamma_inc(a,x) := Gamma(a,x) = int(t^(a-1) * exp(-t), t=x..infinity)

                  = Gamma(a) - gamma(a,x) with:

   gamma(a,x) = int(t^(a-1) * exp(-t), t=0..x).

   The function gamma(a,x) satisfies the Taylor expansions:

   gamma(a,x) = x^a * sum((-x)^k/k!/(a+k), k=0..infinity)

   gamma(a,x) = x^a * exp(-x) * sum(x^k/(a*(a+1)*...*(a+k)), k=0..infinity) */

int
mpfr_gamma_inc (mpfr_ptr y, mpfr_srcptr a, mpfr_srcptr x, mpfr_rnd_t rnd)
{
  MPFR_ASSERTN(mpfr_regular_p (a));
  MPFR_ASSERTN(mpfr_regular_p (x));
  MPFR_ASSERTN(MPFR_SIGN(a) > 0);
  MPFR_ASSERTN(MPFR_SIGN(x) > 0);

  mpfr_prec_t w;
  mpfr_t s, t, u;
  int inex;
  unsigned long k, err;
  mpfr_exp_t e0, e1, e2;
  MPFR_GROUP_DECL(group);
  MPFR_ZIV_DECL(loop);

  w = MPFR_PREC(y) + 13; /* working precision */

  MPFR_GROUP_INIT_2(group, w, s, t);
  mpfr_init2 (u, w); /* u is special (see below) */
  MPFR_ZIV_INIT (loop, w);
  for (;;)
    {
      /* Note: in the error analysis below, theta represents any value of
         absolute value less than 2^(-w) where w is the working precision (two
         instances of theta may represent different values), cf Higham's book.
      */

      /* to ensure that u = a + k is exact, we require that ulp(u) <= 1 */
      if (MPFR_EXP(a) > w)
        mpfr_set_prec (u, MPFR_EXP(u));

      /* estimate Taylor series */
      mpfr_ui_div (t, 1, a, MPFR_RNDZ); /* t = 1/a * (1 + theta) */
      mpfr_set (s, t, MPFR_RNDZ);       /* s = 1/a * (1 + theta) */
      for (k = 1;; k++)
        {
          mpfr_mul (t, t, x, MPFR_RNDZ); /* t = x^k/(a * ... * (a+k-1))
                                          * (1 + theta)^(2k) */
          inex = mpfr_add_ui (u, a, k, MPFR_RNDZ); /* u = a+k exactly */
          MPFR_ASSERTD(inex == 0);
          mpfr_div (t, t, u, MPFR_RNDZ); /* t = x^k/(a * ... * (a+k))
                                          * (1 + theta)^(2k+1) */
          mpfr_add (s, s, t, MPFR_RNDZ);
          /* we stop when |t| < ulp(s) and |x/u| < 1/2, which ensures
             that the tail is at most 2*ulp(s) */
          if (MPFR_EXP(t) + w <= MPFR_EXP(s) && MPFR_EXP(x) + 1 < MPFR_EXP(u))
            break;
        }
      /* since all terms are positive, we have s = S * (1 + theta)^(2k+3)
         with S being the infinite Taylor series */

      /* multiply by exp(-x) */
      mpfr_exp (t, x, MPFR_RNDZ);    /* t = exp(x) * (1+theta) */
      mpfr_div (s, s, t, MPFR_RNDZ); /* s = <exact value> * (1+theta)^(2k+5) */

      /* multiply by x^a */
      mpfr_pow (t, x, a, MPFR_RNDZ); /* t = x^a * (1+theta) */
      mpfr_mul (s, s, t, MPFR_RNDZ); /* s = Gamma(a,x) * (1+theta)^(2k+7) */

      /* Since |theta| < 2^(-w) using the Taylor expansion of log(1+x)
         we have log(1+theta) = theta1 with |theta1| < 1.16*2^(-w) for w >= 2,
         thus (1+theta)^(2k+7) = exp((2k+7)*theta1).
         Assuming 2k+7 = t*2^w for |t| < 0.5, we have
         |(2k+7)*theta1| = |t*2^w*theta1| < 0.58.
         For |u| < 0.58 we have |exp(u)-1| < 1.36*|u|
         thus |(1+theta)^(2k+7) - 1| < 1.36*0.58*(2k+7)/2^w < 0.79*(2k+7)/2^w.
         Since one ulp is at worst a relative error of 2^(1-w),
         the error on s is at most 2*(2k+7) ulps. */

      /* subtract from gamma(a) */
      mpfr_gamma (t, a, MPFR_RNDZ);  /* t = gamma(a) * (1+theta) */
      e0 = MPFR_EXP (t);
      e1 = MPFR_EXP (s);
      mpfr_sub (s, t, s, MPFR_RNDZ);
      e2 = MPFR_EXP (s);
      /* the final error is at most 1 ulp (for the final subtraction)
         + 1 ulp * 2^(e0-e2) # for the error in t
         + 2*(2k+7) ulps * 2^(e1-e2) # for the error in gamma(a,x) */

      e1 += 1 + MPFR_INT_CEIL_LOG2 (2*k+7);
      /* Now the error is <= 1 + 2^(e0-e2) + 2^(e1-e2).
         Assume e0 > e1, then it is <= 1 + 1.5*2^(e0-e2)
                                    <= 2^(e0-e2+1) if e0 > e2
                                    <= 2^2 otherwise */
      if (e0 == e1)
        err = e0 - e2 + 2;
      else
        {
          e0 = (e0 > e1) ? e0 : e1; /* max(e0,e1) */
          err = (e0 > e2) ? e0 - e2 + 1 : 2;
        }

      if (MPFR_LIKELY (MPFR_CAN_ROUND (s, w - err, MPFR_PREC(y), rnd)))
        break;

      MPFR_ZIV_NEXT (loop, w);
      MPFR_GROUP_REPREC_2(group, w, s, t);
    }
  MPFR_ZIV_FREE (loop);
  mpfr_clear (u);

  inex = mpfr_set (y, s, rnd);
  MPFR_GROUP_CLEAR(group);

  return inex;
}
