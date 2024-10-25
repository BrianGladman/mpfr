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
          inex = mpfr_mul (t, x, x, MPFR_RNDN);
          /* t = x^2 * (1 + theta1) with |theta1| < 2^-w */
          inex = inex || mpfr_si_div (t, 1, t, MPFR_RNDN);
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

 end:
  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (y, inex, rnd_mode);
}
