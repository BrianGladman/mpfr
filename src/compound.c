/* mpfr_compound --- compound(x,n) = (1+x)^n

Copyright 2021 Free Software Foundation, Inc.
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
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
https://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#define MPFR_NEED_LONGLONG_H /* needed for MPFR_INT_CEIL_LOG2 */
#include "mpfr-impl.h"

/* assuming |(1+x)^n - 1| < 1/4*ulp(1), return correct rounding,
   where s is the sign of n*log2(1+x) */
static int
mpfr_compound_near_one (mpfr_ptr y, int s, mpfr_rnd_t rnd_mode)
{
  mpfr_set_ui (y, 1, rnd_mode); /* exact */
  if (rnd_mode == MPFR_RNDN || MPFR_RNDF
      || (s > 0 && (rnd_mode == MPFR_RNDZ || rnd_mode == MPFR_RNDD))
      || (s < 0 && (rnd_mode == MPFR_RNDA || rnd_mode == MPFR_RNDU)))
    {
      /* round toward 1 */
      return -s;
    }
  else if (s > 0) /* necessarily RNDA or RNDU */
    {
      /* round toward +Inf */
      mpfr_nextabove (y);
      return +1;
    }
  else /* necessarily s < 0 and RNDZ or RNDD */
    {
      /* round toward 0 */
      mpfr_nextbelow (y);
      return -1;
    }
}

/* put in y the correctly rounded value of (1+x)^n */
int
mpfr_compound (mpfr_ptr y, mpfr_srcptr x, long n, mpfr_rnd_t rnd_mode)
{
  int inexact, compared, k, nloop;
  mpfr_t t;
  mpfr_exp_t e;
  mpfr_prec_t prec;
  MPFR_ZIV_DECL (loop);
  MPFR_SAVE_EXPO_DECL (expo);

  MPFR_LOG_FUNC
    (("x[%Pu]=%.*Rg n=%ld rnd=%d",
      mpfr_get_prec(x), mpfr_log_prec, x, n, rnd_mode),
     ("y[%Pu]=%.*Rg inexact=%d", mpfr_get_prec (y), mpfr_log_prec, y, inexact));

  /* Special cases */
  if (MPFR_IS_SINGULAR (x))
    {
      if (MPFR_IS_NAN (x) || (MPFR_IS_INF (x) && MPFR_SIGN (x) < 0))
        {
          MPFR_SET_NAN (y);
          MPFR_RET_NAN;
        }
      else if (n == 0 || MPFR_IS_ZERO (x))
        {
          /* (1+Inf)^0 = 1 and (1+x)^0 = 1 */
          return mpfr_set_ui (y, 1, rnd_mode);
        }
      else if (MPFR_IS_INF (x)) /* x = +Inf */
        {
          if (n < 0) /* (1+Inf)^n = +0 for n < 0 */
            MPFR_SET_ZERO (y);
          else /* n > 0: (1+Inf)^n = +Inf */
            MPFR_SET_INF (y);
          MPFR_SET_POS (y);
          MPFR_RET (0); /* exact 0 or infinity */
        }
    }

  /* (1+x)^n = NaN for x < -1 */
  compared = mpfr_cmp_si (x, -1);
  if (compared < 0)
    {
      MPFR_SET_NAN (y);
      MPFR_RET_NAN;
    }

  /* compound(x,0) gives 1 for x >= 1 */
  if (n == 0)
    return mpfr_set_ui (y, 1, rnd_mode);

  if (compared == 0)
    {
      if (n < 0)
        {
          /* compound(-1,n) = +Inf with divide-by-zero exception */
          MPFR_SET_INF (y);
          MPFR_SET_POS (y);
          MPFR_SET_DIVBY0 ();
          MPFR_RET (0);
        }
      else
        {
          /* compound(-1,n) = +0 */
          MPFR_SET_ZERO (y);
          MPFR_SET_POS (y);
          MPFR_RET (0);
        }
    }

  MPFR_SAVE_EXPO_MARK (expo);

  prec = MPFR_PREC(y);
  prec += MPFR_INT_CEIL_LOG2 (prec) + 6;

  mpfr_init2 (t, prec);

  k = MPFR_INT_CEIL_LOG2((n > 0) ? n : -n);

  MPFR_ZIV_INIT (loop, prec);
  for (nloop = 0; ; nloop++)
    {
      /* we compute (1+x)^n as 2^(n*log2p1(x)) */
      inexact = mpfr_log2p1 (t, x, MPFR_RNDN);
      e = MPFR_GET_EXP(t);
      /* |t - log2(1+x)| <= 1/2*ulp(t) = 2^(e-prec-1) */
      inexact |= mpfr_mul_ui (t, t, n, MPFR_RNDN);
      /* |t - n*log2(1+x)| <= 2^(e2-prec-1) + n*2^(e-prec-1)
                           <= 2^(e2-prec-1) + 2^(e+k-prec-1) <= 2^(e+k-prec)
                           where n < 2^k, and e2 is the new exponent of t. */
      MPFR_ASSERTD(MPFR_GET_EXP(t) <= e + k);
      e += k;
      /* |t - n*log2(1+x)| <= 2^(e-prec) */
      /* detect overflow */
      if (nloop == 0 && mpfr_cmp_ui (t, __gmpfr_emax) >= 0)
        {
          mpfr_clear (t);
          MPFR_SAVE_EXPO_FREE (expo);
          return mpfr_overflow (y, rnd_mode, 1);
        }
      /* Detect cases where result is 1 or 1+ulp(1) or 1-1/2*ulp(1):
         |2^t - 1| = |exp(t*log(2)) - 1| <= |t|*log(2) < |t| */
      if (nloop == 0 && MPFR_GET_EXP(t) < - (mpfr_exp_t) MPFR_PREC(y))
        {
          int signt = MPFR_SIGN(t);
          /* since ulp(1) = 2^(1-PREC(y)), we have |t| < 1/4*ulp(1) */
          mpfr_clear (t);
          MPFR_SAVE_EXPO_FREE (expo);
          return mpfr_compound_near_one (y, signt, rnd_mode);
        }
      inexact |= mpfr_exp2 (t, t, MPFR_RNDA);
      /* |t - (1+x)^n| <= ulp(t) + |t|*log(2)*2^(e-prec)
                       < 2^(EXP(t)-prec) + 2^(EXP(t)+e-prec) */
      e = (e >= 0) ? e + 1 : 1;
      /* now |t - (1+x)^n| < 2^(EXP(t)+e-prec) */

      if (MPFR_LIKELY (inexact == 0 ||
                       MPFR_CAN_ROUND (t, prec - e, MPFR_PREC(y), rnd_mode)))
        break;

      MPFR_ZIV_NEXT (loop, prec);
      mpfr_set_prec (t, prec);
    }
  inexact = mpfr_set (y, t, rnd_mode);

  MPFR_ZIV_FREE (loop);
  mpfr_clear (t);

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (y, inexact, rnd_mode);
}
