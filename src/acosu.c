/* mpfr_acosu  -- acosu(x)  = acos(x)*u/(2*pi)
   mpfr_acospi -- acospi(x) = acos(x)/pi

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

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* put in y the corrected-rounded value of acos(x)*u/(2*pi) */
int
mpfr_acosu (mpfr_ptr y, mpfr_srcptr x, unsigned long u, mpfr_rnd_t rnd_mode)
{
  mpfr_t tmp, pi;
  mpfr_prec_t prec;
  mpfr_exp_t expx;
  int compared, inexact;
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC
    (("x[%Pu]=%.*Rg u=%lu rnd=%d", mpfr_get_prec(x), mpfr_log_prec, x, u,
      rnd_mode),
     ("y[%Pu]=%.*Rg inexact=%d", mpfr_get_prec (y), mpfr_log_prec, y,
      inexact));

  /* Singular cases */
  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    {
      if (MPFR_IS_NAN (x) || MPFR_IS_INF (x))
        {
          MPFR_SET_NAN (y);
          MPFR_RET_NAN;
        }
      else /* necessarily x=0 */
        {
          MPFR_ASSERTD(MPFR_IS_ZERO(x));
          /* acos(0)=Pi/2 thus acosu(0)=u/4 */
          return mpfr_set_ui_2exp (y, u, -2, rnd_mode);
        }
    }

  compared = mpfr_cmpabs_ui (x, 1);
  if (compared >= 0)
    {
      /* acosu(x) = NaN for |x| > 1 */
      if (compared > 0)
        {
          MPFR_SET_NAN (y);
          MPFR_RET_NAN;
        }
      else /* |x| = 1: acosu(1,u) = 0, acosu(-1,u)=u/2 */
        {
          if (MPFR_SIGN(x) > 0)
            return mpfr_set_ui (y, 0, rnd_mode);
          else
            return mpfr_set_ui_2exp (y, u, -1, rnd_mode);
        }
    }

  /* acos(+/-1/2) = +/-pi/3, thus in this case acos(x,u) is exact when
     u is a multiple of 3 */
  if (mpfr_cmp_si_2exp (x, MPFR_SIGN(x), -1) == 0 && (u % 3) == 0)
    {
      long v = u / 3;
      if (MPFR_IS_NEG (x))
        v = -v;
      return mpfr_set_si_2exp (y, v, -1, rnd_mode);
    }

  prec = MPFR_PREC (y);

  /* For |x|<0.5, we have acos(x) = pi/2 - x*r(x) with |r(x)| < 1.05
     thus acosu(x,u) = u/4*(1 - x*s(x)) with 0 <= s(x) < 1.
     If EXP(x) <= -prec-3, then |u/4*x*s(x)| < u/4*2^(-prec-3) < ulp(u/4)/8
     <= ulp(RN(u/4))/4, thus the result will be u/4, nextbelow(u/4) or
     nextabove(u/4).
     Warning: when u/4 is a power of two, the difference between u/4 and
     nextbelow(u/4) is only 1/4*ulp(u/4).
     We also require x < 2^-64, so that in the case u/4 is not exact,
     the contribution of x*s(x) is smaller compared to the last bit of u. */
  expx = MPFR_GET_EXP(x);
  if (expx <= -64 && expx <= - (mpfr_exp_t) prec - 3)
    {
      MPFR_SET_INEXFLAG ();
      inexact = mpfr_set_ui_2exp (y, u, -2, rnd_mode);
      /* for all rounding modes, if the division u/4 is inexact, it will
         give the correct rounding */
      if (inexact != 0 || rnd_mode == MPFR_RNDF)
        {
          /* For RNDN, if inexact = +/-MPFR_EVEN_INEX, we can't conclude:
             inexact=MPFR_EVEN_INEX x>0: subtract one ulp
             inexact=MPFR_EVEN_INEX x<0: ok
             inexact=-MPFR_EVEN_INEX x>0: ok
             inexact=MPFR_EVEN_INEX x<0: add one ulp */
          if (inexact == MPFR_EVEN_INEX && MPFR_SIGN(x) > 0)
            goto subtract_one_ulp;
          if (inexact == -MPFR_EVEN_INEX && MPFR_SIGN(x) < 0)
            goto add_one_ulp;
          return inexact;
        }
      if (rnd_mode == MPFR_RNDN)
        {
          /* if u/4 is exact, result is u/4-eps for x>0, u/4+eps for x<0 */
          return MPFR_SIGN(x);
        }
      if (MPFR_SIGN(x) > 0) /* exact result is u/4-eps */
        {
          if (MPFR_IS_LIKE_RNDU(rnd_mode,1))
            return +1;
          else /* round down */
            {
            subtract_one_ulp:
              mpfr_nextbelow (y);
              return -1;
            }
        }
      else /* x < 0: exact result is u/4+eps */
        {
          if (MPFR_IS_LIKE_RNDD(rnd_mode,1))
            return -1;
          else /* round up */
            {
            add_one_ulp:
              mpfr_nextabove (y);
              return +1;
            }
        }
    }

  MPFR_SAVE_EXPO_MARK (expo);

  prec += MPFR_INT_CEIL_LOG2(prec) + 10;

  mpfr_init2 (tmp, prec);
  mpfr_init2 (pi, prec);

  MPFR_ZIV_INIT (loop, prec);
  for (;;)
    {
      /* In the error analysis below, each thetax denotes a variable such that
         |thetax| <= 2^-prec */
      mpfr_acos (tmp, x, MPFR_RNDN);
      /* tmp = acos(x) * (1 + theta1) */
      mpfr_const_pi (pi, MPFR_RNDN);
      /* pi = Pi * (1 + theta2) */
      mpfr_div (tmp, tmp, pi, MPFR_RNDN);
      /* tmp = acos(x)/Pi * (1 + theta3)^3 */
      mpfr_mul_ui (tmp, tmp, u, MPFR_RNDN);
      /* tmp = acos(x)*u/Pi * (1 + theta4)^4 */
      mpfr_div_2ui (tmp, tmp, 1, MPFR_RNDN); /* exact */
      /* tmp = acos(x)*u/(2*Pi) * (1 + theta4)^4 */
      /* since |(1 + theta4)^4 - 1| <= 8*|theta4| for prec >= 2,
         the relative error is less than 2^(3-prec) */
      if (MPFR_LIKELY (MPFR_CAN_ROUND (tmp, prec - 3,
                                       MPFR_PREC (y), rnd_mode)))
        break;
      MPFR_ZIV_NEXT (loop, prec);
      mpfr_set_prec (tmp, prec);
      mpfr_set_prec (pi, prec);
    }
  MPFR_ZIV_FREE (loop);

  inexact = mpfr_set (y, tmp, rnd_mode);
  mpfr_clear (tmp);
  mpfr_clear (pi);

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (y, inexact, rnd_mode);
}

int
mpfr_acospi (mpfr_ptr y, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  return mpfr_acosu (y, x, 2, rnd_mode);
}
