/* mpfr_cosu -- cosu(x) = cos(2*pi*x/u)

Copyright 2020 Free Software Foundation, Inc.
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

/* FIXME[VL]: Implement the range reduction in this function.
   That's the whole point of cosu compared to cos. */

/* put in y the corrected-rounded value of cos(2*pi*x/u) */
int
mpfr_cosu (mpfr_ptr y, mpfr_srcptr x, unsigned long u, mpfr_rnd_t rnd_mode)
{
  mpfr_prec_t precy, prec;
  mpfr_exp_t expx, expt, err;
  mpfr_t t;
  int inexact = 0, nloops = 0, underflow = 0;
  MPFR_ZIV_DECL (loop);
  MPFR_SAVE_EXPO_DECL (expo);

  if (u == 0 || MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    {
      /* for u=0, return NaN */
      if (u == 0 || MPFR_IS_NAN (x) || MPFR_IS_INF (x))
        {
          MPFR_SET_NAN (y);
          MPFR_RET_NAN;
        }
      else /* x is zero: cos(0) = 1 */
        {
          MPFR_ASSERTD (MPFR_IS_ZERO (x));
          return mpfr_set_ui (y, 1, rnd_mode);
        }
    }

  MPFR_SAVE_EXPO_MARK (expo);

  precy = MPFR_PREC (y);
  expx = MPFR_GET_EXP (x);
  /* for x large, since argument reduction is expensive, we want to avoid
     any failure in Ziv's strategy, thus we take into account expx too */
  prec = precy + MPFR_INT_CEIL_LOG2 (MAX(precy,expx)) + 8;
  MPFR_ASSERTD(prec >= 2);
  mpfr_init2 (t, prec);
  MPFR_ZIV_INIT (loop, prec);
  for (;;)
    {
      nloops ++;
      /* We first compute an approximation t of 2*pi*x/u, then call cos(t).
         If t = 2*pi*x/u + s, then |cos(t) - cos(2*pi*x/u)| <= |s|. */
      mpfr_set_prec (t, prec);
      mpfr_const_pi (t, MPFR_RNDN); /* t = pi * (1 + theta1) where
                                       |theta1| <= 2^-prec */
      mpfr_mul_2ui (t, t, 1, MPFR_RNDN); /* t = 2*pi * (1 + theta1) */
      mpfr_mul (t, t, x, MPFR_RNDN);     /* t = 2*pi*x * (1 + theta2)^2 where
                                            |theta2| <= 2^-prec */
      mpfr_div_ui (t, t, u, MPFR_RNDN);  /* t = 2*pi*x/u * (1 + theta3)^3 where
                                            |theta3| <= 2^-prec */
      /* if t is zero here, it means the division by u underflowd */
      if (MPFR_UNLIKELY (MPFR_IS_ZERO (t)))
        {
          mpfr_set_ui (y, 1, MPFR_RNDZ);
          if (MPFR_IS_LIKE_RNDZ(rnd_mode,0))
            {
              inexact = -1;
              mpfr_nextbelow (y);
            }
          else
            inexact = 1;  
          goto end;
        }
      /* since prec >= 2, |(1 + theta3)^3 - 1| <= 4*theta3 <= 2^(2-prec) */
      expt = MPFR_GET_EXP (t);
      /* we have |s| <= 2^(expt + 2 - prec) */
      mpfr_cos (t, t, MPFR_RNDN);
      err = expt + 2 - prec;
      expt = MPFR_GET_EXP (t); /* new exponent of t */
      /* the total error is at most 2^err + ulp(t)/2 = 2^err + 2^(expt-prec-1)
         thus if err <= expt-prec-1, it is bounded by 2^(expt-prec),
         otherwise it is bounded by 2^(err+1). */
      err = (err <= expt - prec - 1) ? expt - prec : err + 1;
      /* normalize err for mpfr_can_round */
      err = expt - err;
      if (MPFR_CAN_ROUND (t, err, precy, rnd_mode))
        break;
      /* check exact cases: this can only occur if 2*pi*x/u is a multiple
         of pi/2, i.e., if x/u is a multiple of 1/4 */
      if (nloops == 1)
        {
          inexact = mpfr_div_ui (t, x, u, MPFR_RNDZ);
          mpfr_mul_2ui (t, t, 2, MPFR_RNDZ);
          if (inexact == 0 && mpfr_integer_p (t))
            {
              if (mpfr_odd_p (t))
                /* t is odd: we have kpi+pi/2, thus cosu = 0,
                   for the sign, we always return +0, following IEEE 754-2019:
                   cosPi(n + 1/2) is +0 for any integer n when n + 1/2 is
                   representable. */
                mpfr_set_zero (y, +1);
              else /* t is even: case kpi */
                {
                  mpfr_div_2ui (t, t, 1, MPFR_RNDZ);
                  if (!mpfr_odd_p (t))
                    /* case 2kpi: cosu = 1 */
                    mpfr_set_ui (y, 1, MPFR_RNDZ);
                  else
                    mpfr_set_si (y, -1, MPFR_RNDZ);
                }
              goto end;
            }
        }
      MPFR_ZIV_NEXT (loop, prec);
    }
  MPFR_ZIV_FREE (loop);

  inexact = mpfr_set (y, t, rnd_mode);

 end:
  mpfr_clear (t);
  MPFR_SAVE_EXPO_FREE (expo);
  return underflow ? inexact : mpfr_check_range (y, inexact, rnd_mode);
}
