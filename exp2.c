/* mpfr_exp2 -- power of 2 function 2^y

Copyright 2001, 2002, 2003, 2004, 2005 Free Software Foundation.

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
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <limits.h>

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* The computation of y = 2^z is done by                           *
 *     y = exp(z*log(2)). The result is exact iff z is an integer. */

int
mpfr_exp2 (mpfr_ptr y, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  int inexact;
  MPFR_SAVE_EXPO_DECL (expo);

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    {
      if (MPFR_IS_NAN (x))
        {
          MPFR_SET_NAN (y);
          MPFR_RET_NAN;
        }
      else if (MPFR_IS_INF (x))
        {
          if (MPFR_IS_POS (x))
            MPFR_SET_INF (y);
          else
            MPFR_SET_ZERO (y);
          MPFR_SET_POS (y);
          MPFR_RET (0);
        }
      else /* 2^0 = 1 */
        {
          MPFR_ASSERTD (MPFR_IS_ZERO(x));
          return mpfr_set_ui (y, 1, rnd_mode);
        }
    }

  /* since the smallest representable non-zero float is 1/2*2^__gmpfr_emin,
     if x < __gmpfr_emin - 1, the result is either 1/2*2^__gmpfr_emin or 0 */
  MPFR_ASSERTD (MPFR_EMIN_MIN - 2 >= LONG_MIN);

  if (mpfr_cmp_si_2exp (x, __gmpfr_emin - 1, 0) < 0)
    {
      mp_rnd_t rnd2 = rnd_mode;
      /* in round to nearest mode, round to zero when x <= __gmpfr_emin-2 */
      if (rnd_mode == GMP_RNDN &&
          mpfr_cmp_si_2exp (x, __gmpfr_emin - 2, 0) <= 0)
        rnd2 = GMP_RNDZ;
      return mpfr_underflow (y, rnd2, 1);
    }

  if (mpfr_integer_p (x)) /* we know that x >= 2^(emin-1) */
    {
      long xd;

      MPFR_ASSERTD (MPFR_EMAX_MAX <= LONG_MAX);
      if (mpfr_cmp_si_2exp (x, __gmpfr_emax, 0) > 0)
        return mpfr_overflow (y, rnd_mode, 1);

      xd = mpfr_get_si (x, GMP_RNDN);

      mpfr_set_ui (y, 1, GMP_RNDZ);
      return mpfr_mul_2si (y, y, xd, rnd_mode);
    }

  MPFR_SAVE_EXPO_MARK (expo);

  /* General case */
  {
    /* Declaration of the intermediary variable */
    mpfr_t t;

    /* Declaration of the size variable */
    mp_prec_t Ny = MPFR_PREC(y);              /* target precision */
    mp_prec_t Nt;                             /* working precision */
    mp_exp_t err;                             /* error */
    MPFR_ZIV_DECL (loop);

    /* compute the precision of intermediary variable */
    /* the optimal number of bits : see algorithms.tex */
    Nt = Ny + 5 + MPFR_INT_CEIL_LOG2 (Ny);

    /* initialise of intermediary       variable */
    mpfr_init2 (t, Nt);

    /* First computation */
    MPFR_ZIV_INIT (loop, Nt);
    for (;;)
      {
        /* compute   exp(x*ln(2))*/
        mpfr_const_log2 (t, GMP_RNDU);       /* ln(2) */
        mpfr_mul (t, x, t, GMP_RNDU);        /* x*ln(2) */
        err = Nt - (MPFR_GET_EXP (t) + 2);   /* Estimate of the error */
        mpfr_exp (t, t, GMP_RNDN);           /* exp(x*ln(2))*/

        if (MPFR_LIKELY (MPFR_CAN_ROUND (t, err, Ny, rnd_mode)))
          break;

        /* Actualisation of the precision */
        MPFR_ZIV_NEXT (loop, Nt);
        mpfr_set_prec (t, Nt);
      }
    MPFR_ZIV_FREE (loop);

    inexact = mpfr_set (y, t, rnd_mode);

    mpfr_clear (t);
  }
  MPFR_SAVE_EXPO_FREE (expo);

  return mpfr_check_range (y, inexact, rnd_mode);
}
