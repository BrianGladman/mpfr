/* mpfr_erfc -- The Complementary Error Function of a floating-point number

Copyright 2005 Free Software Foundation, Inc.

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

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* erfc(x) = 1 - erf(x) */

int 
mpfr_erfc (mpfr_ptr y, mpfr_srcptr x, mp_rnd_t rnd)
{
  int inex;
  mpfr_t tmp;
  mp_exp_t te, err;
  mp_prec_t prec;
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC (("x[%#R]=%R rnd=%d", x, x, rnd),
                 ("y[%#R]=%R inexact=%d", y, y, inex));

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    {
      if (MPFR_IS_NAN (x))
        {
          MPFR_SET_NAN (y);
          MPFR_RET_NAN;
        }
      /* erfc(+inf) = 0+, erfc(-inf) = 2 erfc (0) = 1 */
      else if (MPFR_IS_INF (x))
        return mpfr_set_ui (y, MPFR_IS_POS (x) ? 0 : 2, rnd);
      else 
        return mpfr_set_ui (y, 1, rnd);
    }

  /* Init stuff */
  MPFR_SAVE_EXPO_MARK (expo);
  prec = MPFR_PREC (y) + MPFR_INT_CEIL_LOG2 (MPFR_PREC (y)) + 3;
  mpfr_init2 (tmp, prec);

  MPFR_ZIV_INIT (loop, prec);            /* Initialize the ZivLoop controler */
  for (;;)                               /* Infinite loop */
    {
      mpfr_erf (tmp, x, GMP_RNDN);
      MPFR_ASSERTD (!MPFR_IS_SINGULAR (tmp)); /* FIXME: 0 only for x=0 ? */
      te = MPFR_GET_EXP (tmp);
      mpfr_ui_sub (tmp, 1, tmp, GMP_RNDN);
      /* See error analysis of expm1 for details */
      if (MPFR_IS_ZERO (tmp))
        prec *=2;
      else
        {
          err = prec - (MAX (te - MPFR_GET_EXP (tmp), 0) + 1);
          if (MPFR_LIKELY (MPFR_CAN_ROUND (tmp, err, MPFR_PREC (y), rnd)))
            break;
        }
      MPFR_ZIV_NEXT (loop, prec);        /* Increase used precision */
      mpfr_set_prec (tmp, prec);
    }
  MPFR_ZIV_FREE (loop);                  /* Free the ZivLoop Controler */

  inex = mpfr_set (y, tmp, rnd);    /* Set y to the computed value */
  mpfr_clear (tmp);

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (y, inex, rnd);
}
