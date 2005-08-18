/* mpfr_asinh -- inverse hyperbolic sine

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

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* The computation of asinh is done by  *
 *    asinh = ln(x + sqrt(x^2 + 1))     */

int
mpfr_asinh (mpfr_ptr y, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  int inexact;
  int signx, neg;
  mp_prec_t Ny, Nt;
  mpfr_t t; /* auxiliary variables */
  mp_exp_t err;
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC (("x[%#R]=%R rnd=%d", x, x, rnd_mode),
                 ("y[%#R]=%R inexact=%d", y, y, inexact));

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    {
      if (MPFR_IS_NAN (x))
        {
          MPFR_SET_NAN (y);
          MPFR_RET_NAN;
        }
      else if (MPFR_IS_INF (x))
        {
          MPFR_SET_INF (y);
          MPFR_SET_SAME_SIGN (y, x);
          MPFR_RET (0);
        }
      else /* x is necessarily 0 */
        {
          MPFR_ASSERTD (MPFR_IS_ZERO (x));
          MPFR_SET_ZERO (y);   /* asinh(0) = 0 */
          MPFR_SET_SAME_SIGN (y, x);
          MPFR_RET (0);
        }
    }

  /* asinh(x) = x - x^3/6 + ... so the error is < 2^(3*EXP(x)-2) */
  MPFR_FAST_COMPUTE_IF_SMALL_INPUT (y, x, -2*MPFR_GET_EXP (x)+2,0,rnd_mode,);

  Ny = MPFR_PREC (y);   /* Precision of output variable */

  signx = MPFR_SIGN (x);
  neg = MPFR_IS_NEG (x);

  /* General case */

  /* compute the precision of intermediary variable */
  /* the optimal number of bits : see algorithms.tex */
  Nt = Ny + 4 + MPFR_INT_CEIL_LOG2 (Ny);

  MPFR_SAVE_EXPO_MARK (expo);

  /* initialize intermediary variables */
  mpfr_init2 (t, Nt);

  /* First computation of asinh */
  MPFR_ZIV_INIT (loop, Nt);
  for (;;)
    {
      /* compute asinh */
      mpfr_mul (t, x, x, GMP_RNDD);                    /* x^2 */
      mpfr_add_ui (t, t, 1, GMP_RNDD);                 /* x^2+1 */
      mpfr_sqrt (t, t, GMP_RNDN);                      /* sqrt(x^2+1) */
      (neg ? mpfr_sub : mpfr_add) (t, t, x, GMP_RNDN); /* sqrt(x^2+1)+x */
      mpfr_log (t, t, GMP_RNDN);                       /* ln(sqrt(x^2+1)+x)*/

      /* error estimate -- see algorithms.ps */
      err = Nt - (MAX (3 - MPFR_GET_EXP (t), 0) + 1);

      if (MPFR_LIKELY (MPFR_IS_ZERO (t)
                       || MPFR_CAN_ROUND (t, err, Ny, rnd_mode)))
        break;

      /* actualisation of the precision */
      MPFR_ZIV_NEXT (loop, Nt);
      mpfr_set_prec (t, Nt);
    }
  MPFR_ZIV_FREE (loop);

  inexact = mpfr_set4 (y, t, rnd_mode, signx);

  mpfr_clear (t);

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (y, inexact, rnd_mode);
}
