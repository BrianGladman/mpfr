/* mpfr_cosh -- hyperbolic cosine

Copyright 2001, 2002, 2004, 2005 Free Software Foundation.

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
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* The computation of cosh is done by    *
 *  cosh= 1/2[e^(x)+e^(-x)]              */

int
mpfr_cosh (mpfr_ptr y, mpfr_srcptr xt , mp_rnd_t rnd_mode)
{
  mpfr_t x;
  int inexact;
  MPFR_SAVE_EXPO_DECL (expo);

  MPFR_LOG_FUNC (("x[%#R]=%R rnd=%d", xt, xt, rnd_mode),
		 ("y[%#R]=%R inexact=%d", y, y, inexact));

  if (MPFR_UNLIKELY(MPFR_IS_SINGULAR(xt)))
    {
      if (MPFR_IS_NAN(xt))
	{
	  MPFR_SET_NAN(y);
	  MPFR_RET_NAN;
	}
      else if (MPFR_IS_INF(xt))
	{
	  MPFR_SET_INF(y);
	  MPFR_SET_POS(y);
	  MPFR_RET(0);
	}
      else
	{
	  MPFR_ASSERTD(MPFR_IS_ZERO(xt));
	  return mpfr_set_ui (y, 1, rnd_mode); /* cosh(0) = 1 */
	}
    }

  MPFR_SAVE_EXPO_MARK (expo);
  MPFR_TMP_INIT_ABS(x, xt);
  /* General case */
  {
    /* Declaration of the intermediary variable */
    mpfr_t t, te;
    /* Declaration of the size variable */
    mp_prec_t Ny = MPFR_PREC(y);   /* Precision of output variable */
    mp_prec_t Nt;                  /* Precision of the intermediary variable */
    long int err;                  /* Precision of error */
    MPFR_ZIV_DECL (loop);

    /* compute the precision of intermediary variable */
    /* The optimal number of bits : see algorithms.tex */
    Nt = Ny + 3 + MPFR_INT_CEIL_LOG2 (Ny);

    /* initialise of intermediary variables */
    mpfr_init2 (t, Nt);
    mpfr_init2 (te, Nt);

    /* First computation of cosh */
    MPFR_ZIV_INIT (loop, Nt);
    for (;;)
      {
	/* Compute cosh */
        mpfr_clear_flags ();
	mpfr_exp (te, x, GMP_RNDD);         /* exp(x) */
        /* exp can overflow (but not underflow since x>0) */
        /* BUG/TODO/FIXME: exp can overflow but cosh may be representable! */
        if (MPFR_UNLIKELY (mpfr_overflow_p ()))
          {
            inexact = mpfr_overflow (y, rnd_mode, MPFR_SIGN_POS);
            MPFR_SAVE_EXPO_UPDATE_FLAGS (expo, MPFR_FLAGS_OVERFLOW);
            break;
          }
	mpfr_ui_div (t, 1, te, GMP_RNDU);   /* 1/exp(x) */
	mpfr_add (t, te, t, GMP_RNDU);      /* exp(x) + 1/exp(x)*/
	mpfr_div_2ui (t, t, 1, GMP_RNDN);   /* 1/2(exp(x) + 1/exp(x))*/

        /* Estimation of the error */
	err = Nt - 3;
	/* Check if we can round */
	if (MPFR_LIKELY (MPFR_CAN_ROUND (t, err, Ny, rnd_mode)))
          {
            inexact = mpfr_set (y, t, rnd_mode);
            break;
          }

	/* Actualisation of the precision */
	MPFR_ZIV_NEXT (loop, Nt);
        mpfr_set_prec (t, Nt);
        mpfr_set_prec (te, Nt);
      }
    MPFR_ZIV_FREE (loop);

    mpfr_clear (te);
    mpfr_clear (t);
  }

  MPFR_SAVE_EXPO_FREE (expo);
  inexact = mpfr_check_range (y, inexact, rnd_mode);
  return inexact;
}
