/* mpfr_asin -- arc-sinus of a floating-point number

Copyright 2001, 2002, 2003, 2004 Free Software Foundation.

This file is part of the MPFR Library, and was contributed by Mathieu Dutour.

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

#include "mpfr-impl.h"

int
mpfr_asin (mpfr_ptr asin, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  mpfr_t xp;
  int compared, inexact;
  mp_prec_t prec;
  mp_exp_t xp_exp;
  MPFR_SAVE_EXPO_DECL (expo);

  /* Special cases */
  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    {
      if (MPFR_IS_NAN (x) || MPFR_IS_INF (x))
	{
	  MPFR_SET_NAN (asin);
	  MPFR_RET_NAN;
	}
      else /* x = 0 */
	{
          MPFR_ASSERTD (MPFR_IS_ZERO (x));
	  MPFR_SET_ZERO (asin);
	  MPFR_RET (0); /* exact result */
	}
    }

  /* Set x_p=|x| (x is a normal number) */
  mpfr_init2 (xp, MPFR_PREC (x));
  inexact = mpfr_abs (xp, x, GMP_RNDN); 
  MPFR_ASSERTD (inexact == 0);

  compared = mpfr_cmp_ui (xp, 1);

  if (MPFR_UNLIKELY (compared >= 0))
    {
      mpfr_clear (xp);
      if (compared > 0)                  /* asin(x) = NaN for |x| > 1 */
	{
	  MPFR_SET_NAN (asin);
	  MPFR_RET_NAN;
	}
      else                              /* x = 1 or x = -1 */
	{
	  if (MPFR_IS_POS (x)) /* asin(+1) = Pi/2 */
	    inexact = mpfr_const_pi (asin, rnd_mode);
	  else /* asin(-1) = -Pi/2 */
	    {
	      inexact = -mpfr_const_pi (asin, MPFR_INVERT_RND(rnd_mode));
	      MPFR_CHANGE_SIGN (asin);
	    }
	  mpfr_div_2ui (asin, asin, 1, rnd_mode); /* May underflow */
	  return inexact;
	}
    }

  MPFR_SAVE_EXPO_MARK (expo);

  /* Compute exponent of 1 - ABS(x) */
  mpfr_ui_sub (xp, 1, xp, GMP_RNDD);
  MPFR_ASSERTD (MPFR_GET_EXP (xp) <= 0);
  MPFR_ASSERTD (MPFR_GET_EXP (x) <= 0);
  xp_exp = 2 - MPFR_GET_EXP (xp); 

  /* Set up initial prec */
  prec = MPFR_PREC (asin) + 10 + xp_exp;

  /* If x ~ 2^-N, asin(x) ~ x + x^3/6
     If Prec < 2*N, we can't round since x^3/6 won't be count. */
  if (MPFR_PREC (asin) >= MPFR_PREC (x) 
      && prec <= -2*MPFR_GET_EXP (x) + 10)
    prec = -2*MPFR_GET_EXP (x) + 10;

  /* use asin(x) = atan(x/sqrt(1-x^2)) */

  for (;;)
    {
      mpfr_set_prec (xp, prec);
      mpfr_sqr (xp, x, GMP_RNDN);
      mpfr_ui_sub (xp, 1, xp, GMP_RNDN);
      mpfr_sqrt (xp, xp, GMP_RNDN);
      mpfr_div (xp, x, xp, GMP_RNDN);
      mpfr_atan (xp, xp, GMP_RNDN);
      if (mpfr_can_round (xp, prec - xp_exp, GMP_RNDN, GMP_RNDZ,
                          MPFR_PREC (asin) + (rnd_mode == GMP_RNDN)))
	break;
      prec += BITS_PER_MP_LIMB;
    }

  inexact = mpfr_set (asin, xp, rnd_mode);

  mpfr_clear (xp);

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (asin, inexact, rnd_mode);
}
