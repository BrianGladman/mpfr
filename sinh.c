/* mpfr_sinh -- hyperbolic sine

Copyright 2001, 2002, 2003, 2004, 2005 Free Software Foundation, Inc.

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

 /* The computation of sinh is done by
    sinh(x) = 1/2 [e^(x)-e^(-x)]          */

int
mpfr_sinh (mpfr_ptr y, mpfr_srcptr xt, mp_rnd_t rnd_mode)
{  
  mpfr_t x;
  int inexact;

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
	  MPFR_SET_SAME_SIGN(y, xt);
	  MPFR_RET(0);
	}
      else /* xt is zero */
	{
	  MPFR_ASSERTD(MPFR_IS_ZERO(xt));
	  MPFR_SET_ZERO(y);   /* sinh(0) = 0 */
	  MPFR_SET_SAME_SIGN(y, xt);
	  MPFR_RET(0);
	}
    }

  MPFR_TMP_INIT_ABS (x, xt);

  {
    mpfr_t t, ti;
    mp_exp_t d;    
    mp_prec_t Nt;    /* Precision of the intermediary variable */
    long int err;    /* Precision of error */
    MPFR_ZIV_DECL (loop);
    int overflow_p = mpfr_overflow_p ();

    /* compute the precision of intermediary variable */
    Nt = MAX (MPFR_PREC (x), MPFR_PREC (y));
    /* the optimal number of bits : see algorithms.ps */
    Nt = Nt + MPFR_INT_CEIL_LOG2 (Nt) + 4;
    /* If x is near 0, exp(x) - 1/exp(x) = 2*x+x^3/3+O(x^5) */
    if (MPFR_GET_EXP (x) < 0)
      Nt -= 2*MPFR_GET_EXP (x);

    /* initialise of intermediary	variable */
    mpfr_init2 (t, Nt);
    mpfr_init2 (ti, Nt);

    /* First computation of sinh */
    MPFR_ZIV_INIT (loop, Nt);
    for (;;) {
      /* compute sinh */
      mpfr_clear_overflow ();
      mpfr_exp (t, x, GMP_RNDD);        /* exp(x) */
      /* exp(x) can overflow or underflow or return ~1 ! */
      d = MPFR_GET_EXP (t);
      if (MPFR_UNLIKELY (mpfr_overflow_p ())) {
	MPFR_SET_INF (t);
	break;
      }
      mpfr_ui_div (ti, 1, t, GMP_RNDU); /* 1/exp(x) */
      mpfr_sub (t, t, ti, GMP_RNDN);    /* exp(x) - 1/exp(x) */
      mpfr_div_2ui (t, t, 1, GMP_RNDN);  /* 1/2(exp(x) - 1/exp(x)) */
      
      /* it may be that t is zero (in fact, it can only occur when te=1,
	 and thus ti=1 too) */
      err = 0;
      if (!MPFR_IS_ZERO (t))
	{
	  /* calculation of the error */
	  d = d - MPFR_GET_EXP (t) + 2;	  
	  /* error estimate */
	  /* err = Nt-(__gmpfr_ceil_log2(1+pow(2,d)));*/
	  err = Nt - (MAX (d, 0) + 1);
	  
	  if (mpfr_can_round (t, err, GMP_RNDN, GMP_RNDZ,
			      MPFR_PREC (y) + (rnd_mode == GMP_RNDN)))
	    break;
	}
      /* actualisation of the precision */
      Nt += err;
      MPFR_ZIV_NEXT (loop, Nt);
      mpfr_set_prec (t, Nt);
      mpfr_set_prec (ti, Nt);
    }
    MPFR_ZIV_FREE (loop);
    inexact = mpfr_set4 (y, t, rnd_mode, MPFR_SIGN (xt));
    if (overflow_p != 0)
      __gmpfr_flags |= MPFR_FLAGS_OVERFLOW;

    mpfr_clear (t);
    mpfr_clear (ti);
  }

  return inexact;
}
