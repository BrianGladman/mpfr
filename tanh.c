/* mpfr_tanh -- hyperbolic tangent

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
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

int
mpfr_tanh (mpfr_ptr y, mpfr_srcptr xt , mp_rnd_t rnd_mode) 
{
  /****** Declaration ******/
  mpfr_t x;
  int inexact;
  MPFR_SAVE_EXPO_DECL (expo);

  /* Special value checking */
  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (xt)))
    {
      if (MPFR_IS_NAN (xt)) 
	{
	  MPFR_SET_NAN (y); 
	  MPFR_RET_NAN;
	}
      else if (MPFR_IS_INF (xt))
	{
	  /* tanh(inf) = 1 && tanh(-inf) = -1 */
	  return mpfr_set_si (y, MPFR_INT_SIGN (xt), rnd_mode);
	}
      else /* tanh (0) = 0 and xt is zero */
	{
	  MPFR_ASSERTD (MPFR_IS_ZERO(xt));
	  MPFR_SET_ZERO (y);
	  MPFR_SET_SAME_SIGN (y, xt);
	  MPFR_RET (0);
	}
    }
  
  MPFR_SAVE_EXPO_MARK (expo);
  MPFR_TMP_INIT_ABS (x, xt);
  
  /* General case */
  {
    /* Declaration of the intermediary variable */
    mpfr_t t, te;
    mp_exp_t d;
    
    /* Declaration of the size variable */
    mp_prec_t Nx = MPFR_PREC(x);   /* Precision of input variable */
    mp_prec_t Ny = MPFR_PREC(y);   /* Precision of output variable */
    mp_prec_t Nt;                  /* Precision of intermediary variables */
    long int err;                  /* Precision of error */
    
    /* Compute the precision of intermediary variable */
    Nt = MAX (Nx, Ny);
    
    /* The optimal number of bits: see algorithms.ps */
    Nt = Nt + MPFR_INT_CEIL_LOG2 (Nt) + 4;
    
    /* initialise of intermediary variable */
    mpfr_init2 (t, Nt); 
    mpfr_init2 (te, Nt);
    
    if (MPFR_GET_EXP (x) > 10)
      for (;;) {
	/* tanh(x)=1-2/(exp(2x)+1)  */
	mpfr_mul_2ui (t, x, 1, GMP_RNDN);   /* 2x:        err = 0*/
	mpfr_exp (t, t, GMP_RNDN);          /* exp(2x):   err <= ulp(t) */
	mpfr_add_ui (t, t, 1, GMP_RNDD);    /* exp(2x)+1: err <= 2*ulp(t) */
	mpfr_ui_div (t, 1, t, GMP_RNDN);    /* 1/(exp(2x)+1): err <= 8*ulp(t)*/
	mpfr_mul_2ui (t, t, 1, GMP_RNDN);   /* 2/(exp(2x)+1): err <= 8*ulp(t)*/
	d = MPFR_GET_EXP (t);
	mpfr_ui_sub (t, 1, t, GMP_RNDZ);    /*1-2/(exp(2x)+1) */
	
	/* Calculation of the error*/
	/* err (t) <= (1+8*2^(d-EXP(t)))*ulp(t) */
	d = d - MPFR_GET_EXP (t);
	err = Nt - MAX (d + 4, 1);
	
	if (mpfr_can_round (t, err, GMP_RNDZ, GMP_RNDZ,
			    Ny + (rnd_mode == GMP_RNDN)))
	  break;
	
	/* Actualisation of the precision */
	Nt += BITS_PER_MP_LIMB;
	mpfr_set_prec (t, Nt);
	mpfr_set_prec (te, Nt);
      }
    else
      for (;;) {
	/* tanh = (exp(2x)-1)/(exp(2x)+1) */
	mpfr_mul_2ui (te, x, 1, GMP_RNDN);  /* 2x */
	mpfr_exp (te, te, GMP_RNDN);        /* exp(2x) */
	d = MPFR_GET_EXP (te);              /* For Error calculation */
	mpfr_add_ui (t, te, 1, GMP_RNDD);   /* exp(2x) + 1*/
	mpfr_sub_ui (te, te, 1, GMP_RNDU);  /* exp(2x) - 1*/
	mpfr_div (t, te, t, GMP_RNDN);      /* (exp(2x)-1)/(exp(2x)+1)*/
	
	/* Calculation of the error*/
	d = d - MPFR_GET_EXP (t);
	err = Nt - (MAX(d + 1, 3) + 1);
	
	if (mpfr_can_round (t, err, GMP_RNDN, GMP_RNDZ,
			    Ny + (rnd_mode == GMP_RNDN)))
	  break;
	
	/* Actualisation of the precision */
	Nt += BITS_PER_MP_LIMB;
	mpfr_set_prec (t, Nt);
	mpfr_set_prec (te, Nt);
      }
    
    inexact = mpfr_set4 (y, t, rnd_mode, MPFR_SIGN (xt));
    mpfr_clear (te);
    mpfr_clear (t);
  }
  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (y, inexact, rnd_mode);
}

