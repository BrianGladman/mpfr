/* mpfr_expm1 -- Compute exp(x)-1

Copyright 2001, 2002, 2003, 2004 Free Software Foundation.

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

#include "mpfr-impl.h"

 /* The computation of expm1 is done by
    expm1(x)=exp(x)-1
 */

int
mpfr_expm1 (mpfr_ptr y, mpfr_srcptr x , mp_rnd_t rnd_mode) 
{
  int inexact = 0;

  if (MPFR_UNLIKELY( MPFR_IS_SINGULAR(x) ))
    {
      if (MPFR_IS_NAN(x))
	{
	  MPFR_SET_NAN(y);
	  MPFR_RET_NAN;
	}
      /* check for inf or -inf (expm1(-inf)=-1) */
      else if (MPFR_IS_INF(x))
	{ 
	  if (MPFR_IS_POS(x))
	    {
	      MPFR_SET_INF(y);
	      MPFR_SET_POS(y);
	      return 0;
	    }
	  else
	    return mpfr_set_si(y, -1, rnd_mode);
	}
      else if(MPFR_IS_ZERO(x))
	{
	  MPFR_SET_ZERO(y);   /* expm1(+/- 0) = +/- 0 */
	  MPFR_SET_SAME_SIGN(y,x);
	  MPFR_RET(0);
	}
      else
	MPFR_ASSERTN(0);
    }
  /* Useless due to mpfr_set
     MPFR_CLEAR_FLAGS(y);*/

  /* General case */
  {
    /* Declaration of the intermediary variable */
    mpfr_t te,t;

    /* Declaration of the size variable */
    mp_prec_t Nx = MPFR_PREC(x);   /* Precision of input variable */
    mp_prec_t Ny = MPFR_PREC(y);   /* Precision of input variable */
    
    mp_prec_t Nt;   /* Precision of the intermediary variable */
    long int err;  /* Precision of error */
                
    /* compute the precision of intermediary variable */
    Nt=MAX(Nx,Ny);
    /* the optimal number of bits : see algorithms.ps */
    Nt=Nt+5+__gmpfr_ceil_log2(Nt);

    /* initialise of intermediary	variable */
    mpfr_init(t);             
    mpfr_init(te);             

    /* First computation of cosh */
    do
      {

        /* reactualisation of the precision */
        mpfr_set_prec(t,Nt);
        mpfr_set_prec(te,Nt);

        /* compute expm1 */
        mpfr_exp(te,x,GMP_RNDN);        /* exp(x)*/
        mpfr_sub_ui(t,te,1,GMP_RNDN);   /* exp(x)-1 */

        /* error estimate */
        /*err=Nt-(__gmpfr_ceil_log2(1+pow(2,MPFR_EXP(te)-MPFR_EXP(t))));*/
        err = Nt - (MAX (MPFR_GET_EXP (te) - MPFR_GET_EXP (t), 0) + 1);

        /* actualisation of the precision */
        Nt += 10;
      }
    while ((err < 0) || !mpfr_can_round (t, err, GMP_RNDN, GMP_RNDZ,
                                         Ny + (rnd_mode == GMP_RNDN)));
 
    inexact = mpfr_set (y, t, rnd_mode);

    mpfr_clear (t);
    mpfr_clear (te);
  }

  return inexact;
}
