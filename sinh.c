/* mpfr_sinh -- hyperbolic sine

Copyright (C) 2001 Free Software Foundation.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Library General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
License for more details.

You should have received a copy of the GNU Library General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"

 /* The computation of cosh is done by

    cosh= 1/2[e^(x)-e^(-x)]
 */
int mpfr_sinh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));

int
#if __STDC__
mpfr_sinh (mpfr_ptr y, mpfr_srcptr x , mp_rnd_t rnd_mode) 
#else
mpfr_sinh (y, x, rnd_mode)
     mpfr_ptr y;
     mpfr_srcptr x;
     mp_rnd_t rnd_mode;
#endif
{
    
  /****** Declaration ******/

    /* Variable of Intermediary Calculation*/
    mpfr_t t;       

    /* Variable of Intermediary Calculation*/
    mpfr_t te,ti;

    int round;
    int boucle = 1, inexact = 0;

    mp_prec_t Nx;   /* Precision of input variable */
    mp_prec_t Ny;   /* Precision of output variable */
    mp_prec_t Nt;   /* Precision of Intermediary Calculation variable */
    mp_prec_t err;  /* Precision of error */

    
    if (MPFR_IS_NAN(x)) {  MPFR_SET_NAN(y); return 1; }
    MPFR_CLEAR_NAN(y);

    if (MPFR_IS_INF(x)){ 

      MPFR_SET_INF(y);

      if(MPFR_SIGN(x) > 0){
	if (MPFR_SIGN(y) < 0) MPFR_CHANGE_SIGN(y);
      }
      else{
	if (MPFR_SIGN(y) > 0) MPFR_CHANGE_SIGN(y);
      }
      return 0;
    }
    MPFR_CLEAR_INF(y);
  
    if(!MPFR_NOTZERO(x)){
      MPFR_SET_ZERO(y);   /* sinh(0) = 0 */
      if (MPFR_SIGN(y) < 0) MPFR_CHANGE_SIGN(y);
      return 0;
    }
 
        /* Initialisation of the Precision */
	Nx=MPFR_PREC(x);
	Ny=MPFR_PREC(y);

	/* compute the size of intermediary variable */
	if(Ny>=Nx)
	  Nt=Ny+2*(BITS_PER_CHAR); 
	else
	  Nt=Nx+2*(BITS_PER_CHAR); 
	  
	    /* initialise of intermediary	variable */
	    mpfr_init2(t,Nt);             
	    mpfr_init2(te,Nt);             
	    mpfr_init2(ti,Nt);             

	  while (boucle)
	    {

	    /* compute sinh */
	    mpfr_exp(te,x,GMP_RNDN);         /* exp(x) */
	    mpfr_ui_div(ti,1,te,GMP_RNDN);   /* 1/exp(x) */
	    mpfr_sub(t,te,ti,GMP_RNDN);      /* exp(x) - 1/exp(x)*/
	    mpfr_div_2exp(t,t,1,GMP_RNDN);   /* 1/2(exp(x) - 1/exp(x))*/

	    err=Nt-1-((MPFR_EXP(te)-MPFR_EXP(t)));

	    round=mpfr_can_round(t,err,GMP_RNDN,rnd_mode,Ny);

	    if(round == 1){
	      inexact = mpfr_set (y, t, rnd_mode);
	      boucle=0;
	    }
	    else{
	      Nt=Nt+10;
	      /* re-initialise of intermediary	variable */
	      mpfr_set_prec(t,Nt);             
	      mpfr_set_prec(te,Nt);             
	      mpfr_set_prec(ti,Nt);             
	    }
	    
	  }
    
	  mpfr_clear(t);
	  mpfr_clear(ti);
	  mpfr_clear(te);

          return inexact;
}
