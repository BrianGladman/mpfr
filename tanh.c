/* mpfr_tanh -- hyperbolic tangent

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

    tanh= [e^(x)^2-1]/+[e^(x)^2+1]
 */
int mpfr_tanh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));

int
#if __STDC__
mpfr_tanh (mpfr_ptr y, mpfr_srcptr x , mp_rnd_t rnd_mode) 
#else
mpfr_tanh (y, x, rnd_mode)
     mpfr_ptr y;
     mpfr_srcptr x;
     mp_rnd_t rnd_mode;
#endif
{
    
  /****** Declaration ******/

    /* Variable of Intermediary Calculation*/
    mpfr_t t;       

    /* Variable of Intermediary Calculation*/
    mpfr_t te,ta,tb;

    int round;
    int boucle = 1, inexact = 0;

    mp_prec_t Nx;   /* Precision of input variable */
    mp_prec_t Ny;   /* Precision of output variable */
    mp_prec_t Nt;   /* Precision of Intermediary Calculation variable */
    mp_prec_t err;  /* Precision of error */
 

    if (MPFR_IS_NAN(x)) {  MPFR_SET_NAN(y); return 1; }
    MPFR_CLEAR_NAN(y);

    if (MPFR_IS_INF(x)){

      if(MPFR_SIGN(x) > 0)
	mpfr_set_si(y,1,GMP_RNDN); /* tanh(inf) = 1 */
      else
	mpfr_set_si(y,-1,GMP_RNDN); /* tanh(inf) = -1 */

      return 0;
    }

    if(!MPFR_NOTZERO(x)){               /* tanh(0) = 0 */
      	MPFR_CLEAR_INF(y);
	MPFR_SET_ZERO(y);
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
	  
	  boucle=1;

	    /* initialise of intermediary	variable */
	    mpfr_init2(t,Nt);             
	    mpfr_init2(te,Nt);                     
	    mpfr_init2(ta,Nt);             
	    mpfr_init2(tb,Nt);             

	  while (boucle)
	    {
	    /* compute tanh */
	    mpfr_mul_2exp(te,x,1,GMP_RNDN); /* 2x */
	    mpfr_exp(te,te,GMP_RNDN);       /* exp(2x) */
	    mpfr_add_ui(ta,te,1,GMP_RNDN);  /* exp(2x) + 1*/
	    mpfr_sub_ui(tb,te,1,GMP_RNDN);  /* exp(2x) - 1*/
	    mpfr_div(t,tb,ta,GMP_RNDN);     /* (exp(2x)-1)/(exp(2x)+1)*/

	    err=Nt-1-((MPFR_EXP(te)-MPFR_EXP(tb)));

	    round=mpfr_can_round(t,err,GMP_RNDN,rnd_mode,Ny);

	    if(round)
	      {
		inexact = mpfr_set (y, t, rnd_mode);
		boucle=0;
	      }
	    else
	      {
		Nt=Nt+10;
	      /* re-initialise of intermediary variable */
		mpfr_set_prec(t,Nt);             
		mpfr_set_prec(te,Nt);                     
		mpfr_set_prec(ta,Nt);             
		mpfr_set_prec(tb,Nt);             
	      }
	    
	    }
    
	  mpfr_clear(t);
	  mpfr_clear(te);
	  mpfr_clear(ta);
	  mpfr_clear(tb);
          return inexact;
}
