/* mpfr_asinh -- Inverse Hyperbolic Sinus of Unsigned Integer Number

Copyright (C) 1999 Free Software Foundation.

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

 /* The computation of acosh is done by

    asinh= ln(x+sqrt(x^2+1))
 */
int mpfr_asinh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));

int
#if __STDC__
mpfr_asinh (mpfr_ptr y, mpfr_srcptr x , mp_rnd_t rnd_mode) 
#else
mpfr_asinh (y, x, rnd_mode)
     mpfr_ptr y;
     mpfr_srcptr x;
     mp_rnd_t rnd_mode;
#endif
{
    
  /****** Declaration ******/

    /* Variable of Intermediary Calculation*/
    mpfr_t t,xt;       

    /* Variable of Intermediary Calculation*/
    mpfr_t ta,tb;

    int round;
    int boucle;
    int flag_neg;

    mp_prec_t Nx;   /* Precision of input variable */
    mp_prec_t Ny;   /* Precision of output variable */
    mp_prec_t Nt;   /* Precision of Intermediary Calculation variable */
    mp_prec_t err;  /* Precision of error */

    Nx=MPFR_PREC(x);
    mpfr_init2(xt,Nx);
    mpfr_set(xt,x,GMP_RNDN); 
    
    if (MPFR_IS_NAN(xt)) {  MPFR_SET_NAN(y); return 1; }
    MPFR_CLEAR_NAN(y);

    if (MPFR_IS_INF(xt)){ 
      MPFR_SET_INF(y);
      if(MPFR_SIGN(xt)>0){
	if (MPFR_SIGN(y) < 0) MPFR_CHANGE_SIGN(y);
      }
      else{
	if (MPFR_SIGN(y) > 0) MPFR_CHANGE_SIGN(y);
      }
      return 1;
    }

    MPFR_CLEAR_INF(y);

    flag_neg=0;

    if(MPFR_SIGN(xt)<0){
      MPFR_CHANGE_SIGN(xt);
      flag_neg=1;
    }


    if(!MPFR_NOTZERO(xt)){
      MPFR_SET_ZERO(y);   /* asinh(0) = 0 */
      if (MPFR_SIGN(y) < 0) MPFR_CHANGE_SIGN(y);
      return(0);
    }
    else{
      

        /* Initialisation of the Precision */
	Nx=MPFR_PREC(xt);
	Ny=MPFR_PREC(y);

	/* compute the size of intermediary variable */
	if(Ny>=Nx)
	  Nt=Ny+2*(BITS_PER_CHAR); 
	else
	  Nt=Nx+2*(BITS_PER_CHAR); 
	  
	  boucle=1;

	  /* initialise of intermediary	variable */
	  mpfr_init2(t,Nt);             
      	  mpfr_init2(ta,Nt);             
	  mpfr_init2(tb,Nt);             

	  while(boucle==1){


	    /* compute asinh */

	    mpfr_mul(ta,xt,xt,GMP_RNDN);     /* (x^2) */
	    mpfr_add_ui(ta,ta,1,GMP_RNDN);   /* (x^2+1) */
	    mpfr_sqrt(ta,ta,GMP_RNDN);       /* sqrt(x^2+1) */
	    mpfr_add(t,ta,xt,GMP_RNDN);      /* sqrt(x^2+1)+x */
	    mpfr_log(t,t,GMP_RNDN);          /* ln(sqrt(x^2+1)+x) */
	    
	    err=Nt-1-MAX(0,(-MPFR_EXP(t)));

	    round=mpfr_can_round(t,err,GMP_RNDN,rnd_mode,Ny);

	    if(round == 1){
	      if(flag_neg==1)MPFR_CHANGE_SIGN(t);
	      mpfr_set(y,t,rnd_mode);
	      boucle=0;
	    }
	    else{
	      Nt=Nt+10;
	      /* initialise of intermediary	variable */
	      mpfr_set_prec(t,Nt);             
	      mpfr_set_prec(ta,Nt);             
	      mpfr_set_prec(tb,Nt);             

	      boucle=1;
	    }
	    
	  }
    
	  mpfr_clear(t);
	  mpfr_clear(ta);
	  mpfr_clear(tb);
	  mpfr_clear(xt);
          return(1);
      
 
    }
}
