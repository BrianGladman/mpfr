/* mpfr_atanh -- Inverse Hyperbolic Tangente of Unsigned Integer Number

Copyright (C) 1999-2001 Free Software Foundation.

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

#include <limits.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"

 /* The computation of acosh is done by

    atanh= 1/2*ln(x+1)-1/2*ln(1-x)
 */

int mpfr_atanh _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));

int
#if __STDC__
mpfr_atanh (mpfr_ptr y, mpfr_srcptr x , mp_rnd_t rnd_mode) 
#else
mpfr_atanh (y, x, rnd_mode)
     mpfr_ptr y;
     mpfr_srcptr x;
     mp_rnd_t rnd_mode;
#endif
{
    
  /****** Declaration ******/

    /* Variable of Intermediary Calculation*/
    mpfr_t t,xt;       

    /* Variable of Intermediary Calculation*/
    mpfr_t te,ti;

    int round;
    int comp;
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

    comp = mpfr_cmp_ui(xt,0);
    if(comp==0){

      MPFR_SET_ZERO(y);   /* atanh(0) = 0 */
      return(0);

    }
 

    flag_neg=0;

    if(MPFR_SIGN(xt)<0){
      MPFR_CHANGE_SIGN(xt);
      flag_neg=1;
    }

    comp=mpfr_cmp_ui(xt,1);

    if(comp >= 0){

	/*fprintf(stderr,"MPFR atanh function is 
	 only defined for x=[-1,+1]");
	exit(-1);*/

	/* An other srtategy if output is not define for input return NaN*/

	MPFR_SET_NAN(y);return(-1);
    }
    else{
      

        /* Initialization of the precision */
	Nx=MPFR_PREC(xt);
	Ny=MPFR_PREC(y);

	/* compute the size of temporary variable */
	if(Ny>=Nx)
	  Nt=Ny+2*CHAR_BIT+10;
	else
	  Nt=Nx+2*CHAR_BIT+10;
	  
	  boucle=1;

	    /* initialization of a temporary variable */
	    mpfr_init2(t,Nt);             
	    mpfr_init2(ti,Nt);  
	    mpfr_init2(te,Nt);             	    

	  while(boucle){

	    /* compute atanh */


	    /* Good algorithme near 1 but bad near 0 */
	    /*	    
	    mpfr_ui_sub(te,1,xt,GMP_RNDN);   e=1-x
	    mpfr_ui_div(te,2,te,GMP_RNDN);    2/e
	    mpfr_sub_ui(te,te,1,GMP_RNDN);    (2/e)-1
	    mpfr_log(te,te,GMP_RNDN);         ln((2/e)-1)
	    mpfr_div_2exp(t,te,1,GMP_RNDN);   (1/2)*ln((2/e)-1)
	    */


	    mpfr_ui_sub(te,1,xt,GMP_RNDN);   /* (1-xt)*/
	    mpfr_add_ui(ti,xt,1,GMP_RNDN);   /* (xt+1)*/
	    mpfr_div(te,ti,te,GMP_RNDN);     /* (1+xt)/(1-xt)*/
	    mpfr_log(te,te,GMP_RNDN);        /* ln((1+xt)/(1-xt))*/
	    mpfr_div_2exp(t,te,1,GMP_RNDN);  /* (1/2)*ln((1+xt)/(1-xt))*/
	    
	    err = Nt - MAX(0, -MPFR_EXP(t));

	    round=mpfr_can_round(t,err,GMP_RNDN,rnd_mode,Ny);
	    
	    if(round == 1){
	      if(flag_neg==1) MPFR_CHANGE_SIGN(t);
	      mpfr_set(y,t,rnd_mode);
	      boucle=0;
	    }
	    else{
	      Nt=Nt+10;
	      /* initialization of a temporary variable */
	      mpfr_set_prec(t,Nt);
	      mpfr_set_prec(te,Nt);             
	      mpfr_set_prec(ti,Nt);             
	      boucle=1;
	    }
	    
	  }
    
	  mpfr_clear(t);
	  mpfr_clear(te);
	  mpfr_clear(ti);
	  mpfr_clear(xt);
          return(1);
      
 
    }
}
