/* mpfr_cbrt -- power function x^(1/3)

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
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"

 /* The computation of y=x^(1/3) is done by

    Case exp-log y=e^((1/3)*log(x))
    Case Newton y / y_{k+1}=(1/3)[(x/(y_k^2))+2y_k]
 */

int
#if __STDC__
mpfr_cbrt (mpfr_ptr y, mpfr_srcptr x , mp_rnd_t rnd_mode) 
#else
mpfr_cbrt (y, x, rnd_mode)
     mpfr_ptr y;
     mpfr_srcptr x;
     mp_rnd_t rnd_mode;
#endif
{
    
  /****** Declaration ******/

    /* Variable of Intermediary Calculation*/
    mpfr_t t1,t2,t;       

    int round;
    int boucle;
    long int exp_t;
    ldiv_t epsilon;    
    int tau=2;
    int ktau=0;
    int i;

    mp_prec_t Nx;   /* Precision of input variable */
    mp_prec_t Ny;   /* Precision of output variable */
    mp_prec_t Nt;   /* Precision of Intermediary Calculation variable */
    mp_prec_t Ntemp;   /* Precision of Intermediary Calculation variable */
    mp_prec_t err;  /* Precision of error */

    /* Gestion des NaN */
    if (MPFR_IS_NAN(x)) {  MPFR_SET_NAN(y); return 1; }
    MPFR_CLEAR_NAN(y);

    /* Gestion des infinies*/
    if (MPFR_IS_INF(x)){ 
      MPFR_SET_INF(y);
      if(MPFR_SIGN(x) > 0) {
	if (MPFR_SIGN(y) < 0) MPFR_CHANGE_SIGN(y);}
      else{
	if (MPFR_SIGN(y) < 0) MPFR_CHANGE_SIGN(y);}

      return 1;
    }
    MPFR_CLEAR_INF(y);

    /*Gestion du cas 0*/
    if(!MPFR_NOTZERO(x)){
      MPFR_SET_ZERO(y);   /* cbrt(+/- 0) = +/- 0 */

      if(MPFR_SIGN(x) > 0){
	if (MPFR_SIGN(y) < 0) MPFR_CHANGE_SIGN(y);
      }
      else{
	if (MPFR_SIGN(y) > 0) MPFR_CHANGE_SIGN(y);
      }
      return 0;
    }


    /* Initialisation of the Precision */
    Nx=MPFR_PREC(x);
    Ny=MPFR_PREC(y);

    /* compute the size of intermediary variable */
    if(Ny>=Nx)
      Nt=Ny; 
    else
      Nt=Nx;

    /* Calcul du nombre d'iteration necessaire pour newton*/
    /* t0=2, t{k+1}=2.t{k}-1 k / tk>n */ 
    
    while(tau<=Nt){
      tau=2*tau-1;
      ktau++;
    }


    /* Calcul de la taille des variable temporaire */

    Ntemp=0;
    for(i=0;i<ktau;i++){
      Ntemp=10*(Ntemp)+17;
      epsilon=ldiv(Ntemp,3);
      Ntemp=epsilon.quot+1;
    }

    Nt=Nt+(int)_mpfr_ceil_log2((double)Nt)+(int)_mpfr_ceil_log2((double)Ntemp);

    mpfr_init2(t1,Nt);
    mpfr_init2(t2,Nt);
    mpfr_init2(t,Nt);

    mpfr_set(t,x,GMP_RNDN);


    /* normalisation de la valeur de t */
    /* tel que t= (m/2^r)  x 2^(3e') avec e=3e'-r exposant et m mantisse de t*/

    exp_t=(int)MPFR_EXP(t);
    epsilon=ldiv(exp_t,3);
    mpfr_div_2exp(t,t,MPFR_EXP(t),GMP_RNDN);
    mpfr_div_2exp(t,t,(3-epsilon.rem),GMP_RNDN);


    /*Gestion des negatifs*/
    if(MPFR_SIGN(x)<0) MPFR_CHANGE_SIGN(t);
     
    boucle=1;



    while(boucle==1){

      /* compute cbrt */
      /*mpfr_log(t,x,GMP_RNDN);*/         /* ln(x) */
      /*mpfr_div_ui(t,t,3,GMP_RNDN);*/    /* ln(x)/3 */
      /*mpfr_exp(t,t,GMP_RNDN);*/         /* exp(ln(x)/3)*/

      mpfr_set_d(t1,0.75,GMP_RNDN);

      for(i=0;i<ktau;i++){

	mpfr_mul_2exp(t2,t1,1,GMP_RNDN);  /*2x*/
	mpfr_mul(t1,t1,t1,GMP_RNDN);      /*x^2*/
        mpfr_div(t1,t,t1,GMP_RNDN);       /*N/x^2*/
	mpfr_add(t1,t1,t2,GMP_RNDN);      /*2x+N/x^2*/
	mpfr_div_ui(t1,t1,3,GMP_RNDN);    /*(1/3)[2x+N/x^2]*/

      }

      
      err=Nt-1-(int)_mpfr_ceil_log2((double)Nt);

      round=mpfr_can_round(t1,err,GMP_RNDN,rnd_mode,Ny);

      
      if(round == 1){
	/*Gestion des negatifs*/
	if(MPFR_SIGN(x)<0) MPFR_CHANGE_SIGN(t1);
	mpfr_mul_2exp(t1,t1,(epsilon.quot+1),GMP_RNDN);
	mpfr_set(y,t1,rnd_mode);
	boucle=0;
      }
      else{
	Nt=Nt+10; 
	/* re-initialise of intermediary	variable */
	mpfr_set_prec(t1,Nt);             
	mpfr_set_prec(t2,Nt);             
	boucle=1;
      }

   }

   mpfr_clear(t1);
   mpfr_clear(t2);
   mpfr_clear(t);
   return(1);

   
    
}
