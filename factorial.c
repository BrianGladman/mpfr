/* mpfr_facrorial -- Factorial of Unsigned Integer Number

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

 /* The computation of n! is done by

    n!=prod^{n}_{i=1}i
 */

int
#if __STDC__
mpfr_factorial (mpfr_ptr y, unsigned long int x , mp_rnd_t rnd_mode) 
#else
mpfr_factorial (y, x, rnd_mode)
     mpfr_ptr y;
     unsigned long x;
     mp_rnd_t rnd_mode;
#endif
{

  /****** Declaration ******/

    mpfr_t t;       /* Variable of Intermediary Calculation*/
    int i;
    int round;
    int boucle;

    mp_prec_t Nx;   /* Precision of input variable */
    mp_prec_t Ny;   /* Precision of output variable */
    mp_prec_t Nt;   /* Precision of Intermediary Calculation variable */
    mp_prec_t err;  /* Precision of error */

  /***** test x = 0  ******/
  	  
    if(x==0){
      mpfr_set_ui(y,1,GMP_RNDN); /* 0! = 1 */
      return(0);
    }
    else{
 

        /* Initialisation of the Precision */
	Nx=sizeof(unsigned long int)*(BITS_PER_CHAR);
	Ny=MPFR_PREC(y);
	 
	Nt=Ny+2*(int)_mpfr_ceil_log2((double)x)+10; /*compute the size of intermediary variable */

	
	  boucle=1;
	  mpfr_init2(t,Nt);/* initialise of intermediary variable */

	  while(boucle==1){
	   
	    mpfr_set_ui(t, 1, GMP_RNDZ);

	    for(i=2;i<=x;i++)              /* compute factorial */
	      mpfr_mul_ui(t,t,i,GMP_RNDZ);
	    
	    err=Nt-1-(int)_mpfr_ceil_log2((double)Nt);


	    round=mpfr_can_round(t,err,GMP_RNDZ,rnd_mode,Ny);

	    if(round == 1){
	      mpfr_set(y,t,rnd_mode);
	      boucle=0;
	    }
	    else{
	      Nt=Nt+10; 
	      /*initialise of intermediary variable */
	      mpfr_set_prec(t,Nt);
	      boucle=1;
	    }
	    
	  }
   
	  mpfr_clear(t);
          return(1);
      
 
    }
}
