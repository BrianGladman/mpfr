/* mpfr_exp2 -- power of 2 function 2^y 

<<<<<<< exp2.c
Copyright (C) 2001 Free Software Foundation.
=======
Copyright (C) 1999-2001 Free Software Foundation.
>>>>>>> 1.19

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
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"


 /* The computation of y=2^z where z is a long int*/
int mpfr_exp2_si _PROTO((mpfr_ptr, long int n, mp_rnd_t));

int
#if __STDC__
mpfr_exp2_si (mpfr_ptr y, long int n, mp_rnd_t rnd_mode) 
#else
mpfr_exp2_si (y, n, rnd_mode)
     mpfr_ptr y;
     long int n;
     mp_rnd_t rnd_mode;
#endif
{    

  mpfr_set_ui(y,1,rnd_mode);
  MPFR_EXP(y)+=n;
  return 0;
}


 /* The computation of y=pow(2,z) is done by

    y=exp(z*log(2))=2^z
 */
int mpfr_exp2 _PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));

int
#if __STDC__
mpfr_exp2 (mpfr_ptr y, mpfr_srcptr x, mp_rnd_t rnd_mode) 
#else
mpfr_exp2 (y, x, rnd_mode)
     mpfr_ptr y;
     mpfr_srcptr x;
     mp_rnd_t rnd_mode;
#endif
{    

    if (MPFR_IS_NAN(x))
      {
        MPFR_SET_NAN(y); 
        return 1; 
      }

    MPFR_CLEAR_NAN(y);

    if (MPFR_IS_INF(x))
      {
        if (MPFR_SIGN(x) < 0)
          {
            MPFR_SET_ZERO(y);
            if (MPFR_SIGN(y) < 0)
              MPFR_CHANGE_SIGN(y);
            return 0;
          }
        else
          {
            MPFR_SET_INF(y);
            if(MPFR_SIGN(y) < 0) 
              MPFR_CHANGE_SIGN(y);
            return 0;
          }
      }

    /* 2^0 = 1 */
    if(mpfr_cmp_ui(x,0)==0)
      {
        return mpfr_set_ui(y,1,rnd_mode); 
      }

    /* General case */
    {
    /* Declaration of the intermediary variable */
      mpfr_t t, te;

      /* Declaration of the size variable */
      mp_prec_t Nx = MPFR_PREC(x);   /* Precision of input variable */
      mp_prec_t Ny = MPFR_PREC(y);   /* Precision of input variable */

      mp_prec_t Nt;   /* Precision of the intermediary variable */
      mp_prec_t err;  /* Precision of error */
                
      /* compute the precision of intermediary variable */
      Nt=MAX(Nx,Ny);
      /* the optimal number of bits : see algorithms.ps */
      Nt=Nt+5+_mpfr_ceil_log2(Nt);


      /* initialise of intermediary	variable */
      mpfr_init(t);             
      mpfr_init(te);             

      /* First computation of cosh */
      do {

	/* reactualisation of the precision */
	mpfr_set_prec(t,Nt);             
	mpfr_set_prec(te,Nt);             

	/* compute   exp(x*ln(2))*/
        mpfr_const_log2(t,GMP_RNDU);    /* ln(2) */
        mpfr_mul(te,x,t,GMP_RNDU);       /* x*ln(2) */
        mpfr_exp(t,te,GMP_RNDN);         /* exp(x*ln(2))*/

	/* estimation of the error -- see pow function in algorithms.ps*/
	err=Nt-_mpfr_ceil_log2(1+pow(2,MPFR_EXP(te)+1));

	/* actualisation of the precision */
	Nt += 10;

      } while (!mpfr_can_round(t,err,GMP_RNDN,rnd_mode,Ny));
 
      mpfr_set(y,t,rnd_mode);
      mpfr_clear(t);
      mpfr_clear(te);
    }
    return 1;

}
