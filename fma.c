/* mpfr_fma -- Floating multiply-add

Copyright (C) 1999 Free Software Foundation.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute
it and/or modify it under the terms of the GNU Library
General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your
option) any later version.

The MPFR Library is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Library General Public
License for more details.

You should have received a copy of the GNU Library
General Public License along with the MPFR Library; see
the file COPYING.LIB.  If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"


/* The computation of fma of x y and u is done by
    fma(s,x,y,z)= z + x*y = s
*/

int mpfr_fma _PROTO((mpfr_ptr, mpfr_srcptr,mpfr_srcptr,mpfr_srcptr,mp_rnd_t));

int
#if __STDC__
mpfr_fma (mpfr_ptr s, mpfr_srcptr x ,mpfr_srcptr y ,mpfr_srcptr z , mp_rnd_t rnd_mode) 
#else
mpfr_fma (s,x,y,z, rnd_mode)
     mpfr_ptr s;
     mpfr_srcptr x;
     mpfr_srcptr y;
     mpfr_srcptr z;
     mp_rnd_t rnd_mode;
#endif
{
        /* particular cases */
        if (MPFR_IS_NAN(x) || MPFR_IS_NAN(y) ||  MPFR_IS_NAN(z)){
                MPFR_SET_NAN(s); 
                return 1; 
        }
        
        if (MPFR_IS_INF(x) && MPFR_IS_ZERO(y)){ 
                MPFR_SET_NAN(s);
                return 1;
        }
        if (MPFR_IS_INF(y) && MPFR_IS_ZERO(x)){ 
                MPFR_SET_NAN(s);
                return 1;
        }
        if (MPFR_IS_INF(x) && MPFR_IS_INF(z)
                        && MPFR_SIGN(x)!=MPFR_SIGN(z)){
                MPFR_SET_NAN(s);
                return 1;
        }

        MPFR_CLEAR_NAN(s);
        MPFR_CLEAR_INF(s);

        if(!MPFR_NOTZERO(x) || !MPFR_NOTZERO(y)){
                mpfr_set(s,z,rnd_mode);
                return 1;
        }
        if (!MPFR_NOTZERO(z)) {
                mpfr_mul(s,x,y,rnd_mode);
                return 1;
        }

        /* General case */
        /* Detail of the compute */
        /* u <- x*y */
        /* t <- z+u */
        {
                /* Declaration of the intermediary variable */
                mpfr_t t, u;       
                int d;

                /* Declaration of the size variable */
                mp_prec_t Nx = MPFR_PREC(x);   /* Precision of input variable */
                mp_prec_t Ny = MPFR_PREC(y);   /* Precision of input variable */
                mp_prec_t Nz = MPFR_PREC(z);   /* Precision of input variable */
                mp_prec_t Ns = MPFR_PREC(s);   /* Precision of output variable */
                mp_prec_t Nt;   /* Precision of the intermediary variable */
                mp_prec_t err;  /* Precision of error */
                unsigned int first_pass=0; /* temporary precision */
                
                /* compute the precision of intermediary variable */
                Nt=MAX(MAX(Nx,Ny),Nz);

                /* the optimal number of bits is MPFR_EXP(u)-MPFR_EXP(v)+1 */
                /* but u and v are not yet compute, also we take in account */
                /* just one bit */
                Nt=Nt+1+_mpfr_ceil_log2(Nt);

                /* initialise the intermediary variables */
                mpfr_init(u);             
                mpfr_init(t);             
                
                /* First computation of fma */
                do {
                        /* reactualisation of the precision */
                        mpfr_set_prec(u,Nt);             
                        mpfr_set_prec(t,Nt);             

                        /* computations */
                        mpfr_mul(u,x,y,GMP_RNDN);         
                        mpfr_add(t,z,u,GMP_RNDN);        
                        
                        /*Nt=Nt+(d+1)+_mpfr_ceil_log2(Nt); */
                        d = MPFR_EXP(u)-MPFR_EXP(t);

                        /* estimation of the error */
                        err=Nt-(d+1);

                        /* actualisation of the precision */
                        Nt += (1-first_pass)*d + first_pass*10;
                        
                        first_pass=1;

                } while (mpfr_can_round(t,err,GMP_RNDN,rnd_mode,Ns)!=1);
                  
                mpfr_set(s,t,rnd_mode);
                mpfr_clear(t);
                mpfr_clear(u);
        }
        return 1;
}

