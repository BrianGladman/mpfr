/* mpfr_dim -- dim of x, y  

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

 /* The computation of z=dim(x,y)

    x-y if x >  y
    +0    if x <= y
 */

int
#if __STDC__
mpfr_dim (mpfr_ptr z, mpfr_srcptr x ,mpfr_srcptr y , mp_rnd_t rnd_mode) 
#else
mpfr_dim (z, x, y, rnd_mode)
     mpfr_ptr z;
     mpfr_srcptr x;
     mpfr_srcptr y;
     mp_rnd_t rnd_mode;
#endif
{
    
  /****** Declaration ******/

    /* Variable of Intermediary Calculation*/
    mpfr_t t;       

    mp_prec_t Nx;   /* Precision of input variable */
    mp_prec_t Nz;   /* Precision of input variable */
    mp_prec_t Ny;   /* Precision of output variable */
    mp_prec_t Nt;

    int inexact=0;

    if (MPFR_IS_NAN(x) || MPFR_IS_NAN(y) ) 
    {  
      MPFR_SET_NAN(z); 
      return 1; 
    }
    MPFR_CLEAR_NAN(z);

    /* Initialisation of the Precision */
    Nx=MPFR_PREC(x);
    Ny=MPFR_PREC(y);
    Nz=MPFR_PREC(z);

    /* compute the size of intermediary variable */
    Nt=MAX(Nx,MAX(Ny,Nz));
    
    /* initialise of intermediary variable */
    mpfr_init2(t,Nt);             

    if(mpfr_cmp(x,y) > 0)
      mpfr_sub(t,x,y,GMP_RNDN);
    else
    {
      MPFR_SET_ZERO(t);
      if(MPFR_SIGN(t) < 0) 
        CHANGE_SIGN(t);
    }

    inexact = mpfr_set(z,t,rnd_mode);
    
    mpfr_clear(t);
    return inexact;

}
