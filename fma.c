/* mpfr_fma -- Floating multiply-add

Copyright (C) 2001 Free Software Foundation.

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
  int inexact;

        /* particular cases */
        if (MPFR_IS_NAN(x) || MPFR_IS_NAN(y) ||  MPFR_IS_NAN(z))
	  {
	    MPFR_SET_NAN(s); 
	    return 1; 
	  }

        /* cases Inf*0+z, 0*Inf+z, Inf-Inf */
        if ((MPFR_IS_INF(x) && MPFR_IS_ZERO(y)) ||
	    (MPFR_IS_INF(y) && MPFR_IS_ZERO(x)) ||
	    ((MPFR_IS_INF(x) || MPFR_IS_INF(y)) && MPFR_IS_INF(z) &&
	     ((MPFR_SIGN(x) * MPFR_SIGN(y)) != MPFR_SIGN(z))))
	  {
	    MPFR_SET_NAN(s);
	    return 1;
	  }

	MPFR_CLEAR_NAN(s);

        if (MPFR_IS_INF(x) || MPFR_IS_INF(y))
	  {
	    if (MPFR_IS_INF(z)) /* case Inf-Inf already checked above */
	      {
		MPFR_SET_INF(s);
		MPFR_SET_SAME_SIGN(s, z);
		return 0;
	      }
	    else /* z is finite */
	      {
		MPFR_SET_INF(s);
		if (MPFR_SIGN(s) != (MPFR_SIGN(x) * MPFR_SIGN(y)))
		  MPFR_CHANGE_SIGN(s);
		return 0;
	      }
	  }

	/* now x and y are finite */
	if (MPFR_IS_INF(z))
	  {
	    MPFR_SET_INF(s);
	    MPFR_SET_SAME_SIGN(s, z);
	    return 0;
	  }

        MPFR_CLEAR_INF(s);

        if(!MPFR_NOTZERO(x) || !MPFR_NOTZERO(y))
	  return mpfr_set (s, z, rnd_mode);
        if (!MPFR_NOTZERO(z))
	  return mpfr_mul (s, x, y, rnd_mode);
	
        /* General case */
        /* Detail of the compute */
        /* u <- x*y exact */
        /* s <- z+u */
	{
	  mpfr_t u;

	  /* if we take prec(u) >= prec(x) + prec(y), the product
	     u <- x*y is always exact */
	  mpfr_init2 (u, MPFR_PREC(x) + MPFR_PREC(y));
	  
	  mpfr_mul (u, x, y, GMP_RNDN); /* always exact */
	  inexact = mpfr_add (s, z, u, rnd_mode);
	  mpfr_clear(u);
	}
        return inexact;
}

