/* mpfr_set_f -- set a MPFR number from a GNU MPF number

Copyright (C) 1999 PolKA project, Inria Lorraine and Loria

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
#include "longlong.h"
#include "mpfr.h"

void 
#if __STDC__
mpfr_set_f(mpfr_ptr y, mpf_srcptr x, unsigned long prec, char rnd_mode)
#else
mpfr_set_f(y, x, prec, rnd_mode)
     mpfr_ptr y;
     mpf_srcptr x;
     unsigned long prec;
     char rnd_mode; 
#endif
{
  mp_limb_t *my, *mx; unsigned long cnt, sx, sy; 

  sx = SIZ (x); sy = ABSSIZE(y);
  my = MANT(y); mx = MANT(x); 

  count_leading_zeros(cnt, mx[sx - 1]);  

  if (sy < sx)
    {
      if (cnt) 
	{ 
	  mpn_lshift(my, mx + 1, sy, cnt); 
	  my [0] |= mx[0] >> (BITS_PER_MP_LIMB - cnt); 
	  EXP(y) = EXP(x)* BITS_PER_MP_LIMB - cnt; 
	}
      else { MPN_COPY(my, mx, sy); EXP(y) = EXP(x) * BITS_PER_MP_LIMB; }
    }
  else
    {
      if (cnt) 
	{ 
	  mpn_lshift(my, mx, sx, cnt); 
	  EXP(y) = EXP(x)* BITS_PER_MP_LIMB - cnt; 
	}
      else { MPN_COPY(my, mx, sy); EXP(y) = EXP(x) * BITS_PER_MP_LIMB; }
    }
  
  mpfr_round(y, rnd_mode, prec); 
}
