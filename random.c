/* mpfr_random -- generate a random floating-point number

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
#include "mpfr.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "urandom.h"
#include "srandom.h"

/* Computes a random mpfr in [0, 1[ with precision MPFR_PREC */

void
#if __STDC__
mpfr_random(mpfr_ptr x)
#else
mpfr_random(x)
     mpfr_ptr x; 
#endif    
{
  mp_limb_t *xp; unsigned long xn, i, cnt, prec=MPFR_PREC(x); 

  MPFR_CLEAR_FLAGS(x);
  xp = MPFR_MANT(x);
  xn = (prec-1)/BITS_PER_MP_LIMB + 1;

  for (i = 0; i < xn; i++)
    xp[i] = urandom();
  
  count_leading_zeros(cnt, xp[xn - 1]); 
  if (cnt) mpn_lshift(xp, xp, xn, cnt); 
  MPFR_EXP(x) = -cnt; 

  cnt = xn*BITS_PER_MP_LIMB - prec; 
  /* cnt is the number of non significant bits in the low limb */
  xp[0] &= ~((((mp_limb_t)1)<<cnt) - 1);
}

void
#if __STDC__
mpfr_srandom(unsigned long int seed)
#else
mpfr_srandom(seed)
     unsigned long int seed;
#endif
{
  srandom(seed); 
}
