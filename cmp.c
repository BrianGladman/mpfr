/* mpfr_cmp -- compare two floating-point numbers

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
#include "longlong.h"
#include "mpfr.h"
#include "mpfr-impl.h"

/* returns 0 iff b = c
            a positive value iff b > c
            a negative value iff b < c

More precisely, in case b and c are of same sign, the absolute value 
of the result is one plus the absolute difference between the exponents 
of b and c, i.e. one plus the number of bits shifts to align b and c
(this value is useful in mpfr_sub).

*/

/* #define DEBUG */

/* compares b and sign(s)*c */
int 
#if __STDC__
mpfr_cmp3 (mpfr_srcptr b, mpfr_srcptr c, long int s)
#else
mpfr_cmp3 (b, c, s)
     mpfr_srcptr b;
     mpfr_srcptr c; 
     long int s;
#endif
{
   long int diff_exp;
   unsigned long bn, cn;
   mp_limb_t *bp, *cp;

   if (MPFR_IS_NAN(b) || MPFR_IS_NAN(c)) return 1;

   if (MPFR_IS_INF(b)) {
     if (MPFR_IS_INF(c) && (MPFR_SIGN(b) * s * MPFR_SIGN(c) > 0))
       return 0;
     else 
       return MPFR_SIGN(b);
   }

   if (!MPFR_NOTZERO(b)) {
     if (!MPFR_NOTZERO(c)) return 0; else return -(s*MPFR_SIGN(c));
   }
   else if (!MPFR_NOTZERO(c)) return MPFR_SIGN(b);

   s = s * MPFR_SIGN(b) * MPFR_SIGN(c);
   if (s<0) return(MPFR_SIGN(b));

   /* now signs are equal */
   diff_exp = MPFR_EXP(b)-MPFR_EXP(c);
   s = (MPFR_SIGN(b) > 0) ? 1 : -1;

   if (diff_exp>0) return(s*(1+diff_exp));
   else if (diff_exp<0) return(s*(-1+diff_exp));
   /* both signs and exponents are equal */

   bn = (MPFR_PREC(b)-1)/BITS_PER_MP_LIMB+1;
   cn = (MPFR_PREC(c)-1)/BITS_PER_MP_LIMB+1;
   bp = MPFR_MANT(b); cp = MPFR_MANT(c);

   while (bn && cn) {
     if (bp[--bn] != cp[--cn])
       return((bp[bn]>cp[cn]) ? s : -s);
   }

   if (bn) { while (bn) if (bp[--bn]) return(s); }
   else if (cn) while (cn) if (cp[--cn]) return(-s);

   return 0;
}

