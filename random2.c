/* mpf_random2 -- Generate a positive random mpf_t of specified size, with
   long runs of consecutive ones and zeros in the binary representation.
   Intended for testing of other MP routines.

Copyright (C) 1995, 1996 Free Software Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Library General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
License for more details.

You should have received a copy of the GNU Library General Public License
along with the GNU MP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

#if defined (__hpux) || defined (__alpha)  || defined (__svr4__) || defined (__SVR4)
/* HPUX lacks random().  DEC OSF/1 1.2 random() returns a double.  */
long mrand48 ();
static inline long
random ()
{
  return mrand48 ();
}
#else
long random ();
#endif

void
#if __STDC__
mpfr_random2 (mpfr_ptr x, mp_size_t size, mp_exp_t exp)
#else
mpfr_random2 (x, size, exp)
     mpfr_ptr x;
     mp_size_t size;
     mp_exp_t exp;
#endif
{
  mp_size_t xn; unsigned long cnt; mp_ptr xp = MANT(x); 
  mp_size_t prec = (PREC(x) - 1)/BITS_PER_MP_LIMB; 

  xn = ABS (size);
  if (xn != 0)
    {
      if (xn > prec + 1)
	xn = prec + 1;

      mpn_random2 (xp, xn);
    }

  if (exp != 0)
    exp = random () % (2 * exp) - exp;

  count_leading_zeros(cnt, xp[xn - 1]); 
  if (cnt) mpn_lshift(xp, xp, xn, cnt); 
  EXP(x) = exp-cnt; 
  cnt = xn*BITS_PER_MP_LIMB - prec; 
  /* cnt is the number of non significant bits in the low limb */
  xp[0] &= ~((((mp_limb_t)1)<<cnt) - 1);
}
