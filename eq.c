/* mpfr_eq -- Compare two floats up to a specified bit #.

Copied from mpf_eq.

Copyright (C) 1993, 1995, 1996 Free Software Foundation, Inc.

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
#include "mpfr.h"

int
#if __STDC__
mpfr_eq (mpfr_srcptr u, mpfr_srcptr v, unsigned long int n_bits)
#else
mpfr_eq (u, v, n_bits)
     mpfr_srcptr u;
     mpfr_srcptr v;
     unsigned long int n_bits;
#endif
{
  mp_srcptr up, vp;
  mp_size_t usize, vsize, size, i;
  mp_exp_t uexp, vexp;
  int usign, k;

  uexp = u->_mp_exp;
  vexp = v->_mp_exp;

  usize = (PREC(u)-1)/BITS_PER_MP_LIMB + 1;
  vsize = (PREC(v)-1)/BITS_PER_MP_LIMB + 1;

  /* 1. Are the signs different?  */
  if (MPFR_SIGN(u) ==  MPFR_SIGN(v))
    {
      /* U and V are both non-negative or both negative.  */
      if (!NOTZERO(u))
	return !NOTZERO(v); 
      if (!NOTZERO(v))
	return !NOTZERO(u);

      /* Fall out.  */
    }
  else
    {
      /* Either U or V is negative, but not both.  */
      return 0;
    }

  /* U and V have the same sign and are both non-zero.  */

  usign = MPFR_SIGN(u); 

  /* 2. Are the exponents different?  */
  if (uexp > vexp)
    return 0;			/* ??? handle (uexp = vexp + 1)   */
  if (vexp > uexp)
    return 0;			/* ??? handle (vexp = uexp + 1)   */

  usize = ABS (usize);
  vsize = ABS (vsize);

  up = u->_mp_d;
  vp = v->_mp_d;

  if (usize > vsize)
    {
      if (vsize * BITS_PER_MP_LIMB < n_bits)
	{
	  k = usize - vsize - 1; 
	  while (k >= 0 && !up[k]) --k; 
	  if (k >= 0) 
	    return 0;		/* surely too different */
	}
      size = vsize;
    }
  else if (vsize > usize)
    {
      if (usize * BITS_PER_MP_LIMB < n_bits)
	{
	  k = vsize - usize - 1; 
	  while (k >= 0 && !vp[k]) --k; 
	  if (k >= 0) 
	    return 0;		/* surely too different */
	}
      size = usize;
    }
  else
    {
      size = usize;
    }

  if (size > (n_bits + BITS_PER_MP_LIMB - 1) / BITS_PER_MP_LIMB)
    size = (n_bits + BITS_PER_MP_LIMB - 1) / BITS_PER_MP_LIMB;

  up += usize - size;
  vp += vsize - size;

  for (i = size - 1; i > 0; i--)
    {
      if (up[i] != vp[i])
	return 0;
    }

  if (n_bits & (BITS_PER_MP_LIMB - 1))
    return (up[i] >> (BITS_PER_MP_LIMB - (n_bits & (BITS_PER_MP_LIMB - 1))) == 
	    vp[i] >> (BITS_PER_MP_LIMB - (n_bits & (BITS_PER_MP_LIMB - 1)))); 
  else
    return (up[i] == vp[i]); 
}
