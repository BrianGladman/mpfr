/* mpfr_sqrt -- square root of a floating-point number

Copyright (C) 1999-2001 Free Software Foundation.
Contributed by the Spaces project.

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
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"

/* #define DEBUG */

int
#if __STDC__
mpfr_sqrt (mpfr_ptr r, mpfr_srcptr u, mp_rnd_t rnd_mode)
#else
mpfr_sqrt (r, u, rnd_mode)
     mpfr_ptr r;
     mpfr_srcptr u;
     mp_rnd_t rnd_mode;
#endif
{
  mp_ptr up, rp, tmp, remp;
  mp_size_t usize, rrsize;
  mp_size_t rsize;
  mp_size_t prec, err;
  mp_limb_t q_limb;
  long rw, nw, k; 
  int inexact = 0, t;
  unsigned long cc = 0; 
  char can_round = 0; 
  TMP_DECL(marker0);
  {
  TMP_DECL (marker);

  if (MPFR_IS_NAN(u)) {
    MPFR_SET_NAN(r);
    return 1; /* NaN is always inexact */
  }

  if (MPFR_SIGN(u) < 0) {
    if (MPFR_IS_INF(u) || MPFR_NOTZERO(u)) {
      MPFR_SET_NAN(r);
      return 1; /* NaN is always inexact */
    }
    else { /* sqrt(-0) = -0 */
      MPFR_SET_ZERO(r);
      if (MPFR_SIGN(r) > 0) MPFR_CHANGE_SIGN(r);
      return 0; /* zero is exact */
    }
  }

  MPFR_CLEAR_NAN(r);

  if (MPFR_SIGN(r) < 0) MPFR_CHANGE_SIGN(r);
  if (MPFR_IS_INF(u)) 
    { 
      MPFR_SET_INF(r);
      return 0; /* infinity is exact */
    }

  MPFR_CLEAR_INF(r);

  prec = MPFR_PREC(r);

  if (!MPFR_NOTZERO(u))
    {
      MPFR_EXP(r) = 0;
      rsize = (prec-1)/BITS_PER_MP_LIMB + 1;
      MPN_ZERO(MPFR_MANT(r), rsize);
      return 0; /* zero is exact */
    }

  up = MPFR_MANT(u);
  usize = (MPFR_PREC(u) - 1)/BITS_PER_MP_LIMB + 1; 

#ifdef DEBUG
      printf("Entering square root : "); 
      for(k = usize - 1; k >= 0; k--) { printf("%lu ", up[k]); }
      printf(".\n"); 
#endif

  /* Compare the mantissas */
  
  rsize = ((MPFR_PREC(r) + 2 + (MPFR_EXP(u) & 1))/BITS_PER_MP_LIMB + 1) << 1; 
  rrsize = (MPFR_PREC(r) + 2 + (MPFR_EXP(u) & 1))/BITS_PER_MP_LIMB + 1;
  /* One extra bit is needed in order to get the square root with enough
     precision ; take one extra bit for rrsize in order to solve more 
     easily the problem of rounding to nearest.
     Need to have 2*rrsize = rsize...
     Take one extra bit if the exponent of u is odd since we shall have
     to shift then.
  */

  TMP_MARK(marker0); 
  if (MPFR_EXP(u) & 1) /* Shift u one bit to the right */
    {
      if (MPFR_PREC(u) & (BITS_PER_MP_LIMB - 1))
	{
	  up = TMP_ALLOC(usize*BYTES_PER_MP_LIMB);
	  mpn_rshift(up, MPFR_MANT(u), usize, 1); 
	}
      else
	{
	  up = TMP_ALLOC((usize + 1)*BYTES_PER_MP_LIMB);
	  if (mpn_rshift(up + 1, MPFR_MANT(u), usize, 1))
	    up [0] = ((mp_limb_t) 1) << (BITS_PER_MP_LIMB - 1); 
	  else up[0] = 0; 
	  usize++; 
	}
    }

  MPFR_EXP(r) = ((MPFR_EXP(u) + (MPFR_EXP(u) & 1)) / 2) ;  
  
  do
    {
      TMP_MARK (marker);

      err = rsize*BITS_PER_MP_LIMB; 
      if (rsize < usize) { err--; }
      if (err > rrsize * BITS_PER_MP_LIMB) 
	{ err = rrsize * BITS_PER_MP_LIMB; }
      
      tmp = (mp_ptr) TMP_ALLOC (rsize * BYTES_PER_MP_LIMB);  
      rp = (mp_ptr) TMP_ALLOC (rrsize * BYTES_PER_MP_LIMB); 
      remp = (mp_ptr) TMP_ALLOC (rsize * BYTES_PER_MP_LIMB); 

      if (usize >= rsize) { 
	MPN_COPY (tmp, up + usize - rsize, rsize);
      }
      else { 
	MPN_COPY (tmp + rsize - usize, up, usize);
	MPN_ZERO (tmp, rsize - usize); 
      }

      /* Do the real job */
 
#ifdef DEBUG
      printf("Taking the sqrt of : "); 
      for(k = rsize - 1; k >= 0; k--) { 
	printf("+%lu*2^%lu",tmp[k],k*BITS_PER_MP_LIMB); }
      printf(".\n"); 
#endif

      q_limb = mpn_sqrtrem_new (rp, remp, tmp, rsize);

#ifdef DEBUG
      printf ("The result is : \n"); 
      printf ("sqrt : "); 
      for (k = rrsize - 1; k >= 0; k--)
	printf ("%lu ", rp[k]);
      printf ("(inexact = %lu)\n", q_limb); 
#endif
      
      can_round = mpfr_can_round_raw(rp, rrsize, 1, err, 
				     GMP_RNDZ, rnd_mode, MPFR_PREC(r));

      /* If we used all the limbs of both the dividend and the divisor, 
	 then we have the correct RNDZ rounding */

      if (!can_round && (rsize < 2*usize)) 
	{ 
#ifdef DEBUG
	  printf("Increasing the precision.\n"); 
#endif
	  TMP_FREE(marker); 
	}
    }
  while (!can_round && (rsize < 2*usize) 
	 && (rsize += 2) && (rrsize ++)); 
#ifdef DEBUG
  printf ("can_round = %d\n", can_round);
#endif

  /* This part may be deplaced upper to avoid a few mpfr_can_round_raw */
  /* when the square root is exact. It is however very unprobable that */
  /* it would improve the behaviour of the present code on average.    */

  if (!q_limb) /* the sqrtrem call was exact, possible exact square root */
    {
      /* if we have taken into account the whole of up */
      for (k = usize - rsize - 1; k >= 0; k ++)
	if (up[k]) break; 
      
      if (k < 0)
	goto fin; /* exact square root ==> inexact = 0 */
    }

  if (can_round) 
    {
      cc = mpfr_round_raw (rp, rp, err, 0, MPFR_PREC(r), rnd_mode, &inexact);
      if (!inexact) /* exact high part: inexact flag depends from remainder */
	inexact = -q_limb;
      rrsize = (MPFR_PREC(r) - 1)/BITS_PER_MP_LIMB + 1; 
    }
  else
    /* Use the return value of sqrtrem to decide of the rounding         */
    /* Note that at this point the sqrt has been computed                */
    /* EXACTLY.                                                          */
    switch (rnd_mode)
      {
      case GMP_RNDZ : 
      case GMP_RNDD :
	inexact = -1; /* result is truncated */
	break; 

      case GMP_RNDN : 
	/* Not in the situation ...0 111111 */
	rw = (MPFR_PREC(r) + 1) & (BITS_PER_MP_LIMB - 1);
	if (rw)
	  {
	    rw = BITS_PER_MP_LIMB - rw;
	    nw = 0;
	  }
	else
	  nw = 1; 
	if ((rp[nw] >> rw) & 1 &&                     /* Not 0111111111 */
	    (q_limb ||                                /* Nonzero remainder */
	    (rw ? (rp[nw] >> (rw + 1)) & 1 : 
	     (rp[nw] >> (BITS_PER_MP_LIMB - 1)) & 1))) /* or even rounding */ 
	  {
	    cc = mpn_add_1 (rp + nw, rp + nw, rrsize, ((mp_limb_t)1) << rw);
	    inexact = 1;
	  }
	else
	  inexact = -1;
	break;
 
      case GMP_RNDU:
	/* we should arrive here only when the result is inexact,
	   i.e. either q_limb > 0 (the remainder from mpn_sqrtrem is non-zero)
	   or up[0..usize-rsize-1] is non zero, thus we have to add one
	   ulp, and inexact = 1 */
	inexact = 1;
	t = MPFR_PREC(r) & (BITS_PER_MP_LIMB - 1); 
	rsize = (MPFR_PREC(r) - 1)/BITS_PER_MP_LIMB + 1;
	if (t) 
	    cc = mpn_add_1 (rp + rrsize - rsize, rp + rrsize - rsize, rsize, 1 << (BITS_PER_MP_LIMB - t));
	else
	    cc = mpn_add_1 (rp + rrsize - rsize, rp + rrsize - rsize, rsize, 1);
      }

  if (cc) {
    mpn_rshift(rp, rp, rrsize, 1);
    rp[rrsize-1] |= (mp_limb_t) 1 << (BITS_PER_MP_LIMB-1);
    MPFR_EXP(r)++; 
  }

 fin:
  rsize = rrsize; 
  rrsize = (MPFR_PREC(r) - 1)/BITS_PER_MP_LIMB + 1;  
  MPN_COPY(MPFR_MANT(r), rp + rsize - rrsize, rrsize); 

  if (MPFR_PREC(r) & (BITS_PER_MP_LIMB - 1))
    MPFR_MANT(r) [0] &= ~(((mp_limb_t)1 << (BITS_PER_MP_LIMB - 
				   (MPFR_PREC(r) & (BITS_PER_MP_LIMB - 1)))) - 1) ; 
  
  TMP_FREE (marker);
  }
  TMP_FREE(marker0);
  return inexact;
}
