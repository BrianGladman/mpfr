/* mpfr_div -- divide two floating-point numbers

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
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"
#include "mpfr-impl.h"

// #define DEBUG

#ifdef DEBUG
void 
affiche_mp(mp_ptr z, mp_size_t length)
{
  int k; 

  if (length == 1) { printf("[%lu]\n", *z); return; }

  printf("[%lu, ", z[length - 1]); 
  for(k = length - 2; k >= 1; k--) { printf("%lu, ", z[k]); }
  printf("%lu]\n", z[0]); 
}
#endif

void
#if __STDC__
mpfr_div (mpfr_ptr r, mpfr_srcptr u, mpfr_srcptr v, mp_rnd_t rnd_mode)
#else
mpfr_div (r, u, v, rnd_mode)
     mpfr_ptr r;
     mpfr_srcptr u;
     mpfr_srcptr v;
     mp_rnd_t rnd_mode;
#endif
{
  mp_srcptr up, vp;
  mp_ptr rp, tp, tp0, tmp, tmp2;
  mp_size_t usize, vsize, drsize, dsize;
  mp_size_t oldrsize, rsize, sign_quotient, err;
  mp_limb_t q_limb;
  mp_exp_t rexp;
  long k, t; 
  unsigned long cc = 0, rw, nw; 
  int can_round = 0, sh; 
  TMP_DECL (marker);


  /**************************************************************************
   *                                                                        *
   *              This part of the code deals with special cases            *
   *                                                                        *
   **************************************************************************/

  if (MPFR_IS_NAN(u) || MPFR_IS_NAN(v)) { MPFR_SET_NAN(r); return; }

  MPFR_CLEAR_NAN(r);

  if (MPFR_IS_INF(u)) 
    { 
      if (MPFR_IS_INF(v)) 
	MPFR_SET_NAN(r);
      else
	{ 
	  MPFR_SET_INF(r); 
	  if (MPFR_SIGN(r) != MPFR_SIGN(u) * MPFR_SIGN(v)) 
	    MPFR_CHANGE_SIGN(r);
	}
      return;
    }
  else 
    if (MPFR_IS_INF(v)) 
      {
	MPFR_CLEAR_INF(r);
	MPFR_SET_ZERO(r); 
	if (MPFR_SIGN(r) != MPFR_SIGN(u) * MPFR_SIGN(v)) 
	  MPFR_CHANGE_SIGN(r);
	return; 
      }

  MPFR_CLEAR_INF(r); /* clear Inf flag */

  if (!MPFR_NOTZERO(v))
    {
      if (!MPFR_NOTZERO(u)) 
	{ MPFR_SET_NAN(r); return; }
      else
	{
	  MPFR_SET_INF(r); 
	  if (MPFR_SIGN(r) != MPFR_SIGN(v) * MPFR_SIGN(u)) 
	    MPFR_CHANGE_SIGN(r); 
	  return;
	}
    }
  
  if (!MPFR_NOTZERO(u)) { MPFR_SET_ZERO(r); return; }

  /**************************************************************************
   *                                                                        *
   *              End of the part concerning special values.                *
   *                                                                        *
   **************************************************************************/


#ifdef DEBUG
  printf("u = "); mpfr_out_str(stdout, 2, 0, u, GMP_RNDN); printf("\n"); 
  printf("v = "); mpfr_out_str(stdout, 2, 0, v, GMP_RNDN); printf("\n"); 
#endif

  up = MPFR_MANT(u);
  vp = MPFR_MANT(v);

  TMP_MARK (marker);

  usize = (MPFR_PREC(u) - 1)/BITS_PER_MP_LIMB + 1; 
  vsize = (MPFR_PREC(v) - 1)/BITS_PER_MP_LIMB + 1; 

#ifdef DEBUG
  printf("Entering division : "); 
  affiche_mp(up, usize); 
  affiche_mp(vp, vsize); 
#endif
      
  /**************************************************************************
   *                                                                        *
   *   First try to use only part of u, v. If this is not sufficient,       *
   *   use the full u and v, to avoid long computations eg. in the case     *
   *   u = v.                                                               *
   *                                                                        *
   **************************************************************************/

  dsize = (MPFR_PREC(r) + 3)/BITS_PER_MP_LIMB + 1; 
  drsize = MPFR_PREC(r)/BITS_PER_MP_LIMB + 1 + dsize; 

  /* Compute effective dividend */
  if (vsize < dsize) { dsize = vsize; }
  tmp = (mp_ptr) vp + vsize - dsize; 

  /* Compute effective divisor. One extra bit allowed for RNDN. */
  
  tp0 = (mp_ptr) TMP_ALLOC(drsize * BYTES_PER_MP_LIMB); 
  if (usize >= drsize) 
    MPN_COPY (tp0, up + usize - drsize, drsize); 
  else {
    MPN_COPY (tp0 + drsize - usize, up, usize); 
    MPN_ZERO (tp0, drsize - usize); 
  }

  /* Allocate limbs for quotient. */
  rp = (mp_ptr) TMP_ALLOC ((drsize - dsize + 1) * BYTES_PER_MP_LIMB);
      
#ifdef DEBUG
  printf("Dividing : "); 
  affiche_mp(tp0, drsize); 
  printf(" by "); 
  affiche_mp(tmp, dsize); 
  printf(".\n"); 
#endif
      
#if (__GNU_MP_VERSION < 3)
  q_limb = mpn_divrem (rp, 0, tp0, drsize, tmp, dsize);
  tp = tp0;
#else 
  tmp2 = (mp_ptr) TMP_ALLOC (dsize * BYTES_PER_MP_LIMB);
  mpn_tdiv_qr(rp, tmp2, 0, tp0, drsize, tmp, dsize);
  q_limb = rp[drsize - dsize];
  tp = tmp2;
#endif
  
  /* Estimate number of correct bits. */
  err = (drsize - dsize) * BITS_PER_MP_LIMB; 
  if (drsize < usize) err --; 
  if (dsize < vsize) err -= 2; 

#ifdef DEBUG
  printf("The result is : \n"); 
  printf("Quotient : "); 

  affiche_mp(rp, drsize - dsize + 1); 
 
  printf("Remainder : "); 
  affiche_mp(tp, dsize); 

  printf("Number of correct bits = %lu\n", err); 
#endif
      
  /* Compute sign and exponent. Correction might occur just below. */
  sign_quotient = ((MPFR_SIGN(u) * MPFR_SIGN(v) > 0) ? 1 : -1); 
  rexp = MPFR_EXP(u) - MPFR_EXP(v);
  rsize = drsize - dsize;
  
  /* We want to check if rounding is possible. Have to mimic normalization */
  if (q_limb) { sh = -1; } 
  else { count_leading_zeros(sh, rp[rsize - 1]); } 
  
  can_round = mpfr_can_round_raw(rp, rsize + 1, sign_quotient, err + sh + 
				 BITS_PER_MP_LIMB, GMP_RNDN, rnd_mode, 
				 MPFR_PREC(r) + sh + BITS_PER_MP_LIMB); 

  if (!can_round && (drsize < usize || dsize < vsize)) 
    {
      mp_size_t ulrsize; 
      int b = 0; 
      mp_ptr rp2, ulorem, ulorem2; 

  /**************************************************************************
   *                                                                        *
   *   The attempt to use only part of u and v failed. We first compute a   *
   *   correcting term, then perform the full division.                     *
   *   Put u = uhi + ulo, v = vhi + vlo. We have uhi = vhi * rp + tp,       *
   *   thus u - v * rp = tp + ulo - rp*vlo, that we shall divide by v.      *
   *                                                                        *
   **************************************************************************/
#ifdef DEBUG
      printf("Using the full u and v.\n");
#endif
      
      if (usize > drsize)
	  ulrsize = usize - drsize + dsize; 
	  // store the low part of u and the remainder. 
      else ulrsize = dsize; // just store the remainder. 

      if (vsize > dsize) ulrsize++; 

      ulorem = TMP_ALLOC(ulrsize * BYTES_PER_MP_LIMB);
      ulorem2 = TMP_ALLOC(ulrsize * BYTES_PER_MP_LIMB);
      ulorem[ulrsize - 1] = ulorem2 [ulrsize - 1] = 0; 

      if (dsize < vsize)
	{
	  mp_size_t lg = vsize + rsize - dsize + 1; 

	  if (rsize > vsize - dsize)
	    mpn_mul(ulorem + ulrsize - lg, rp, rsize, vp, vsize - dsize);
	  else
	    mpn_mul(ulorem + ulrsize - lg, vp, vsize - dsize, rp, rsize);

	  MPN_ZERO(ulorem, ulrsize - lg); 
	}
      else MPN_ZERO(ulorem, ulrsize); 

#ifdef DEBUG
      printf("vlo * q: "); 
      affiche_mp(ulorem, ulrsize); 
#endif 
     
      MPN_COPY(ulorem2 + ulrsize - 1 - dsize, tp, dsize);      
      if (drsize < usize) 
	MPN_COPY(ulorem2, up, usize - drsize); 

#ifdef DEBUG
      printf("(b = %d) ulo + r: ", b); 
      affiche_mp(ulorem2, ulrsize); 
#endif
      
      if (mpn_cmp(ulorem2, ulorem, ulrsize) > 0)
	mpn_sub_n(ulorem, ulorem2, ulorem, ulrsize); 
      else 
	{ 
	  ulorem2[ulrsize - 1] = 
	    mpn_add_n(ulorem2 + ulrsize - vsize - 1, 
		      ulorem2 + ulrsize - vsize - 1, vp, vsize); 

	  if (mpn_cmp(ulorem2, ulorem, ulrsize) < 0)
	    {
	      b = 2; 
	      ulorem2[ulrsize - 1] += 
		mpn_add_n(ulorem2 + ulrsize - vsize - 1, 
			  ulorem2 + ulrsize - vsize - 1, vp, vsize); 
	    }
	  else b = 1; 

	  mpn_sub_n(ulorem, ulorem2, ulorem, ulrsize); 
	}
	  
#ifdef DEBUG
  printf("(b = %d) ulo + r - vlo * q: ", b); 
  affiche_mp(ulorem, ulrsize); 
#endif
     
      rp2 = (mp_ptr) TMP_ALLOC((ulrsize - vsize + 1) * BYTES_PER_MP_LIMB); 

#if (__GNU_MP_VERSION < 3)
      q_limb = mpn_divrem (rp2, 0, ulorem, ulrsize, vp, vsize);
      tp = ulorem; 
#else 
      tmp2 = (mp_ptr) TMP_ALLOC (vsize * BYTES_PER_MP_LIMB);
      mpn_tdiv_qr(rp2, tmp2, 0, ulorem, ulrsize, vp, vsize);
      tp = tmp2; 
#endif

      if (!b)
	mpn_add_1(rp, rp, rsize, rp2[ulrsize - vsize]);
      else
	mpn_sub_1(rp, rp, rsize, b); 

      dsize = vsize; /* The length of the remainder, in the end of the code */
    }

  /**************************************************************************
   *                                                                        *
   *                       Final stuff (rounding and so.)                   *
   *                                                                        *
   **************************************************************************/

  if (q_limb)
    {
      mpn_rshift(rp, rp, rsize, 1);
      rp[rsize - 1] |= (mp_limb_t)1 << (BITS_PER_MP_LIMB - 1); rexp ++; 
    }
  else
    if (sh) { mpn_lshift(rp, rp, rsize, sh); rexp -= sh; }
  
  oldrsize = rsize;
  rsize = (MPFR_PREC(r) - 1)/BITS_PER_MP_LIMB + 1;

  if (can_round) /* Lazy case. */
    {
      cc = mpfr_round_raw(rp, rp, err, (sign_quotient == -1 ? 1 : 0),
			  MPFR_PREC(r), rnd_mode, NULL);
    }
  else {
    rp += oldrsize-rsize;
    if ((rnd_mode == GMP_RNDD && sign_quotient == -1) 
	|| (rnd_mode == GMP_RNDU && sign_quotient == 1)
	|| (rnd_mode == GMP_RNDN))
      {
	/* We cannot round, so that the last bits of the quotient
	   have to be zero; just look if the remainder is nonzero */

	k = dsize - 1; 
	while (k >= 0) { if (tp[k]) break; k--; }
	if (k >= 0) 
	  {
	    t = MPFR_PREC(r) & (BITS_PER_MP_LIMB - 1); 
	    if (t)
	      cc = mpn_add_1(rp, rp, rsize, 
			     (mp_limb_t)1 << (BITS_PER_MP_LIMB - t)); 
	    else
	      cc = mpn_add_1(rp, rp, rsize, 1); 
	  }
	else
	  {
	    if (rnd_mode == GMP_RNDN) /* even rounding */
	      {
		rw = (MPFR_PREC(r) + 1) & (BITS_PER_MP_LIMB - 1);
		if (rw) { rw = BITS_PER_MP_LIMB - rw; nw = 0; } else nw = 1; 
		if ((rw ? (rp[nw] >> (rw + 1)) & 1 : 
		     (rp[nw] >> (BITS_PER_MP_LIMB - 1)) & 1))
		  {
		    cc = mpn_add_1(rp + nw, rp + nw, rsize, 
				   ((mp_limb_t)1) << rw); 
		  }
	      }
	/* cas 0111111 */
	  }
      }
  }

  if (sign_quotient * MPFR_SIGN(r) < 0) { MPFR_CHANGE_SIGN(r); } 
  MPFR_EXP(r) = rexp;

  if (cc) {
    mpn_rshift(rp, rp, rsize, 1);
    rp[rsize-1] |= (mp_limb_t) 1 << (BITS_PER_MP_LIMB-1);
    MPFR_EXP(r)++; 
  }
    
  rw = rsize * BITS_PER_MP_LIMB - MPFR_PREC(r);
  MPN_COPY(MPFR_MANT(r), rp, rsize); 
  MPFR_MANT(r)[0] &= ~(((mp_limb_t)1 << rw) - 1);
#ifdef DEBUG
  printf("r = "); mpfr_out_str(stdout, 2, 0, r, GMP_RNDN); printf("\n"); 
#endif
  TMP_FREE (marker);
}
