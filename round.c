/* mpfr_round_raw2, mpfr_round_raw, mpfr_round, mpfr_can_round, 
   mpfr_can_round_raw -- various rounding functions

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
#include "mpfr.h"
#include "mpfr-impl.h"

#if (BITS_PER_MP_LIMB & (BITS_PER_MP_LIMB - 1))
#error "BITS_PER_MP_LIMB must be a power of 2"
#endif

/* returns 0 if round(sign*xp[0..xn-1], prec, rnd) = 
   round(sign*xp[0..xn-1], prec, GMP_RNDZ), 1 otherwise,
   where sign=1 if neg=0, sign=-1 otherwise.

   Does *not* modify anything.
*/
int 
#if __STDC__
mpfr_round_raw2 (mp_limb_t *xp, mp_prec_t xn, 
		int neg, mp_rnd_t rnd, mp_prec_t prec)
#else
mpfr_round_raw2 (xp, xn, neg, rnd, prec)
     mp_limb_t *xp; 
     mp_prec_t xn; 
     int neg; 
     mp_rnd_t rnd;
     mp_prec_t prec; 
#endif
{
  mp_prec_t nw; long wd; char rw; short l; mp_limb_t mask;

  nw = prec / BITS_PER_MP_LIMB; rw = prec & (BITS_PER_MP_LIMB - 1); 
  if (rw) nw++; 
  if (rnd==GMP_RNDZ || xn<nw || (rnd==GMP_RNDU && neg)
      || (rnd==GMP_RNDD && neg==0)) return 0;

  if (rw) 
    mask = ~((((mp_limb_t)1)<<(BITS_PER_MP_LIMB - rw)) - 1);
  else mask = ~((mp_limb_t)0); 

  switch (rnd)
    {
    case GMP_RNDU:
    case GMP_RNDD:
      if (xp[xn - nw] & ~mask) return 1;
      for (l = nw + 1;l <= xn; l++)
	if (xp[xn - l]) break;
      return (l <= xn);

    case GMP_RNDN:
    /* First check if we are just halfway between two representable numbers */
      wd = xn - nw;
      if (!rw)
	{
	  if (!wd) /* all bits are significative */ return 0; 
	  wd--;
	  if (xp[wd] == ((mp_limb_t)1 << (BITS_PER_MP_LIMB - 1)))
	    {
	      do wd--; while (wd > 0 && !xp[wd]);
	      if (!wd) { return 1; } else return xp[xn - nw] & 1;
	    }

	  return xp[wd]>>(BITS_PER_MP_LIMB - 1);
	}
      else
      if (rw + 1 < BITS_PER_MP_LIMB)
	{
	  if ((xp[wd] & (~mask)) == (((mp_limb_t)1) << (BITS_PER_MP_LIMB - rw - 1)))
	      do { wd--; } while (wd >= 0 && !xp[wd]);
	  else return ((xp[wd]>>(BITS_PER_MP_LIMB - rw - 1)) & 1);
	  
	  /* first limb was in the middle, and others down to wd+1 were 0 */
	  if (wd>=0) return 1;
	  else
	      return ((xp[xn - nw] & mask) >> (BITS_PER_MP_LIMB - rw)) & 1;
	}
      else
	/* Modified PZ, 27 May 1999:
	   rw, i.e. the number of bits stored in xp[xn-nw], is 
	   BITS_PER_MP_LIMB-1, i.e. there is exactly one non significant bit. 
	   We are just halfway iff xp[wd] has its low significant bit 
	   set and all limbs xp[0]...xp[wd-1] are zero */
	{
	  if (xp[wd] & 1)
	      do wd--; while (wd >= 0 && !xp[wd]);
	  return ((wd<0) ? xp[xn-nw]>>1 : xp[xn-nw]) & 1;
	}
    default: return 0;
    }
}

/* puts in y the value of xp (with precision xprec and sign 1 if negative=0,
   -1 otherwise) rounded to precision yprec and direction rnd_mode 
   Supposes x is not zero nor NaN nor +/- Infinity (i.e. *xp != 0).
*/
int
#if __STDC__
mpfr_round_raw (mp_limb_t *yp, mp_limb_t *xp, mp_prec_t xprec, int neg,
                mp_prec_t yprec, mp_rnd_t rnd_mode, int *inexp)
#else
mpfr_round_raw (yp, xp, xprec, neg, yprec, rnd_mode, inexp)
     mp_limb_t *yp;
     mp_limb_t *xp;
     mp_prec_t xprec;
     int neg;
     mp_prec_t yprec;
     mp_rnd_t rnd_mode;
     int *inexp;
#endif
{
  mp_size_t xsize, nw;
  mp_limb_t himask, lomask;
  int rw, carry = 0;

  xsize = (xprec-1)/BITS_PER_MP_LIMB + 1;

  nw = yprec / BITS_PER_MP_LIMB;
  rw = yprec & (BITS_PER_MP_LIMB - 1);
  if (rw)
  {
    nw++;
    lomask = ((((mp_limb_t)1)<<(BITS_PER_MP_LIMB - rw)) - (mp_limb_t)1);
    himask = ~lomask;
  }
  else
  {
    lomask = -1;
    himask = -1;
  }
  MPFR_ASSERTN(nw >= 1);

  if (xprec <= yprec)
  { /* No rounding is necessary. */
    /* if yp=xp, maybe an overlap: MPN_COPY_DECR is ok when src <= dst */
    MPFR_ASSERTN(nw >= xsize);
    MPN_COPY_DECR(yp + (nw - xsize), xp, xsize);
    MPN_ZERO(yp, nw - xsize); /* PZ 27 May 99 */
    if (inexp) *inexp = 0;
  }
  else
  {
    mp_limb_t rb, sb;

    if ((rnd_mode == GMP_RNDU && neg) ||
        (rnd_mode == GMP_RNDD && !neg))
      rnd_mode = GMP_RNDZ;

    if (inexp || rnd_mode != GMP_RNDZ)
    {
      mp_size_t k;
      mp_limb_t rbmask;

      k = xsize - nw;
      if (!rw) k--;
      MPFR_ASSERTN(k >= 0);
      sb = xp[k] & lomask;  /* First non-significant bits */
      if (rnd_mode == GMP_RNDN)
      {
        rbmask = ((mp_limb_t)1) << (BITS_PER_MP_LIMB - rw - 1);
        rb = sb & rbmask;
        sb &= ~rbmask;
      }
      while (sb == 0 && k > 0)
        sb = xp[--k];
      if (rnd_mode == GMP_RNDN && (rb != 0 || sb != 0))
      {
        sb = sb == 0 ? xp[xsize - nw] & (himask ^ (himask << 1)) : rb;
        if (inexp)
          *inexp = ((neg != 0) ^ (sb != 0)) ? 1 : -1;
      }
      else if (inexp)
        *inexp = sb == 0 ? 0
          : (((neg != 0) ^ (rnd_mode != GMP_RNDZ)) ? 1 : -1);
    }
    else
      sb = 0;

    if (sb != 0 && rnd_mode != GMP_RNDZ)
      carry = mpn_add_1(yp, xp + xsize - nw, nw,
                        rw ? ((mp_limb_t)1) << (BITS_PER_MP_LIMB - rw) : 1);
    else
      MPN_COPY_INCR(yp, xp + xsize - nw, nw);

    yp[0] &= himask;
  }

  return carry;
}

void
#if __STDC__
mpfr_round (mpfr_ptr x, mp_rnd_t rnd_mode, mp_prec_t prec)
#else
mpfr_round (x, rnd_mode, prec)
     mpfr_ptr x; 
     mp_rnd_t rnd_mode; 
     mp_prec_t prec; 
#endif
{
  mp_limb_t *tmp;
  int carry, neg;
  mp_prec_t nw;
  TMP_DECL(marker);

  if (MPFR_IS_NAN(x) || MPFR_IS_INF(x)) return; 
  nw = 1 + (prec - 1) / BITS_PER_MP_LIMB; /* needed allocated limbs */
  neg = MPFR_SIGN(x) < 0;

  /* check if x has enough allocated space for the mantissa */
  if (nw > MPFR_ABSSIZE(x)) {
    MPFR_MANT(x) = (mp_ptr) (*__gmp_reallocate_func) 
      (MPFR_MANT(x), MPFR_ABSSIZE(x)*BYTES_PER_MP_LIMB, nw * BYTES_PER_MP_LIMB);
    MPFR_SIZE(x) = nw; /* new number of allocated limbs */
    if (neg) MPFR_CHANGE_SIGN(x);
  }

  TMP_MARK(marker); 
  tmp = TMP_ALLOC (nw * BYTES_PER_MP_LIMB);
  carry = mpfr_round_raw(tmp, MPFR_MANT(x), MPFR_PREC(x), neg, prec, rnd_mode,
                         NULL);

  if (carry)
    {      
      mpn_rshift (tmp, tmp, nw, 1); 
      tmp [nw-1] |= (((mp_limb_t)1) << (BITS_PER_MP_LIMB - 1)); 
      MPFR_EXP(x)++; 
    }

  MPFR_PREC(x) = prec; 
  MPN_COPY(MPFR_MANT(x), tmp, nw); 
  TMP_FREE(marker);
}

/* hypotheses : BITS_PER_MP_LIMB est une puissance de 2 
                dans un test, la premiere partie du && est evaluee d'abord */


/* assuming b is an approximation of x in direction rnd1 
   with error at most 2^(MPFR_EXP(b)-err), returns 1 if one is 
   able to round exactly x to precision prec with direction rnd2,
   and 0 otherwise.

   Side effects: none.
*/

int 
#if __STDC__
mpfr_can_round (mpfr_ptr b, mp_prec_t err, mp_rnd_t rnd1, 
		mp_rnd_t rnd2, mp_prec_t prec)
#else
mpfr_can_round (b, err, rnd1, rnd2, prec) 
     mpfr_ptr b; 
     mp_prec_t err;
     mp_rnd_t rnd1;
     mp_rnd_t rnd2; 
     mp_prec_t prec;
#endif
{
  return mpfr_can_round_raw (MPFR_MANT(b),
			     (MPFR_PREC(b) - 1)/BITS_PER_MP_LIMB + 1, 
			     MPFR_SIGN(b), err, rnd1, rnd2, prec); 
}

int
#if __STDC__
mpfr_can_round_raw (mp_limb_t *bp, mp_prec_t bn, int neg, mp_prec_t err,
		    mp_rnd_t rnd1, mp_rnd_t rnd2, mp_prec_t prec)
#else
mpfr_can_round_raw (bp, bn, neg, err, rnd1, rnd2, prec)
     mp_limb_t *bp;
     mp_prec_t bn;
     int neg; 
     mp_prec_t err;
     mp_rnd_t rnd1;
     mp_rnd_t rnd2; 
     mp_prec_t prec; 
#endif
{
  int k, k1, l, l1, tn;
  mp_limb_t cc, cc2, *tmp;
  TMP_DECL(marker); 

  if (err<=prec) return 0;
  neg = (neg > 0 ? 0 : 1);

  /* if the error is smaller than ulp(b), then anyway it will propagate
     up to ulp(b) */
  if (err > bn * BITS_PER_MP_LIMB)
    err = bn * BITS_PER_MP_LIMB;

  /* warning: if k = m*BITS_PER_MP_LIMB, consider limb m-1 and not m */
  k = (err-1)/BITS_PER_MP_LIMB;
  l = err % BITS_PER_MP_LIMB; if (l) l = BITS_PER_MP_LIMB-l;
  /* the error corresponds to bit l in limb k, the most significant limb
   being limb 0 */
  k1 = (prec-1)/BITS_PER_MP_LIMB;
  l1 = prec%BITS_PER_MP_LIMB; if (l1) l1 = BITS_PER_MP_LIMB-l1;

  /* the last significant bit is bit l1 in limb k1 */

  /* don't need to consider the k1 most significant limbs */
  k -= k1; bn -= k1; prec -= k1*BITS_PER_MP_LIMB; k1=0;

  if (rnd1==GMP_RNDU) { if (neg) rnd1=GMP_RNDZ; }
  if (rnd1==GMP_RNDD) { if (neg) rnd1=GMP_RNDU; else rnd1=GMP_RNDZ; }

  /* in the sequel, RNDU = towards infinity, RNDZ = towards zero */

  TMP_MARK(marker);
  tn = bn;
  k++; /* since we work with k+1 everywhere */

  switch (rnd1) {
    
  case GMP_RNDZ: /* b <= x <= b+2^(MPFR_EXP(b)-err) */
    tmp = TMP_ALLOC(tn*BYTES_PER_MP_LIMB); 
    cc = (bp[bn-1]>>l1) & 1;
    cc ^= mpfr_round_raw2(bp, bn, neg, rnd2, prec);

    /* now round b+2^(MPFR_EXP(b)-err) */
    cc2 = mpn_add_1(tmp+bn-k, bp+bn-k, k, (mp_limb_t)1<<l);
    /* if carry, then all bits up to err were to 1, and we can round only
       if cc==0 and mpfr_round_raw2 returns 0 below */
    if (cc2 && cc) { TMP_FREE(marker); return 0; }
    cc2 = (tmp[bn-1]>>l1) & 1; /* gives 0 when carry */
    cc2 ^= mpfr_round_raw2(tmp, bn, neg, rnd2, prec);

    TMP_FREE(marker); 
    return (cc == cc2);

  case GMP_RNDU: /* b-2^(MPFR_EXP(b)-err) <= x <= b */
    tmp = TMP_ALLOC(tn*BYTES_PER_MP_LIMB); 
    /* first round b */
    cc = (bp[bn-1]>>l1) & 1;
    cc ^= mpfr_round_raw2(bp, bn, neg, rnd2, prec);

    /* now round b-2^(MPFR_EXP(b)-err) */
    cc2 = mpn_sub_1(tmp+bn-k, bp+bn-k, k, (mp_limb_t)1<<l);
    /* if borrow, then all bits up to err were to 0, and we can round only
       if cc==0 and mpfr_round_raw2 returns 1 below */
    if (cc2 && cc) { TMP_FREE(marker); return 0; }
    cc2 = (tmp[bn-1]>>l1) & 1; /* gives 1 when carry */
    cc2 ^= mpfr_round_raw2(tmp, bn, neg, rnd2, prec);

    TMP_FREE(marker); 
    return (cc == cc2);

  case GMP_RNDN: /* b-2^(MPFR_EXP(b)-err-1) <= x <= b+2^(MPFR_EXP(b)-err-1) */
    if (l==0) tn++; 
    tmp = TMP_ALLOC(tn*BYTES_PER_MP_LIMB); 

    /* this case is the same than GMP_RNDZ, except we first have to
       subtract 2^(MPFR_EXP(b)-err-1) from b */

    if (l) { 
      l--; /* tn=bn */
      mpn_sub_1(tmp+tn-k, bp+bn-k, k, (mp_limb_t)1<<l);
    } 
    else { 
      MPN_COPY(tmp+1, bp, bn); *tmp=0; /* extra limb to add or subtract 1 */
      k++; l=BITS_PER_MP_LIMB-1; 
      mpn_sub_1(tmp+tn-k, tmp+tn-k, k, (mp_limb_t)1<<l);
    }

    /* round b-2^(MPFR_EXP(b)-err-1) */
    /* we can disregard borrow, since we start from tmp in 2nd case too */
    cc = (tmp[tn-1]>>l1) & 1;
    cc ^= mpfr_round_raw2(tmp, tn, neg, rnd2, prec);

    if (l==BITS_PER_MP_LIMB-1) { l=0; k--; } else l++;

    /* round b+2^(MPFR_EXP(b)-err-1) = b-2^(MPFR_EXP(b)-err-1) + 2^(MPFR_EXP(b)-err) */    
    cc2 = mpn_add_1(tmp+tn-k, tmp+tn-k, k, (mp_limb_t)1<<l);
    /* if carry, then all bits up to err were to 1, and we can round only
       if cc==0 and mpfr_round_raw2 returns 0 below */
    if (cc2 && cc) { TMP_FREE(marker); return 0; }
    cc2 = (tmp[tn-1]>>l1) & 1; /* gives 0 when carry */
    cc2 ^= mpfr_round_raw2(tmp, tn, neg, rnd2, prec);

    TMP_FREE(marker); 
    return (cc == cc2);
  }
  return 0;
}
