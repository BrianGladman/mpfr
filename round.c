#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

#ifdef Exp
#include "longlong.h"
#endif

/* returns 0 if round(sign*xp[0..xn-1], prec, rnd) = 
   round(sign*xp[0..xn-1], prec, GMP_RNDZ), 1 otherwise,
   where sign=1 if neg=0, sign=-1 otherwise.

   Does *not* modify anything.
*/
int mpfr_round_raw2(mp_limb_t *xp, unsigned long xn, char neg, char rnd, 
		    unsigned long prec)
{
  unsigned long mask, nw; long wd; char rw; short l;

  nw = prec / BITS_PER_MP_LIMB; rw = prec & (BITS_PER_MP_LIMB - 1); 
  if (rw) nw++; 
  if (rnd==GMP_RNDZ || xn<nw || (rnd==GMP_RNDU && neg)
      || (rnd==GMP_RNDD && neg==0)) return 0;
  mask = ~((1UL<<(BITS_PER_MP_LIMB - rw)) - 1);
  switch (rnd)
    {
    case GMP_RNDU:
    case GMP_RNDD:
      if (xp[xn-nw] & ~mask) return 1;
      for (l=nw+1;l<=xn;l++)
	if (xp[xn-l]) break;
      return (l<=xn);
    case GMP_RNDN:
    /* First check if we are just halfway between two representable numbers */
      wd = xn - nw;
      if (!rw)
	{
	  wd--; 
	  if ((xp[wd] & (~mask)) == (1UL << (BITS_PER_MP_LIMB - rw - 1)))
	      do wd--; while (wd > 0 && !xp[wd]);

	  if (wd)
	      return ((xp[xn - nw - 1]>>(BITS_PER_MP_LIMB - 1)) & 1);
	  else
	      return xp[xn - nw] & 1;
	}
      else
      if (rw + 1 < BITS_PER_MP_LIMB)
	{
	  if ((xp[wd] & (~mask)) == (1UL << (BITS_PER_MP_LIMB - rw - 1)))
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

/* puts in y the value of xp (which has xsize limbs) rounded to precision
   xprec and direction RND_MODE */
int
mpfr_round_raw(mp_limb_t *y, mp_limb_t *xp, char RND_MODE, 
	       long xsize, unsigned long xprec)
{
  unsigned long nw, mask;
  char negative = 0, rw, carry = 0;

  if (xsize == 0) { return 0; }
  /* warning: -xsize is not correct since some bits are used for Nan, exact */
  if (xsize < 0) { negative = 1; xsize = xsize ^ (1<<31); }

#ifdef Exp
  count_leading_zeros(flag, xp[xsize-1]); 
  xprec += flag; 
#endif

  nw = xprec / BITS_PER_MP_LIMB; rw = xprec & (BITS_PER_MP_LIMB - 1); 
  if (rw) nw++; 
  /* number of words needed to represent x */
  /* precision is larger than the size of x. No rounding is necessary. 
     Just add zeroes at the end */
  if (xsize < nw) { 
    MPN_COPY(y + nw - xsize, xp, xsize); 
    MPN_ZERO(y, nw - xsize); /* PZ 27 May 99 */
    return 0; 
  }
  
  mask = ~((1UL<<(BITS_PER_MP_LIMB - rw)) - 1); 

  if (mpfr_round_raw2(xp, xsize, negative, RND_MODE, xprec))
    carry = mpn_add_1(y, xp + xsize - nw, nw,
                          1UL << (BITS_PER_MP_LIMB - rw));
  else MPN_COPY(y, xp + xsize - nw, nw);

  y[0] &= mask;
  return carry; 
}

void
mpfr_round(mpfr_t x, char RND_MODE, unsigned long prec)
{
  mp_limb_t *tmp; int carry; unsigned long nw; 

  nw = prec / BITS_PER_MP_LIMB; 
  if (prec & (BITS_PER_MP_LIMB - 1)) nw++;
  tmp = (mp_ptr) (*_mp_allocate_func) (nw * BYTES_PER_MP_LIMB);
  carry = mpfr_round_raw(tmp, x->_mp_d, RND_MODE, x->_mp_size, prec); 

  if (carry)
    {      
      mpn_rshift(tmp, tmp, nw, 1); 
      tmp [nw-1] |= (1UL << (BITS_PER_MP_LIMB - 1)); 
      EXP(x)++; 
    }

  (*_mp_free_func) (x->_mp_d, ABSSIZE(x) * BYTES_PER_MP_LIMB);
  if (SIGN(x)<0) { SIZE(x) = nw; CHANGE_SIGN(x); }
  else SIZE(x) = nw;
  PREC(x) = prec; 
  MANT(x) = tmp; 
}

/* hypotheses : BITS_PER_MP_LIMB est une puissance de 2 
                dans un test, la premiere partie du && est evaluee d'abord */


/* assuming b is an approximation of x in direction rnd1 
   with error at most 2^(EXP(b)-err), returns 1 if one is 
   able to round exactly x to precision prec with direction rnd2,
   and 0 otherwise.

   Side effects: none.
*/
int mpfr_can_round(b, err, rnd1, rnd2, prec) mpfr_t b; unsigned long err;
unsigned char rnd1, rnd2; unsigned long prec;
{
  return mpfr_can_round_raw(MANT(b), (PREC(b) - 1)/BITS_PER_MP_LIMB + 1, 
			    SIGN(b), err, rnd1, rnd2, prec); 
}

int
mpfr_can_round_raw(mp_limb_t *bp, unsigned long bn, char neg, 
		   unsigned long err, unsigned char rnd1, unsigned char rnd2, 
		   unsigned long prec)
{
  int k, k1, l, l1; mp_limb_t cc, cc2, *tmp;
  TMP_DECL(marker); 

  if (err<=prec) return 0;
  neg = (neg > 0 ? 0 : 1); 

  /* warning: if k = m*BITS_PER_MP_LIMB, consider limb m-1 and not m */
  k = (err-1)/BITS_PER_MP_LIMB;
  l = err % BITS_PER_MP_LIMB; if (l) l = BITS_PER_MP_LIMB-l;
  /* the error corresponds to bit l in limb k */
  k1 = (prec-1)/BITS_PER_MP_LIMB;
  l1 = prec%BITS_PER_MP_LIMB; if (l1) l1 = BITS_PER_MP_LIMB-l1;

  /* the last significant bit is bit l1 in limb k1 */

  if (rnd1==GMP_RNDU) { if (neg) rnd1=GMP_RNDZ; }
  if (rnd1==GMP_RNDD) { if (neg) rnd1=GMP_RNDU; else rnd1=GMP_RNDZ; }

  /* in the sequel, RNDU = towards infinity, RNDZ = towards zero */

  TMP_MARK(marker); 
  tmp = TMP_ALLOC(bn*BYTES_PER_MP_LIMB); 
  MPN_COPY(tmp, bp, bn); 

  switch (rnd1) {

  case GMP_RNDZ: /* b <= x <= b+2^(EXP(b)-err) */
    cc = (tmp[bn-k1-1]>>l1) & 1;
    cc ^= mpfr_round_raw2(tmp, bn, neg, rnd2, prec);

    /* now round b+2^(EXP(b)-err) */
    mpn_add_1(tmp+bn-k-1, tmp+bn-k-1, k+1, (mp_limb_t)1<<l);
    cc2 = (tmp[bn-k1-1]>>l1) & 1;
    cc2 ^= mpfr_round_raw2(tmp, bn, neg, rnd2, prec);

    /* if parity of cc and cc2 equals, then one is able to round */
    TMP_FREE(marker); 
    return (cc == cc2);

  case GMP_RNDU: /* b-2^(EXP(b)-err) <= x <= b */
    /* first round b */
    cc = (tmp[bn-k1-1]>>l1) & 1;
    cc ^= mpfr_round_raw2(tmp, bn, neg, rnd2, prec);

    /* now round b-2^(EXP(b)-err) */
    cc2 = mpn_sub_1(tmp+bn-k-1, tmp+bn-k-1, k+1, (mp_limb_t)1<<l);
    if (cc2) { TMP_FREE(marker); return 0; }
    cc2 = (tmp[bn-k1-1]>>l1) & 1;
    cc2 ^= mpfr_round_raw2(tmp, bn, neg, rnd2, prec);

    /* if parity of cc and cc2 equals, then one is able to round */
    TMP_FREE(marker); 
    return (cc == cc2);

  case GMP_RNDN: /* b-2^(EXP(b)-err-1) <= x <= b+2^(EXP(b)-err-1) */
    
    if (l) { l--; } else { k++; l=BITS_PER_MP_LIMB-1; }

    /* first round b+2^(EXP(b)-err-1)*/    
    cc = mpn_add_1(tmp+bn-k-1, tmp+bn-k-1, k+1, (mp_limb_t)1<<l);
    if (cc) { TMP_FREE(marker); return 0; }
    cc = (tmp[bn-k1-1]>>l1) & 1;
    cc ^= mpfr_round_raw2(tmp, bn, neg, rnd2, prec);

    MPN_COPY(tmp, bp, bn); 

    /* now round b-2^(EXP(b)-err-1) */
    cc2 = mpn_sub_1(tmp+bn-k-1, tmp+bn-k-1, k+1, (mp_limb_t)1<<l);
    if (cc2) { TMP_FREE(marker); return 0; }
    cc2 = (tmp[bn-k1-1]>>l1) & 1;
    cc2 ^= mpfr_round_raw2(tmp, bn, neg, rnd2, prec);

    TMP_FREE(marker); 
    return (cc == cc2);

  default:
    printf("rnd1=%d not yet implemented in mpfr_round2\n",rnd1);
    exit(1);
  }
}
