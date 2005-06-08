/* mpfr_div -- divide two floating-point numbers

Copyright 1999, 2001, 2002, 2003, 2004, 2005 Free Software Foundation.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA. */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

#ifdef DEBUG
#define mpn_print(ap,n) mpn_print3(ap,n,MPFR_LIMB_ZERO)
static void
mpn_print3 (mp_ptr ap, mp_size_t n, mp_limb_t cy)
{
  mp_size_t i;
  for (i = 0; i < n; i++)
    printf ("+%lu*2^%u", ap[i], BITS_PER_MP_LIMB * i);
  if (cy)
    printf ("+2^%u", BITS_PER_MP_LIMB * n);
  printf ("\n");
}
#endif

/* check if {ap, an} is zero */
static int
mpn_cmpzero (mp_ptr ap, mp_size_t an)
{
  while (an > 0)
    if (MPFR_LIKELY(ap[--an] != MPFR_LIMB_ZERO))
      return 1;
  return 0;
}

/* compare {ap, an} and {bp, bn} >> extra,
   aligned by the more significant limbs.
   Takes into account bp[0] for extra=1.
*/
static int
mpn_cmp_aux (mp_ptr ap, mp_size_t an, mp_ptr bp, mp_size_t bn, int extra)
{
  int cmp = 0;
  mp_size_t k;
  mp_limb_t bb;

  if (an >= bn)
    {
      k = an - bn;
      while (cmp == 0 && bn > 0)
	{
	  bn --;
	  bb = (extra) ? ((bp[bn+1] << (BITS_PER_MP_LIMB - 1)) | (bp[bn] >> 1))
	    : bp[bn];
	  cmp = (ap[k + bn] > bb) ? 1 : ((ap[k + bn] < bb) ? -1 : 0);
	}
      bb = (extra) ? bp[0] << (BITS_PER_MP_LIMB - 1) : MPFR_LIMB_ZERO;
      while (cmp == 0 && k > 0)
	{
	  k--;
	  cmp = (ap[k] > bb) ? 1 : ((ap[k] < bb) ? -1 : 0);
	  bb = MPFR_LIMB_ZERO; /* ensure we consider only once bp[0] & 1 */
	}
      if (cmp == 0 && bb != MPFR_LIMB_ZERO)
	cmp = -1;
    }
  else /* an < bn */
    {
      k = bn - an;
      while (cmp == 0 && an > 0)
        {
          an --;
	  bb = (extra) ? ((bp[k+an+1] << (BITS_PER_MP_LIMB - 1)) | (bp[k+an] >> 1))
            : bp[k+an];
          if (ap[an] > bb)
            cmp = 1;
          else if (ap[an] < bb)
            cmp = -1;
        }
      while (cmp == 0 && k > 0)
	{
	  k--;
	  bb = (extra) ? ((bp[k+1] << (BITS_PER_MP_LIMB - 1)) | (bp[k] >> 1))
            : bp[k];
	  cmp = (bb != MPFR_LIMB_ZERO) ? -1 : 0;
	}
      if (cmp == 0 && extra && (bp[0] & MPFR_LIMB_ONE))
	cmp = -1;
    }
  return cmp;
}

/* {ap, n} <- {ap, n} - {bp, n} >> extra - cy, with cy=0 or 1 */
static mp_limb_t
mpn_sub_aux (mp_ptr ap, mp_ptr bp, mp_size_t n, mp_limb_t cy, int extra)
{
  mp_limb_t bb, rp;
  while (n--)
    {
      bb = (extra) ? ((bp[1] << (BITS_PER_MP_LIMB-1)) | (bp[0] >> 1)) : bp[0];
      rp = ap[0] - bb - cy;
      cy = ((ap[0] < bb) || (cy && ~rp == MPFR_LIMB_ZERO)) ? MPFR_LIMB_ONE : MPFR_LIMB_ZERO;
      ap[0] = rp;
      ap ++;
      bp ++;
    }
  return cy;
}

int
mpfr_div (mpfr_ptr q, mpfr_srcptr u, mpfr_srcptr v, mp_rnd_t rnd_mode)
{
  mp_size_t q0size = MPFR_LIMB_SIZE(q); /* number of limbs of destination */
  mp_size_t usize = MPFR_LIMB_SIZE(u);
  mp_size_t vsize = MPFR_LIMB_SIZE(v);
  mp_size_t qsize; /* number of limbs of the computed quotient */
  mp_size_t qqsize;
  mp_size_t k, l;
  mp_ptr q0p = MPFR_MANT(q), qp;
  mp_ptr up = MPFR_MANT(u);
  mp_ptr vp = MPFR_MANT(v);
  mp_ptr ap;
  mp_ptr bp;
  mp_limb_t qh;
  mp_limb_t sticky_u = MPFR_LIMB_ZERO;
  mp_limb_t low_u;
  mp_limb_t sticky_v = MPFR_LIMB_ZERO;
  mp_limb_t sticky;
  mp_limb_t sticky3;
  mp_limb_t round_bit = MPFR_LIMB_ZERO;
  mp_exp_t qexp;
  int sign_quotient;
  int extra_bit;
  int sh, sh2;
  int inex;
  MPFR_TMP_DECL(marker);

  /**************************************************************************
   *                                                                        *
   *              This part of the code deals with special cases            *
   *                                                                        *
   **************************************************************************/

  if (MPFR_UNLIKELY(MPFR_ARE_SINGULAR(u,v)))
    {
      if (MPFR_IS_NAN(u) || MPFR_IS_NAN(v))
	{
	  MPFR_SET_NAN(q);
	  MPFR_RET_NAN;
	}
      sign_quotient = MPFR_MULT_SIGN( MPFR_SIGN(u) , MPFR_SIGN(v) );
      MPFR_SET_SIGN(q, sign_quotient);
      if (MPFR_IS_INF(u))
	{
	  if (MPFR_IS_INF(v))
	    {
	      MPFR_SET_NAN(q);
	      MPFR_RET_NAN;
	    }
	  else
	    {
	      MPFR_SET_INF(q);
	      MPFR_RET(0);
	    }
	}
      else if (MPFR_IS_INF(v))
	{
	  MPFR_SET_ZERO (q);
	  MPFR_RET (0);
	}
      else if (MPFR_IS_ZERO (v))
	{
	  if (MPFR_IS_ZERO (u))
	    {
	      MPFR_SET_NAN(q);
	      MPFR_RET_NAN;
	    }
	  else
	    {
	      MPFR_SET_INF(q);
	      MPFR_RET(0);
	    }
	}
      else
	{
	  MPFR_ASSERTD (MPFR_IS_ZERO (u));
	  MPFR_SET_ZERO (q);
	  MPFR_RET (0);
	}
    }
  MPFR_CLEAR_FLAGS (q);

  /**************************************************************************
   *                                                                        *
   *              End of the part concerning special values.                *
   *                                                                        *
   **************************************************************************/

  MPFR_TMP_MARK(marker);

  /* set sign */
  sign_quotient = MPFR_MULT_SIGN( MPFR_SIGN(u) , MPFR_SIGN(v) );
  MPFR_SET_SIGN(q, sign_quotient);

  /* determine if an extra bit comes from the division, i.e. if the
     significand of u (as a fraction in [1/2, 1[) is larger than that
     of v */
  if (MPFR_LIKELY(up[usize - 1] != vp[vsize - 1]))
    extra_bit = (up[usize - 1] > vp[vsize - 1]) ? 1 : 0;
  else /* most significant limbs are equal, must look at further limbs */
    {
      k = usize - 1;
      l = vsize - 1;
      while (k != 0 && l != 0 && up[--k] == vp[--l]);
      /* now k=0 or l=0 or up[k] != vp[l] */
      if (up[k] > vp[l])
	extra_bit = 1;
      else if (up[k] < vp[l])
	extra_bit = 0;
      /* now up[k] = vp[l], thus either k=0 or l=0 */
      else if (l == 0) /* no more divisor limb */
	extra_bit = 1;
      else /* k=0: no more dividend limb */
        extra_bit = mpn_cmpzero (vp, l) == 0;
    }
#ifdef DEBUG
  printf ("extra_bit=%u\n", extra_bit);
#endif

  /* set exponent */
  qexp = MPFR_GET_EXP (u) - MPFR_GET_EXP (v) + extra_bit;

  MPFR_UNSIGNED_MINUS_MODULO(sh, MPFR_PREC(q));

  if (MPFR_UNLIKELY(rnd_mode == GMP_RNDN && sh == 0))
    { /* we compute the quotient with one more limb, in order to get
	 the round bit in the quotient, and the remainder only contains
	 sticky bits */
      qsize = q0size + 1;
      /* need to allocate memory for the quotient */
      qp = (mp_ptr) MPFR_TMP_ALLOC (qsize*sizeof(mp_limb_t));
    }
  else
    {
      qsize = q0size;
      qp = q0p; /* directly put the quotient in the destination */
    }
  qqsize = qsize + qsize;

  /* prepare the dividend */
  ap = (mp_ptr) MPFR_TMP_ALLOC (qqsize*sizeof(mp_limb_t));
  if (MPFR_LIKELY(qqsize > usize)) /* use the full dividend */
    {
      k = qqsize - usize; /* k > 0 */
      MPN_ZERO(ap, k);
      if (extra_bit)
	ap[k - 1] = mpn_rshift (ap + k, up, usize, 1);
      else
	MPN_COPY(ap + k, up, usize);
    }
  else /* truncate the dividend */
    {
      k = usize - qqsize;
      if (extra_bit)
	sticky_u = mpn_rshift (ap, up + k, qqsize, 1);
      else
	MPN_COPY(ap, up + k, qqsize);
      sticky_u = sticky_u || mpn_cmpzero (up, k);
    }
  low_u = sticky_u;
  
  /* now sticky_u is non-zero iff the truncated part of u is non-zero */

  /* prepare the divisor */
  if (MPFR_LIKELY(vsize >= qsize))
    {
      k = vsize - qsize;
      if (qp != vp)
	bp = vp + k; /* avoid copying the divisor */
      else /* need to copy, since mpn_divrem doesn't allow overlap
	      between quotient and divisor, necessarily k = 0
	      since quotient and divisor are the same mpfr variable */
	{
	  bp = (mp_ptr) MPFR_TMP_ALLOC (qsize * sizeof(mp_limb_t));
	  MPN_COPY(bp, vp, vsize);
	}
      sticky_v = sticky_v || mpn_cmpzero (vp, k);
    }
  else /* vsize < qsize */
    {
      k = qsize - vsize;
      bp = (mp_ptr) MPFR_TMP_ALLOC (qsize * sizeof(mp_limb_t));
      MPN_COPY(bp + k, vp, vsize);
      MPN_ZERO(bp, k);
    }

  /* we now can perform the division */
  qh = mpn_divrem (qp, 0, ap, qqsize, bp, qsize);
  /* warning: qh may be 1 if u1 == v1, but u < v */
#ifdef DEBUG
  printf ("q="); mpn_print (qp, qsize);
  printf ("r="); mpn_print (ap, qsize);
#endif

  k = qsize;
  sticky_u = sticky_u || mpn_cmpzero (ap, k);

  sticky = sticky_u | sticky_v;

  /* now sticky is non-zero iff one of the following holds:
     (a) the truncated part of u is non-zero
     (b) the truncated part of v is non-zero
     (c) the remainder from division is non-zero */

  if (MPFR_LIKELY(qsize == q0size))
    {
      sticky3 = qp[0] & MPFR_LIMB_MASK(sh); /* does nothing when sh=0 */
      sh2 = sh;
    }
  else /* qsize = q0size + 1: only happens when rnd_mode=GMP_RNDN and sh=0 */
    {
      MPN_COPY (q0p, qp + 1, q0size);
      sticky3 = qp[0];
      sh2 = BITS_PER_MP_LIMB;
    }
  qp[0] ^= sticky3;
  /* sticky3 contains the truncated bits from the quotient,
     including the round bit, and 1 <= sh2 <= BITS_PER_MP_LIMB
     is the number of bits in sticky3 */
  inex = (sticky != MPFR_LIMB_ZERO) || (sticky3 != MPFR_LIMB_ZERO);
#ifdef DEBUG
  printf ("sticky=%lu sticky3=%lu inex=%d\n", sticky, sticky3, inex);
#endif

  if (sign_quotient < 0)
    rnd_mode = MPFR_INVERT_RND(rnd_mode);

  /* to round, we distinguish two cases:
     (a) vsize <= qsize: we used the full divisor
     (b) vsize > qsize: the divisor was truncated
  */

#ifdef DEBUG
  printf ("vsize=%u qsize=%u\n", vsize, qsize);
#endif
  if (MPFR_LIKELY(vsize <= qsize)) /* use the full divisor */
    {
      if (MPFR_LIKELY(rnd_mode == GMP_RNDN))
	{
	  round_bit = sticky3 & (MPFR_LIMB_ONE << (sh2 - 1));
	  sticky = (sticky3 ^ round_bit) | sticky_u;
	}
      else if (rnd_mode == GMP_RNDZ || rnd_mode == GMP_RNDD || inex == MPFR_LIMB_ZERO)
	sticky = (inex == 0) ? MPFR_LIMB_ZERO : MPFR_LIMB_ONE;
      else /* rnd_mode = GMP_RNDU */
	sticky = MPFR_LIMB_ONE;
      goto case_1;
    }
  else /* vsize > qsize: need to truncate the divisor */
    {
      if (inex == MPFR_LIMB_ZERO)
	goto truncate;
      else
	{
	  /* we can round except when sticky3 is 000...000 or 000...001
	     for directed rounding, and 100...000 or 100...001 for rounding
	     to nearest. (For rounding to nearest, we cannot determine the
	     inexact flag for 000...000 or 000...001.)
	  */
	  mp_limb_t sticky3orig = sticky3;
	  if (rnd_mode == GMP_RNDN)
	    {
	      round_bit = sticky3 & (MPFR_LIMB_ONE << (sh2 - 1));
	      sticky3   = sticky3 ^ round_bit;
#ifdef DEBUG
              printf ("rb=%lu sb=%lu\n", round_bit, sticky3);
#endif
	    }
	  if (sticky3 != MPFR_LIMB_ZERO && sticky3 != MPFR_LIMB_ONE)
	    {
	      sticky = sticky3;
	      goto case_1;
	    }
	  else /* hard case: we have to compare q1 * v0 and r + low(u),
		 where q1 * v0 has qsize + (vsize-qsize) = vsize limbs, and
		 r + low(u) has qsize + (usize-2*qsize) = usize-qsize limbs */
	    {
	      mp_size_t l;
	      mp_ptr sp;
	      int cmp_s_r;

	      sp = (mp_ptr) MPFR_TMP_ALLOC (vsize*sizeof(mp_limb_t));
	      k = vsize - qsize;
	      /* sp <- {qp, qsize} * {vp, vsize-qsize} */
	      qp[0] ^= sticky3orig; /* restore original quotient */
	      if (qsize >= k)
		mpn_mul (sp, qp, qsize, vp, k);
	      else
		mpn_mul (sp, vp, k, qp, qsize);
              if (qh)
                mpn_add_n (sp + qsize, sp + qsize, vp, k);
	      qp[0] ^= sticky3orig; /* restore truncated quotient */

	      /* compare {sp, vsize = k + qsize} to {ap, qsize} + low(u) */
	      cmp_s_r = mpn_cmp (sp + k, ap, qsize);
	      if (cmp_s_r == 0) /* compare {sp, k} and low(u) */
                {
                  cmp_s_r = (usize >= qqsize)
                             ? mpn_cmp_aux (sp, k, up, usize-qqsize, extra_bit)
                             : mpn_cmpzero (sp, k);
                }
#ifdef DEBUG
              printf ("cmp(q*v0,r+u0)=%d\n", cmp_s_r);
#endif
	      /* now cmp_s_r > 0 if {sp, vsize} > {ap, qsize} + low(u)
		     cmp_s_r = 0 if {sp, vsize} = {ap, qsize} + low(u) 
		     cmp_s_r < 0 if {sp, vsize} < {ap, qsize} + low(u) */
	      if (cmp_s_r <= 0) /* quotient is in [q1, q1+1) */
		{
		  sticky = (cmp_s_r == 0) ? sticky3 : MPFR_LIMB_ONE;
		  goto case_1;
		}
	      else /* cmp_s_r > 0, quotient is < q1 */
		{
		  mp_limb_t cy = MPFR_LIMB_ZERO;
		  /* subtract low(u)>>extra_bit if non-zero */
		  if (low_u != MPFR_LIMB_ZERO)
		    {
		      mp_size_t m;
		      l = usize - qqsize; /* number of low limbs in u */
		      m = (l > k) ? l - k : 0;
		      cy = (extra_bit) ? (up[m] & MPFR_LIMB_ONE) : MPFR_LIMB_ZERO;
		      if (l >= k) /* u0 has more limbs */
			{
                          cy = cy || mpn_cmpzero (up, m);
			  low_u = cy;
			  cy = mpn_sub_aux (sp, up + l - k, k,
					    (cy) ? MPFR_LIMB_ONE : MPFR_LIMB_ZERO, extra_bit);
			}
		      else /* l < k: s has more limbs than u0 */
			{
			  low_u = MPFR_LIMB_ZERO;
			  if (cy != MPFR_LIMB_ZERO)
			    cy = mpn_sub_1 (sp + k - l - 1, sp + k - l - 1, 1, MPFR_LIMB_HIGHBIT);
			  cy = mpn_sub_aux (sp + k - l, up, l, cy, extra_bit);
			}
		    }
		  cy = mpn_sub_1 (sp + k, sp + k, qsize, cy);
		  /* subtract r */
		  cy = mpn_sub_nc (sp + k, sp + k, ap, qsize, cy);
		  /* now compare {sp, ssize} to v */
		  cmp_s_r = mpn_cmp (sp, vp, vsize);
		  if (cmp_s_r == 0 && low_u != MPFR_LIMB_ZERO)
		    cmp_s_r = 1; /* since in fact we subtracted less than 1 */
#ifdef DEBUG
                  printf ("cmp(q*v0-(r+u0),v)=%d\n", cmp_s_r);
#endif
		  if (cmp_s_r <= 0) /* q1-1 <= u/v < q1 */
		    {
		      if (sticky3 == MPFR_LIMB_ONE)
			{ /* q1-1 is either representable (directed rounding),
                             or the middle of two numbers (nearest) */
                          sticky = (cmp_s_r) ? MPFR_LIMB_ONE : MPFR_LIMB_ZERO;
                          goto case_1;
			}
                      /* now necessarily sticky3=0 */
		      else if (round_bit == MPFR_LIMB_ZERO)
			{ /* round_bit=0, sticky3=0: q1-1 is exact only
                             when sh=0 */
			  inex = (cmp_s_r || sh) ? -1 : 0;
			  if ((rnd_mode == GMP_RNDU && inex != 0)
			      || rnd_mode == GMP_RNDN)
                            {
                              inex = 1;
                              goto truncate_check_qh;
                            }
			  else /* round down */
                            goto sub_one_ulp;
			}
		      else /* sticky3=0, round_bit=1 ==> rounding to nearest */
			{
			  inex = cmp_s_r;
			  goto truncate;
			}
		    }
		  else /* q1-2 < u/v < q1-1 */
		    {
		      /* if rnd=GMP_RNDU, the result is up(q1-1),
			 which is q1 unless sh = 0, where it is q1-1 */
		      if (rnd_mode == GMP_RNDU)
			{
			  inex = 1;
			  if (sh > 0)
			    goto truncate_check_qh;
			  else /* sh = 0 */
			    goto sub_one_ulp;
			}
		      /* if rnd=GMP_RNDN, the result is q1 when
			 q1-2 >= q1-2^(sh-1), i.e. sh >= 2,
			 otherwise (sh=1) it is q1-2 */
		      else if (rnd_mode == GMP_RNDN) /* sh > 0 */
			{
			  /* Case sh=1: sb=0 always, and q1-rb is exactly
			     representable, like q1-rb-2.
			     rb action
			     0  subtract two ulps, inex=-1
                             1  truncate, inex=1

			     Case sh>1: one ulp is 2^(sh-1) >= 2
			     rb sb action
                             0  0  truncate, inex=1
                             0  1  truncate, inex=1
                             1  x  truncate, inex=-1
			   */
			  if (sh == 1)
			    {
			      if (round_bit == MPFR_LIMB_ZERO)
				{
				  inex = -1;
				  sh = 0;
				  goto sub_two_ulp;
				}
			      else
				{
				  inex = 1;
				  goto truncate_check_qh;
				}
			    }
			  else /* sh > 1 */
			    {
                              inex = (round_bit == MPFR_LIMB_ZERO) ? 1 : -1;
                              goto truncate_check_qh;
			    }
			}
		      else /* round down */
			{
			  /* the result is down(q1-2), i.e. subtract one
			     ulp if sh > 0, and two ulps if sh=0 */
			  inex = -1;
			  if (sh > 0)
			    goto sub_one_ulp;
			  else
			    goto sub_two_ulp;
			}
		    }
		}
	    }
	}
    }

 case_1: /* quotient is in [q1, q1+1),
	    round_bit is the round_bit (0 for directed rounding),
	    sticky the sticky bit */
  if (rnd_mode == GMP_RNDZ || rnd_mode == GMP_RNDD || 
      (round_bit == MPFR_LIMB_ZERO && sticky == MPFR_LIMB_ZERO))
    {
      inex = (round_bit == MPFR_LIMB_ZERO && sticky == MPFR_LIMB_ZERO) ? 0 : -1;
      goto truncate;
    }
  else if (rnd_mode == GMP_RNDN) /* sticky <> 0 or round <> 0 */
    {
      if (round_bit == MPFR_LIMB_ZERO) /* necessarily sticky <> 0 */
	{
	  inex = -1;
	  goto truncate;
	}
      /* round_bit = 1 */
      else if (sticky != MPFR_LIMB_ZERO)
	goto add_one_ulp; /* inex=1 */
      else /* round_bit=1, sticky=0 */
	goto even_rule;
    }
  else /* rnd_mode = GMP_RNDU, sticky <> 0 */
    goto add_one_ulp; /* with inex=1 */

 sub_two_ulp:
  /* we cannot subtract MPFR_LIMB_MPFR_LIMB_ONE << (sh+1) since this is 
     undefined for sh = BITS_PER_MP_LIMB */
  qh -= mpn_sub_1 (q0p, q0p, q0size, MPFR_LIMB_ONE << sh);
  /* go through */

 sub_one_ulp:
  qh -= mpn_sub_1 (q0p, q0p, q0size, MPFR_LIMB_ONE << sh);
  /* go through truncate_check_qh */

 truncate_check_qh:
  if (qh)
    {
      qexp ++;
      q0p[q0size - 1] = MPFR_LIMB_HIGHBIT;
    }
  goto truncate;

 even_rule: /* has to set inex */
  inex = (q0p[0] & (MPFR_LIMB_ONE << sh)) ? 1 : -1;
  if (inex < 0)
    goto truncate;
  /* else go through add_one_ulp */

 add_one_ulp:
  inex = 1; /* always here */
  if (mpn_add_1 (q0p, q0p, q0size, MPFR_LIMB_ONE << sh))
    {
      qexp ++;
      q0p[q0size - 1] = MPFR_LIMB_HIGHBIT;
    }
  
 truncate: /* inex already set */

  MPFR_TMP_FREE(marker);

  /* check for underflow/overflow */
  if (MPFR_UNLIKELY(qexp > __gmpfr_emax))
    return mpfr_overflow (q, rnd_mode, sign_quotient);
  else if (MPFR_UNLIKELY(qexp < __gmpfr_emin))
    {
      if (rnd_mode == GMP_RNDN && ((qexp < __gmpfr_emin - 1) ||
                                   (inex == 0 && mpfr_powerof2_raw (q))))
        rnd_mode = GMP_RNDZ;
      return mpfr_underflow (q, rnd_mode, sign_quotient);
    }
  MPFR_SET_EXP(q, qexp);

  MPFR_RET (inex*sign_quotient);
}
