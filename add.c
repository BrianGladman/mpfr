/* mpfr_add -- add two floating-point numbers

Copyright (C) 1999-2001 Free Software Foundation.
Contributed by the Spaces project, INRIA Lorraine.

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

#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"

#define ONE ((mp_limb_t) 1)

extern void mpfr_sub1 _PROTO((mpfr_ptr, mpfr_srcptr, mpfr_srcptr,
                              mp_rnd_t, mp_exp_unsigned_t));
void mpfr_add1 _PROTO((mpfr_ptr, mpfr_srcptr, mpfr_srcptr,
                       mp_rnd_t, mp_exp_unsigned_t));

/* signs of b and c are supposed equal,
   diff_exp is the difference between the exponents of b and c,
   which is supposed >= 0 */

void
#if __STDC__
mpfr_add1 (mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c,
	   mp_rnd_t rnd_mode, mp_exp_unsigned_t diff_exp)
#else
mpfr_add1 (a, b, c, rnd_mode, diff_exp)
     mpfr_ptr a;
     mpfr_srcptr b;
     mpfr_srcptr c;
     mp_rnd_t rnd_mode;
     mp_exp_unsigned_t diff_exp;
#endif
{
  mp_limb_t *ap, *bp, *cp;
  mp_prec_t aq, bq, cq;
  mp_size_t an, bn, cn, k;
  int sh;
  int nulp = 1;
  TMP_DECL(marker);

  MPFR_ASSERTN(MPFR_NOTZERO(c));

  TMP_MARK(marker);
  ap = MPFR_MANT(a);
  bp = MPFR_MANT(b);
  cp = MPFR_MANT(c);

  if (ap == bp)
  {
    bp = (mp_ptr) TMP_ALLOC(MPFR_ABSSIZE(b) * BYTES_PER_MP_LIMB);
    MPN_COPY (bp, ap, MPFR_ABSSIZE(b));
    if (ap == cp) { cp = bp; }
  }
  else if (ap == cp)
  {
    cp = (mp_ptr) TMP_ALLOC (MPFR_ABSSIZE(c) * BYTES_PER_MP_LIMB);
    MPN_COPY(cp, ap, MPFR_ABSSIZE(c));
  }

  aq = MPFR_PREC(a);
  bq = MPFR_PREC(b);
  cq = MPFR_PREC(c);

  an = (aq-1)/BITS_PER_MP_LIMB + 1; /* number of significant limbs of a */
  sh = an*BITS_PER_MP_LIMB - aq;    /* non-significant bits in low limb */
  bn = (bq-1)/BITS_PER_MP_LIMB + 1; /* number of significant limbs of b */
  cn = (cq-1)/BITS_PER_MP_LIMB + 1; /* number of significant limbs of c */

  MPFR_EXP(a) = MPFR_EXP(b);
  MPFR_SET_SAME_SIGN(a, b);

  if (aq <= diff_exp)
  {

    /* diff_exp>=MPFR_PREC(a): c does not overlap with a */
    /* either MPFR_PREC(b)<=MPFR_PREC(a), and we can copy the mantissa of b
       directly into that of a, or MPFR_PREC(b)>MPFR_PREC(a) and we have to
       round b+c */

    if (bq <= aq)
    {

      MPN_COPY(ap + (an - bn), bp, bn);

      /* fill low significant limbs with zero */
      MPN_ZERO(ap, an - bn);

      /* now take c into account (c != 0) */
      if (rnd_mode == GMP_RNDN)
      {
	/* to nearest */
	/* if diff_exp > MPFR_PREC(a), no change */

	if (aq == diff_exp)
	{
	  /* as c is normalized, we have to add one to the lsb of a
             if c>1/2, or c=1/2 and lsb(a)=1 (round to even) */

          mp_limb_t cc;
          mp_limb_t *cp2 = cp + (cn-1);  /* highest limb */

          cc = *cp2 - MP_LIMB_T_HIGHBIT;
          while (cc == 0 && cp2 > cp) cc = *--cp2;

          if (cc || (ap[0] & (ONE<<sh))) goto add_one_ulp;
          /* mant(c) != 1/2 or mant(c) = 1/2: add 1 iff lsb(a)=1 */
        } /* aq == diff_exp */
      } /* rnd_mode == GMP_RNDN */
      else if ((rnd_mode == GMP_RNDU && MPFR_ISNONNEG(b)) ||
               (rnd_mode == GMP_RNDD && MPFR_ISNEG(b)))  /* round up */
        goto add_one_ulp;
      /* in the other cases (round to zero, or up/down with sign -/+),
         nothing to do */
    } /* bq <= aq */
    else /* bq > aq */
    {
      mp_limb_t inv, bb = 0, cc = 0;
      int difs, r = 0;
      mp_exp_t difw;  /* mp_exp_t may be larger than mp_size_t */

      /* MPFR_PREC(b)>MPFR_PREC(a) : we have to round b+c */
      k = bn - an;

      /* first copy the 'an' most significant limbs of b to a */
      MPN_COPY(ap, bp + k, an);

      difs = diff_exp % BITS_PER_MP_LIMB;
      difw = an - (mp_exp_t) (diff_exp / BITS_PER_MP_LIMB);
      MPFR_ASSERTN(difw <= 1);

      nulp = -1;

      inv = rnd_mode == GMP_RNDN ? (sh ? ONE<<(sh-1) : MP_LIMB_T_HIGHBIT) : 0;
      if (sh)
      {
        bb = *ap & ((ONE<<sh)-1);
        *ap &= ~bb; /* truncate last bits */
        if (inv)
        {
          r = bb >> (sh-1);
          bb ^= inv;
          inv = 0;
        }
        if (difw > 0)
        {
          /* c overlaps with the lowest limb of a */
          MPFR_ASSERTN(difs > 0);
          cc = cp[--cn];
          bb += cc >> difs;
        }
        if (bb >= ONE<<sh)
        {
          nulp = 1;
          bb -= ONE<<sh;
        }
        else if (bb < (ONE<<sh)-1)
        {
          nulp = 0;
        }
      }

      if (nulp < 0 || (nulp > 0 && bb == 0))
        while (1) /* look for carry and/or sticky bit */
        {
          if (k == 0)
          {
            if (nulp < 0)
              nulp = 0;
            else /* nulp == 1 */
            {
              if (difw > 0)
                bb = cc << (BITS_PER_MP_LIMB - difs);
              while (bb == 0 && cn)
                bb = cp[--cn];
            }
            break;
          }
          bb = bp[--k];
          if (inv)
          {
            r = bb >> (BITS_PER_MP_LIMB - 1);
            bb ^= inv;
            inv = 0;
          }
          if (difw >= 0)
          {
            mp_limb_t b2;
            if (difw > 0)
            {
              b2 = bb;
              bb = b2 + (cc << (BITS_PER_MP_LIMB - difs));
              if (bb < b2)
                nulp = 1;
            }
            if (cn == 0)
            {
              if (nulp < 0)
                nulp = 0;
              else /* nulp == 1 */
                while (bb == 0 && k)
                  bb = bp[--k];
              break;
            }
            cc = cp[--cn];
            b2 = bb;
            if (difs)
            {
              difw = 1;
              bb = b2 + (cc >> difs);
            }
            else
            {
              bb = b2 + cc;
            }
            if (bb < b2)
              nulp = 1;
          }
          else
            difw++;
          if (nulp > 0 && bb != 0)
            break;
          if (nulp < 0 && bb != (mp_limb_t) (-1))
          {
            nulp = 0;
            break;
          }
        } /* while */

      /* if nulp == 1, sticky bit = bb != 0 */

      if (rnd_mode == GMP_RNDN)
      { /* r: old rounding bit of b
           ulp(a)/2 were added if r = 0
           ulp(a)/2 were subtracted if r = 1
             r    carry  stbit  -->  nulp to add
             0      0      1              0
             0      1      0         0 or 1 (to even)
             0      1      1              1
             1      0      1              1
             1      1      0         1 or 2 (to even)
             1      1      1              2 [*]
           [*] or only 1 if adding the first ulp changes the exponent */

        MPFR_ASSERTN(inv == 0);  /* r has been initialized */
        if (nulp == 0)
          bb = 1;
        nulp += r;
        if (!bb && (ap[0] & (ONE<<sh)))
          nulp++;  /* sticky bit is 0 --> round to even */
      }
      else
      {
        if ((nulp == 0 || bb != 0) &&
               ((rnd_mode == GMP_RNDU && MPFR_ISNONNEG(b))
             || (rnd_mode == GMP_RNDD && MPFR_ISNEG(b))))
          nulp++;
      }
      goto add_one_ulp;
    } /* bq > aq */
  } /* aq <= diff_exp */
  else /* aq > diff_exp */
  { /* diff_exp < MPFR_PREC(a) : c overlaps with a */

    mp_exp_t dif;
    mp_limb_t cc, c2 = 0, c3 = 0;
    unsigned int overlap;

    /* first copy upper part of c into a (after shift) */
    dif = aq - diff_exp;
    /* dif is the number of bits of c which overlap with a */
    MPFR_ASSERTN(dif > 0);
    k = (dif-1)/BITS_PER_MP_LIMB + 1; /* only the highest k limbs from c
					 have to be considered */
    MPN_ZERO(ap+k, an-k);

    overlap = dif <= cq;
    if (overlap)
    {
      /* c has to be truncated */
      dif = dif % BITS_PER_MP_LIMB;
      dif = dif ? BITS_PER_MP_LIMB - dif - sh : -sh;

      /* we have to shift by dif bits to the right */

      if (dif > 0)
        mpn_rshift(ap, cp+(cn-k), k, dif);
      else if (dif < 0)
      {
	ap[k] = mpn_lshift(ap, cp+(cn-k), k, -dif);

	/* put the non-significant bits in low limb for further rounding */

	if (cn >= k+1)
	  ap[0] += cp[cn-k-1]>>(BITS_PER_MP_LIMB+dif);
      }
      else
        MPN_COPY(ap, cp+(cn-k), k);
    }
    else /* dif > cq */
    {
      /* c is not truncated, but we have to fill low limbs with 0 */

      int shift;

      k = diff_exp / BITS_PER_MP_LIMB;
      shift = diff_exp % BITS_PER_MP_LIMB;

      /* warning: a shift of zero bit is not allowed */
      MPN_ZERO(ap, an-k-cn);
      if (shift)
      {
	cc=mpn_rshift(ap + (an-k-cn), cp, cn, shift);
	if (an-k-cn > 0)
          ap[an-k-cn-1] = cc;
      }
      else
        MPN_COPY(ap + (an-k-cn), cp, cn);
    }

    /* here overlap=1 iff ulp(c)<ulp(a) */
    /* then put high limbs to zero */
    /* now add 'an' upper limbs of b in place */

    if (bq <= aq)
    {
      overlap += 2;
      cc = mpn_add_n(ap+(an-bn), ap+(an-bn), bp, bn);
    }
    else
    {
      /* MPFR_PREC(b) > MPFR_PREC(a): we have to truncate b */
      cc = mpn_add_n(ap, ap, bp+(bn-an), an);
    }

    if (cc) {

      /* shift one bit to the right */

      c3 = (ap[0]&1) && (MPFR_PREC(a)%BITS_PER_MP_LIMB==0);
      mpn_rshift(ap, ap, an, 1);
      ap[an-1] += MP_LIMB_T_HIGHBIT;
      MPFR_EXP(a)++;
    }

    /* remains to do the rounding */

    if (rnd_mode==GMP_RNDN || (MPFR_ISNONNEG(b) && rnd_mode==GMP_RNDU)
	|| (MPFR_ISNEG(b) && rnd_mode==GMP_RNDD)) {

      int kc;

      /* four cases: overlap =
         (0) MPFR_PREC(b) > MPFR_PREC(a) and diff_exp+MPFR_PREC(c) <= MPFR_PREC(a)
         (1) MPFR_PREC(b) > MPFR_PREC(a) and diff_exp+MPFR_PREC(c) > MPFR_PREC(a)
         (2) MPFR_PREC(b) <= MPFR_PREC(a) and diff_exp+MPFR_PREC(c) <= MPFR_PREC(a)
         (3)  MPFR_PREC(b) <= MPFR_PREC(a) and diff_exp+MPFR_PREC(c) > MPFR_PREC(a) */

      switch (overlap)
	{ mp_limb_t cout;
        case 1: /* both b and c to round */
	  kc = cn-k; /* remains kc limbs from c */
	  k = bn-an; /* remains k limbs from b */

	  /* truncate last bits and store the difference with 1/2*ulp in cc */

	  cc = *ap & ((ONE<<sh)-1);
	  *ap &= ~cc; /* truncate last bits */
	  if (rnd_mode==GMP_RNDN)
	    cout = -mpn_sub_1(&cc, &cc, 1, ONE<<(sh-1));
	  else cout=0;
	  if ((~cout==0) && (~cc)) break;
	  cout = cc;
	  while ((cout==0 || cout==-1) && k!=0 && kc!=0) {
	    kc--;
	    cout += mpn_add_1(&cc, bp+(--k), 1,(cp[kc+1]<<(BITS_PER_MP_LIMB-dif))
			    +(cp[kc]>>dif));
	    if (cout==0 || (~cout==0)) cout=cc;
	  }
	  if (kc==0 && dif) {
	    /* it still remains cp[0]<<(BITS_PER_MP_LIMB-dif) */
	    if (k!=0) cout += mpn_add_1(&cc, bp+(--k), 1,
				      cp[0]<<(BITS_PER_MP_LIMB-dif));
	    else cc = cp[0]<<(BITS_PER_MP_LIMB-dif);
	    if ((cout==0 && cc==0) || (~cout==0 && ~cc==0)) cout=cc;
	  }
	  if ((long)cout>0 || (cout==0 && cc)) goto add_one_ulp;
	  else if ((long)cout<0)
	    { TMP_FREE(marker); return; /* no carry possible any more */ }
	  else if (kc==0) {
	    while (k && cout==0) cout=bp[--k];
	    if ((~cout) && (cout || (rnd_mode==GMP_RNDN && (*ap & (ONE<<sh)))))
	      goto add_one_ulp;
	    else goto end_of_add;
	  }

	  /* else round c: go through */

	case 3: /* only c to round */
	  bp=cp; k=cn-k; bn=cn;
	  goto to_nearest;

	case 0: /* only b to round */
	  k=bn-an; dif=0;
	  goto to_nearest;

	  /* otherwise the result is exact: nothing to do */
	}
    }
    /* else round towards zero:
       - if low(b)+low(c) >= ulp(a) then add one
         (this can happen only if the last sh bits from a are all 1)
       - else truncate
     */
    else
      {
	mp_limb_t carry, mask, tp, lastc;

	/*           
	   <--------|---b----|-------->
	     <--------|---c----|--------|--------> (dif bits to the right)
	*/
	mask = (ONE << sh) - 1;
	carry = *ap & mask;
	*ap -= carry;
	if (carry == mask) /* all last sh bits from a are 1 */
	  {
	    bn -= an;
	    cn -= k;
	    carry = ~((mp_limb_t) 0);
	    lastc = (dif) ? (cp[cn] << (BITS_PER_MP_LIMB - dif)) : 0;
	    while ((bn || cn) && (~carry == 0))
	      {
		tp = (bn) ? bp[--bn] : 0; /* current limb from b */
		if (cn)
		  lastc += cp[--cn] >> dif; /* corresponding limb from c */
		carry += lastc;
		/* if carry < lastc : b[i] + c[i] >= BASE ==> add one ulp */
		if (carry < lastc)
		  goto add_one_ulp;
	      }
	  }
      }
    goto end_of_add;

  to_nearest: /* 0 <= sh < BITS_PER_MP_LIMB : number of bits of a to truncate
                 bp[k] : last significant limb from b
		 bn : number of limbs of b
	      */
        /* c3=1 whenever b+c gave a carry out in most significant limb
	   and the least significant bit (shifted right) was 1.
	   This can occur only when BITS_PER_MP_LIMB divides MPFR_PREC(a),
	   i.e. sh=0.
	 */
        if (sh) {
	  cc = *ap & ((ONE<<sh)-1);
	  *ap &= ~cc; /* truncate last bits */
	  c2 = (rnd_mode==GMP_RNDN) ? ONE<<(sh-1) : 0;
	}
	else /* sh=0: no bit to truncate */
	  {
	    if (k) cc = bp[--k]; else cc = 0;
	    c2 = (rnd_mode==GMP_RNDN) ? ONE<<(BITS_PER_MP_LIMB-1) : 0;
	    if (c3 && (cc || c2==0)) cc=c2+1; /* will force adding one ulp */
	  }
	if (cc>c2) goto add_one_ulp; /* trunc(b)>1/2*lsb(a) -> round up */
	else if (cc==c2) {
	  /* special case of rouding c shifted to the right */
	  if (dif>0 && k<bn) cc=bp[k]<<(BITS_PER_MP_LIMB-dif);
	  else cc=0;
	  while (k && cc==0) cc=bp[--k];
	  /* now if the truncated part of b = 1/2*lsb(a), check whether c=0 */
	  if (cc || (rnd_mode==GMP_RNDN && (*ap & (ONE<<sh))))
	    goto add_one_ulp;
	}
  }
  goto end_of_add;

  add_one_ulp: /* add one unit in last place to a */
    while (nulp--)
    {
      if (mpn_add_1(ap, ap, an, ONE<<sh))
      {
        mp_exp_t exp = MPFR_EXP(a);
        if (exp == __mpfr_emax)
        {
          mpfr_set_overflow(a, rnd_mode, MPFR_SIGN(a));
          break;
        }
        else
        {
          MPFR_EXP(a)++;
          ap[an-1] = MP_LIMB_T_HIGHBIT;
          if (rnd_mode == GMP_RNDN)
            break;  /* because ulp is doubled */
        }
      }
    }

 end_of_add:
  TMP_FREE(marker);
  return;
}

void
#if __STDC__
mpfr_add (mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mp_rnd_t rnd_mode)
#else
mpfr_add (a, b, c, rnd_mode)
     mpfr_ptr a;
     mpfr_srcptr b;
     mpfr_srcptr c;
     mp_rnd_t rnd_mode;
#endif
{
  if (MPFR_IS_NAN(b) || MPFR_IS_NAN(c))
  {
    MPFR_SET_NAN(a);
    return;
  }

  MPFR_CLEAR_NAN(a);

  if (MPFR_IS_INF(b))
  {
    if (!MPFR_IS_INF(c) || MPFR_SIGN(b) == MPFR_SIGN(c))
    {
      MPFR_SET_INF(a);
      MPFR_SET_SAME_SIGN(a, b);
    }
    else
      MPFR_SET_NAN(a);
    return;
  }
  else
    if (MPFR_IS_INF(c))
    {
      MPFR_SET_INF(a);
      MPFR_SET_SAME_SIGN(a, c);
      return;
    }

  MPFR_ASSERTN(MPFR_IS_FP(b) && MPFR_IS_FP(c));

  if (MPFR_IS_ZERO(b))
  {
    if (MPFR_IS_ZERO(c))
    {
      if (MPFR_SIGN(a) !=
          (rnd_mode != GMP_RNDD ?
           ((MPFR_SIGN(b) < 0 && MPFR_SIGN(c) < 0) ? -1 : 1) :
           ((MPFR_SIGN(b) > 0 && MPFR_SIGN(c) > 0) ? 1 : -1)))
        MPFR_CHANGE_SIGN(a);
      MPFR_CLEAR_INF(a);
      MPFR_SET_ZERO(a);
      return;
    }
    mpfr_set(a, c, rnd_mode);
    return;
  }

  if (MPFR_IS_ZERO(c))
  {
    mpfr_set(a, b, rnd_mode);
    return;
  }

  MPFR_CLEAR_INF(a); /* clear Inf flag */

  if (MPFR_SIGN(b) != MPFR_SIGN(c))
  { /* signs differ, it's a subtraction */
    if (MPFR_EXP(b) < MPFR_EXP(c))
    {
      mpfr_sub1(a, c, b, rnd_mode,
                (mp_exp_unsigned_t) MPFR_EXP(c) - MPFR_EXP(b));
    }
    else if (MPFR_EXP(b) > MPFR_EXP(c))
    {
      mpfr_sub1(a, b, c, rnd_mode,
                (mp_exp_unsigned_t) MPFR_EXP(b) - MPFR_EXP(c));
    }
    else
    { /* MPFR_EXP(b) == MPFR_EXP(c) */
      int d = mpfr_cmp_abs(b,c);
      if (d == 0)
      {
        if (rnd_mode == GMP_RNDD)
          MPFR_SET_NEG(a);
        else
          MPFR_SET_POS(a);
        MPFR_SET_ZERO(a);
      }
      else if (d > 0)
        mpfr_sub1(a, b, c, rnd_mode, 0);
      else
        mpfr_sub1(a, c, b, rnd_mode, 0);
    }
  }
  else
  { /* signs are equal, it's an addition */
    if (MPFR_EXP(b) < MPFR_EXP(c))
    {
      mpfr_add1(a, c, b, rnd_mode,
                (mp_exp_unsigned_t) MPFR_EXP(c) - MPFR_EXP(b));
    }
    else
    {
      mpfr_add1(a, b, c, rnd_mode,
                (mp_exp_unsigned_t) MPFR_EXP(b) - MPFR_EXP(c));
    }
  }
}
