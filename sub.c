/* mpfr_sub -- subtract two floating-point numbers

Copyright (C) 2001 Free Software Foundation.
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

/* #define DEBUG */

#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"

#define ONE ((mp_limb_t) 1)

extern void mpfr_add1 _PROTO((mpfr_ptr, mpfr_srcptr, mpfr_srcptr,
			      mp_rnd_t, mp_exp_unsigned_t));
int mpfr_sub1 _PROTO ((mpfr_ptr, mpfr_srcptr, mpfr_srcptr,
                        mp_rnd_t, mp_exp_unsigned_t));

/* signs of b and c differ, abs(b) > abs(c), 
   diff_exp = EXP(b) - EXP(c).
   Returns 0 iff result is exact,
   a negative value when the result is less than the exact value,
   a positive value otherwise.
*/
int
#if __STDC__
mpfr_sub1 (mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c,
	   mp_rnd_t rnd_mode, mp_exp_unsigned_t diff_exp)
#else
mpfr_sub1 (a, b, c, rnd_mode, diff_exp)
     mpfr_ptr a;
     mpfr_srcptr b;
     mpfr_srcptr c;
     mp_rnd_t rnd_mode;
     mp_exp_unsigned_t diff_exp;
#endif
{
  unsigned long cancel, cancel1, sh, k;
  long int cancel2, an, bn, cn, cn0;
  mp_limb_t *ap, *bp, *cp, carry, bb, cc, borrow = 0;
  int inexact = 0, shift_b, shift_c, maybe_exact = 0, down = 0;
  TMP_DECL(marker);

#ifdef DEBUG
  printf("\nenter mpfr_sub, rnd_mode=%s:\n", mpfr_print_rnd_mode(rnd_mode));
  printf("b="); if (MPFR_SIGN(b)>0) putchar(' '); mpfr_print_raw(b); putchar('\n');
  printf("c="); if (MPFR_SIGN(c)>0) putchar(' '); for (k=0; k<diff_exp; k++) putchar(' '); mpfr_print_raw(c); putchar('\n');
  printf("PREC(a)=%u PREC(b)=%u PREC(c)=%u\n", MPFR_PREC(a), MPFR_PREC(b),
	 MPFR_PREC(c));
#endif
  TMP_MARK(marker);
  ap = MPFR_MANT(a);
  an = 1 + (MPFR_PREC(a) - 1) / BITS_PER_MP_LIMB;

  cancel = mpfr_cmp2 (b, c);

  /* reserve a space to store b aligned with the result, i.e. shifted by
     (-cancel) % BITS_PER_MP_LIMB to the right */
  bn = 1 + (MPFR_PREC(b) - 1) / BITS_PER_MP_LIMB;
  shift_b = cancel % BITS_PER_MP_LIMB;
  if (shift_b)
    shift_b = BITS_PER_MP_LIMB - shift_b;
  cancel1 = (cancel + shift_b) / BITS_PER_MP_LIMB;
  /* the high cancel1 limbs from b should not be taken into account */
  if (shift_b == 0)
    bp = MPFR_MANT(b); /* no need of an extra space */
  else
    {
      bp = TMP_ALLOC ((bn + 1) * BYTES_PER_MP_LIMB);
      bp[0] = mpn_rshift (bp + 1, MPFR_MANT(b), bn++, shift_b);
    }

  /* reserve a space to store c aligned with the result, i.e. shifted by
     (diff_exp-cancel) % BITS_PER_MP_LIMB to the right */
  cn = 1 + (MPFR_PREC(c) - 1) / BITS_PER_MP_LIMB;
  shift_c = diff_exp - (cancel % BITS_PER_MP_LIMB);
  shift_c = (shift_c + BITS_PER_MP_LIMB) % BITS_PER_MP_LIMB;
  if (shift_c == 0)
    cp = MPFR_MANT(c);
  else
    {
      cp = TMP_ALLOC ((cn + 1) * BYTES_PER_MP_LIMB);
      cp[0] = mpn_rshift (cp + 1, MPFR_MANT(c), cn++, shift_c);
    }

#ifdef DEBUG
  printf("shift_b=%u shift_c=%u\n", shift_b, shift_c);
#endif

  /* ensure ap != bp and ap != cp */
  if (ap == bp)
    {
      bp = (mp_ptr) TMP_ALLOC(bn * BYTES_PER_MP_LIMB);
      MPN_COPY (bp, ap, bn);
      /* ap == cp cannot occur since we would have b=c, which is detected
	 in mpfr_add or mpfr_sub */
    }
  else if (ap == cp)
    {
      cp = (mp_ptr) TMP_ALLOC (cn * BYTES_PER_MP_LIMB);
      MPN_COPY(cp, ap, cn);
    }

  cancel2 = (long int) (cancel + shift_c - diff_exp) / BITS_PER_MP_LIMB;
  /* the high cancel2 limbs from b should not be taken into account */
#ifdef DEBUG
  printf("cancel=%u cancel1=%u cancel2=%d\n", cancel, cancel1, cancel2);
#endif

  /* adjust exponent and sign of result */
  MPFR_EXP(a) = MPFR_EXP(b) - cancel;
  if (MPFR_SIGN(a)*MPFR_SIGN(b) < 0)
    MPFR_CHANGE_SIGN(a);

  /*               ap[an-1]        ap[0]
             <----------------+-----------|---->
             <----------PREC(a)----------><-sh->
 cancel1
 limbs        bp[bn-cancel1-1]
 <--...-----><----------------+-----------+----------->
  cancel2
  limbs       cp[cn-cancel2-1]                                    cancel2 >= 0
    <--...--><----------------+----------------+---------------->
                (-cancel2)                                        cancel2 < 0
                   limbs      <----------------+---------------->
  */

  /* first part: put in ap[0..an-1] the value of high(b) - high(c),
     where high(b) consists of the high an+cancel1 limbs of b,
     and high(c) consists of the high an+cancel2 limbs of c.
   */

  /* copy high(b) into a */
  if (an + cancel1 <= bn) /* a: <----------------+-----------|---->
		         b: <-----------------------------------------> */
      MPN_COPY (ap, bp + bn - (an + cancel1), an);
  else  /* a: <----------------+-----------|---->
       b: <-------------------------> */
    if (cancel1 < bn) /* otherwise b does not overlap with a */
      {
	MPN_ZERO (ap, an + cancel1 - bn);
	MPN_COPY (ap + an + cancel1 - bn, bp, bn - cancel1);
      }
    else
      MPN_ZERO (ap, an);

#ifdef DEBUG
  printf("after copying high(b), a="); mpfr_print_raw(a); putchar('\n');
#endif

  /* subtract high(c) */
  if (an + cancel2 > 0) /* otherwise c does not overlap with a */
    {
      mp_limb_t *ap2;

      if (cancel2 >= 0)
	{
	  if (an + cancel2 <= cn) /* a: <----------------------------->
			      c: <-----------------------------------------> */
	    mpn_sub_n (ap, ap, cp + cn - (an + cancel2), an);
	  else /* a: <---------------------------->
	      c: <-------------------------> */
	    {
	      ap2 = ap + an + cancel2 - cn;
	      if (cn > cancel2)
		mpn_sub_n (ap2, ap2, cp, cn - cancel2);
	    }
	}
      else /* cancel2 < 0 */
	{
	  if (an + cancel2 <= cn) /* a: <----------------------------->
			                  c: <-----------------------------> */
	      borrow = mpn_sub_n (ap, ap, cp + cn - (an + cancel2), an + cancel2);
	  else /* a: <---------------------------->
	                c: <----------------> */
	    {
	      ap2 = ap + an + cancel2 - cn;
	      borrow = mpn_sub_n (ap2, ap2, cp, cn);
	    }
	  ap2 = ap + an + cancel2;
	  mpn_sub_1 (ap2, ap2, -cancel2, borrow);
	}
    }

#ifdef DEBUG
  printf("after subtracting high(c), a="); mpfr_print_raw(a); putchar('\n');
#endif

  /* now perform rounding */
  sh = an * BITS_PER_MP_LIMB - MPFR_PREC(a); /* last unused bits from a */
  carry = ap[0] & ((ONE << sh) - 1);
  ap[0] -= carry;

  if (rnd_mode == GMP_RNDN)
    {
      maybe_exact = (sh == 0) || (carry == 0);
      if (sh)
	{
	  /* can decide except when carry = 2^(sh-1) [middle]
	     or carry = 0 [truncate, but cannot decide inexact flag] */
	  if (carry > (ONE << (sh - 1)))
	    goto add_one_ulp;
	  else if ((0 < carry) && (carry < (ONE << (sh - 1))))
	    {
	      inexact = -1; /* result if smaller than exact value */
	      goto truncate;
	    }
	}
    }
  else /* directed rounding: set rnd_mode to RNDZ iff towards zero */
    {
      if (((rnd_mode == GMP_RNDD) && (MPFR_SIGN(b) > 0)) ||
	  ((rnd_mode == GMP_RNDU) && (MPFR_SIGN(b) < 0)))
	rnd_mode = GMP_RNDZ;

      if (carry)
	{
	  if (rnd_mode == GMP_RNDZ)
	    {
	      inexact = -1;
	      goto truncate;
	    }
	  else /* round away */
	    goto add_one_ulp;
	}
    }

  /* we have to consider the low (bn - (an+cancel1)) limbs from b,
     and the (cn - (an+cancel2)) limbs from c. */
  bn -= an + cancel1;
  cn0 = cn;
  cn -= (long int) an + cancel2;
#ifdef DEBUG
  printf("last %u bits from a are %lu, bn=%ld, cn=%ld\n", sh, carry, bn, cn);
#endif

  for (k=0; (bn > 0) || (cn > 0); k++)
    {
      bb = (bn > 0) ? bp[--bn] : 0;
      if ((cn > 0) && (cn-- <= cn0))
	cc = cp[cn];
      else
	cc = 0;

#ifdef DEBUG
      printf("k=%u bb=%lu cc=%lu\n", k, bb, cc);
#endif
      if ((rnd_mode == GMP_RNDN) && !k && !sh && !(maybe_exact = (bb == cc)))
	{
	  mp_limb_t half = ONE << (BITS_PER_MP_LIMB - 1);

	  /* add one ulp if bb > cc + half
	     truncate if cc - half < bb < cc + half
	     sub one ulp if bb < cc - half
	  */
	  if ((down = (bb < cc)))
	    {
	      if (cc >= half)
		cc -= half;
	      else
		bb += half;
	    }
	  else /* bb > cc */
	    {
	      if (cc < half)
		cc += half;
	      else
		bb -= half;
	    }
	}

#ifdef DEBUG
      printf("    bb=%lu cc=%lu\n", bb, cc);
#endif
      if (bb < cc)
	{
	  if (rnd_mode == GMP_RNDZ)
	    goto sub_one_ulp;
	  else if (rnd_mode != GMP_RNDN) /* round away */
	    {
	      inexact = 1;
	      goto truncate;
	    }
	  else /* round to nearest */
	    {
	      if (maybe_exact)
		{
		  inexact = 1;
		  goto truncate;
		}
	      else if (down)
		goto sub_one_ulp;
	      else
		{
		  inexact = -1;
		  goto truncate;
		}
	    }
	}
      else if (bb > cc)
	{
	  if (rnd_mode == GMP_RNDZ)
	    {
	      inexact = -1;
	      goto truncate;
	    }
	  else if (rnd_mode != GMP_RNDN) /* round away */
	      goto add_one_ulp;
	  else /* round to nearest */
	    {
	      if (maybe_exact)
		{
		  inexact = -1;
		  goto truncate;
		}
	      else if (down)
		{
		  inexact = 1;
		  goto truncate;
		}
	      else
		goto add_one_ulp;
	    }
	}
    }

  if ((rnd_mode == GMP_RNDN) && !maybe_exact)
    {
      /* even rounding rule */
      if ((ap[0] >> sh) & 1)
	goto add_one_ulp;
      else
	inexact = -1;
    }
  else
    inexact = 0;
  goto truncate;
  
 sub_one_ulp: /* add one unit in last place to a */
  mpn_sub_1 (ap, ap, an, ONE << sh);
  inexact = -1;
  goto end_of_sub;

 add_one_ulp: /* add one unit in last place to a */
  if (mpn_add_1 (ap, ap, an, ONE << sh)) /* result is a power of two */
    {
      ap[an-1] |= ONE << (BITS_PER_MP_LIMB - 1);
      MPFR_EXP(a)++;
    }
  inexact = 1; /* result larger than exact value */

 truncate:
  if ((ap[an-1] >> (BITS_PER_MP_LIMB - 1)) == 0) /* case 1 - epsilon */
    {
      ap[an-1] = ONE << (BITS_PER_MP_LIMB - 1);
      MPFR_EXP(a) ++;
    }

 end_of_sub:
  TMP_FREE(marker);
#ifdef DEBUG
  printf ("result is a="); mpfr_print_raw(a); putchar('\n');
#endif
  /* check that result is msb-normalized */
  ASSERT_ALWAYS(ap[an-1] > ~ap[an-1]);
  return inexact * MPFR_SIGN(b);
}

int
#if __STDC__
mpfr_sub (mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mp_rnd_t rnd_mode)
#else
mpfr_sub (a, b, c, rnd_mode)
     mpfr_ptr a;
     mpfr_srcptr b;
     mpfr_srcptr c;
     mp_rnd_t rnd_mode;
#endif
{
  int inexact;

  if (MPFR_IS_NAN(b) || MPFR_IS_NAN(c))
  {
    MPFR_SET_NAN(a);
    return 1; /* a NaN result is inexact */
  }

  MPFR_CLEAR_NAN(a);

  if (MPFR_IS_INF(b))
  {
    if (!MPFR_IS_INF(c) || MPFR_SIGN(b) != MPFR_SIGN(c))
    {
      MPFR_SET_INF(a);
      MPFR_SET_SAME_SIGN(a, b);
      return 0; /* +/-infinity is exact */
    }
    else
      {
	MPFR_SET_NAN(a);
	return 1; /* a NaN result is inexact */
      }
  }
  else
    if (MPFR_IS_INF(c))
    {
      MPFR_SET_INF(a);
      if (MPFR_SIGN(c) == MPFR_SIGN(a)) 
	MPFR_CHANGE_SIGN(a);
      return 0; /* +/-infinity is exact */
    }

  MPFR_ASSERTN(MPFR_IS_FP(b) && MPFR_IS_FP(c));

  if (MPFR_IS_ZERO(b))
  {
    if (MPFR_IS_ZERO(c))
    {
      if (MPFR_SIGN(a) !=
          (rnd_mode != GMP_RNDD ?
           ((MPFR_SIGN(b) < 0 && MPFR_SIGN(c) > 0) ? -1 : 1) :
           ((MPFR_SIGN(b) > 0 && MPFR_SIGN(c) < 0) ? 1 : -1)))
        MPFR_CHANGE_SIGN(a);
      MPFR_CLEAR_INF(a);
      MPFR_SET_ZERO(a);
      return 0; /* 0 - 0 is exact */
    }
    mpfr_neg (a, c, rnd_mode);
    return 0; /* 0 - c is exact */
  }

  if (MPFR_IS_ZERO(c))
  {
    mpfr_set (a, b, rnd_mode);
    return 0; /* b - 0 is exact */
  }

  MPFR_CLEAR_INF(a);

  if (MPFR_SIGN(b) == MPFR_SIGN(c))
  { /* signs are equal, it's a real subtraction */
    if (MPFR_EXP(b) < MPFR_EXP(c))
    { /* exchange rounding modes towards +/- infinity */
      if (rnd_mode == GMP_RNDU)
        rnd_mode = GMP_RNDD;
      else if (rnd_mode == GMP_RNDD)
        rnd_mode = GMP_RNDU;
      inexact = -mpfr_sub1(a, c, b, rnd_mode,
			   (mp_exp_unsigned_t) MPFR_EXP(c) - MPFR_EXP(b));
      MPFR_CHANGE_SIGN(a);
    }
    else if (MPFR_EXP(b) > MPFR_EXP(c))
      inexact = mpfr_sub1(a, b, c, rnd_mode,
			  (mp_exp_unsigned_t) MPFR_EXP(b) - MPFR_EXP(c));
    else
    { /* MPFR_EXP(b) == MPFR_EXP(c) */
      int d = mpfr_cmp_abs (b, c);

      if (d == 0)
	{
	  if (rnd_mode == GMP_RNDD)
	    MPFR_SET_NEG(a);
	  else
	    MPFR_SET_POS(a);
	  MPFR_SET_ZERO(a);
	  inexact = 0;
      }
      else if (d > 0)
        inexact = mpfr_sub1 (a, b, c, rnd_mode, 0);
      else
      { /* exchange rounding modes towards +/- infinity */
        if (rnd_mode == GMP_RNDU)
          rnd_mode = GMP_RNDD;
        else if (rnd_mode == GMP_RNDD)
          rnd_mode = GMP_RNDU;
	inexact = -mpfr_sub1 (a, c, b, rnd_mode, 0);
	MPFR_CHANGE_SIGN(a);
      }
    }
  }
  else
  { /* signs differ, it's an addition */
    if (MPFR_EXP(b) < MPFR_EXP(c))
    { /* exchange rounding modes towards +/- infinity */
      if (rnd_mode == GMP_RNDU)
        rnd_mode = GMP_RNDD;
      else if (rnd_mode == GMP_RNDD)
        rnd_mode = GMP_RNDU;
      mpfr_add1(a, c, b, rnd_mode,
                (mp_exp_unsigned_t) MPFR_EXP(c) - MPFR_EXP(b));
      MPFR_CHANGE_SIGN(a);
    }
    else
    {
      mpfr_add1(a, b, c, rnd_mode,
                (mp_exp_unsigned_t) MPFR_EXP(b) - MPFR_EXP(c));
    }
  }
  return inexact;
}
