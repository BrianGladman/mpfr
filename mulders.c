/* Mulder's MulHigh function
   
Copyright 2005 Free Software Foundation.

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
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */


#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* Don't MPFR_MULHIGH_SIZE since it is handled by tuneup */
#ifdef MPFR_MULHIGH_TAB_SIZE
static short mulhigh_ktab[MPFR_MULHIGH_TAB_SIZE];
#else
static short mulhigh_ktab[] = {MPFR_MULHIGH_TAB};
#define MPFR_MULHIGH_TAB_SIZE (sizeof(mulhigh_ktab) / sizeof(mulhigh_ktab[0]))
#endif


/* Put in  rp[n-1..2n-1] an approximation of the n+1 high limbs
   of {mp, n} * {np, n}. 
   The error is at worst of ln(n) for rp[n] and rp[n-1] is totally wrong */
static void
mpfr_mulhigh_n_basecase (mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n)
{
  mp_size_t i;

  rp += n-1;
  umul_ppmm (rp[1], rp[0], up[n-1], vp[0]);
  for (i = 1 ; i < n ; i++)
    rp[i+1] = mpn_addmul_1 (rp, up + (n - i - 1), i+1, vp[i]);
}

void
mpfr_mulhigh_n (mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
  mp_size_t k;

  k = MPFR_LIKELY (n < MPFR_MULHIGH_TAB_SIZE) ? mulhigh_ktab[n] : 2*n/3;
  if (k < 0)
    mpn_mul_basecase (rp, np, n, mp, n);
  else if (k == 0)
    mpfr_mulhigh_n_basecase (rp, np, mp, n);
  else
    {
      mp_size_t l = n - k;
      mp_limb_t cy;

      mpn_mul_n (rp + 2 * l, np + l, mp + l, k); /* fills rp[2l..2n-1] */
      mpfr_mulhigh_n (rp, np + k, mp, l);          /* fills rp[l-1..2l-1] */
      cy = mpn_add_n (rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
      mpfr_mulhigh_n (rp, np, mp + k, l);          /* fills rp[l-1..2l-1] */
      cy += mpn_add_n (rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
      mpn_add_1 (rp + n + l, rp + n + l, k, cy); /* propagate carry */
    }
}

#if 0
int mpfr_mul () 
{
#define MPFR_MUL_BASECASE_THREEHOLD 5

  /* multiplies two mantissa in temporary allocated space */
  b1 = MPFR_LIKELY (bn >= cn)
    ? mpn_mul (tmp, MPFR_MANT (b), bn, MPFR_MANT (c), cn)
    : mpn_mul (tmp, MPFR_MANT (c), cn, MPFR_MANT (b), bn);


  /* TO REPLACE WITH -->; */
  
  if (MPFR_UNLIKELY (bn < cn)) 
    {
      mpfr_ptr  tmp = b;
      mp_size_t tn  = bn;
      b = c;
      c = tmp;
      bn = cn;
      cn = tn;
    }
  if (MPFR_UNLIKELY (bn > MPFR_MUL_BASECASE_THREEHOLD))
    {
      mp_size_t cancel;
      mp_prec_t prec_cn;

      prec_cn = cn*BITS_PER_MP_LIMB-MPFR_INT_CEIL_LOG2 (cn);
      /* prec_cn is the expected precision of mulhigh */

      MPFR_ASSERTD (bn >= cn);
      /* FIXME: Find best guard bits to add */
      if (MPFR_UNLIKELY (MPFR_PREC (a) > prec_cn - 4))
	/* MulHigh can't produce a roundable result.
	   Do the full multiply instead. */
	goto full_multiply;
      cancel = 0;
      if (MPFR_UNLIKELY (MPFR_PREC (a) < prec_cn - 4 -  BITS_PER_MP_LIMB))
	{
	  /* MulHigh will computes too much bits */
	  cancel = (prec_cn - 4 - MPFR_PREC (a)) / BITS_PER_MP_LIMB;
	  MPFR_ASSERTD (cancel >= 1);
	}
      mpfr_mulhigh_n (tmp+2*cancel, MPFR_MANT (b) + cancel, 
		      MPFR_MANT (c) + cancel, cn-cancel);
      /* FIXME: tn or k? */
      if (MPFR_LIKELY (mpfr_can_round_raw (tmp, k, sign, prec_cn,
					   GMP_RNDN, GMP_RNDZ,
					   MPFR_PREC(a)+(rnd_mode==GMP_RNDN))))
	b1 = tmp[2*bn-1];
      else
	goto full_multiply;
      }
    }
  else
    {
    full_multiply:
      b1 = mpn_mul (tmp, MPFR_MANT (b), bn, MPFR_MANT (c), cn);      
    }

#endif
