/* This file is a new (faster) implementation of mpn_rootrem.

   Copyright 2005 Torbjo"rn Granlund, Paul Zimmermann.
*/

#include <stdio.h>		/* for printf */
#include <stdlib.h>		/* for exit, strtoul. */

#include "gmp.h"

/* replacement for internal gmp macros/functions */
#define ASSERT_ALWAYS(x)
#ifndef TMP_ALLOC_LIMBS
#define TMP_ALLOC_LIMBS(x) TMP_ALLOC(x * BYTES_PER_MP_LIMB)
#endif
#ifndef mpn_incr_u
#define mpn_incr_u(p,incr)           \
  do {                               \
    mp_limb_t cy = incr;             \
    mp_ptr q = p;                    \
    do {                             \
       cy = mpn_add_1 (q, q, 1, cy); \
       q ++;                         \
    } while (cy != 0);                 \
  } while (0)
#endif
#ifndef mpn_decr_u
#define mpn_decr_u(p,incr)           \
  do {                               \
    mp_limb_t cy = incr;             \
    mp_ptr q = p;                    \
    do {                             \
       cy = mpn_sub_1 (q, q, 1, cy); \
       q ++;                         \
    } while (cy != 0);                 \
  } while (0)
#endif

/* #define DEBUG */

#define MPN_RSHIFT(cy,rp,up,un,cnt) \
  do {									\
    if ((cnt) != 0)							\
      cy = mpn_rshift (rp, up, un, cnt);				\
    else								\
      {									\
	MPN_COPY_INCR (rp, up, un);					\
	cy = 0;								\
      }									\
  } while (0)

#define MPN_LSHIFT(cy,rp,up,un,cnt) \
  do {									\
    if ((cnt) != 0)							\
      cy = mpn_lshift (rp, up, un, cnt);				\
    else								\
      {									\
	MPN_COPY_DECR (rp, up, un);					\
	cy = 0;								\
      }									\
  } while (0)

#define MAXOFF 100
int off[MAXOFF];

static mp_size_t
mpfr_mpn_rootrem (mp_ptr rootp, mp_ptr remp,
	     mp_srcptr up, mp_size_t un, mp_limb_t k)
{
  mp_ptr pp, qp, rp, sp, tp, wp;
  mp_size_t pn, qn, rn, sn, tn, wn, nl, bn;
  mp_limb_t save, save2, cy;
  unsigned long int unb, xnb;
  unsigned int cnt;
  unsigned long b, kk;
  unsigned long sizes[GMP_NUMB_BITS];
  int ni, i;
  int c;
  long B;
  TMP_DECL (marker);

  ASSERT_ALWAYS (un > 0);
  ASSERT_ALWAYS (up[un - 1] != 0);
  ASSERT_ALWAYS (k > 1);

  TMP_MARK (marker);

#define PP_ALLOC (2 + (mp_size_t) (un*1.585))
  pp = TMP_ALLOC_LIMBS (PP_ALLOC);
  qp = TMP_ALLOC_LIMBS (PP_ALLOC);
  rp = TMP_ALLOC_LIMBS (PP_ALLOC);
  sp = TMP_ALLOC_LIMBS (PP_ALLOC);
  wp = TMP_ALLOC_LIMBS (PP_ALLOC);

  count_leading_zeros (cnt, up[un - 1]);
  unb = un * GMP_NUMB_BITS - cnt + GMP_NAIL_BITS;

  xnb = (unb - 1) / k + 1;	/* ceil (unb / k) */
#ifdef DEBUG
  printf ("unb=%u xnb=%u\n", unb, xnb);
#endif

  if (xnb == 1)
    {
      if (remp == NULL)
	remp = pp;
      mpn_sub_1 (remp, up, un, (mp_limb_t) 1);
      MPN_NORMALIZE (remp, un);	/* There should be at most one zero limb,
				   if we demand u to be normalized  */
      rootp[0] = 1;
      TMP_FREE (marker);
      return un;
    }

  tp = TMP_ALLOC_LIMBS (PP_ALLOC);
  kk = k * (xnb - 1);
  MPN_RSHIFT (cy, rp, up + kk / GMP_NUMB_BITS, un - kk / GMP_NUMB_BITS, kk % GMP_NUMB_BITS);
  rn = un - kk / GMP_NUMB_BITS;
  mpn_sub_1 (rp, rp, rn, 1);
  sp[0] = 1;
  sn = 1;
  b = 1;
  B = 1;

  b = xnb - 1;
  ni = 0;
  while (b != 1)
    {
      sizes[ni] = b;
      b = (b + 1) / 2;
      ni ++;
    }
  sizes[ni++] = 1;
  sizes[ni] = 0;

  for (i = ni; i != 0; i--)
    {
      b = sizes[i - 1] - sizes[i];

      /* Reinsert a low zero limb if we normalized away the entire remainder. */
      if (rn == 0)
	{
	  rp[0] = 0;
	  rn = 1;
	}

      MPN_LSHIFT (cy, rp + b / GMP_NUMB_BITS, rp, rn, b % GMP_NUMB_BITS);
      rn = rn + b / GMP_NUMB_BITS;
      rp[rn] = cy;
      rn += cy != 0;
      kk = kk - b;

      /* Now insert more bits from U.  The lowest bit is kk, then get b higher bits. */
      bn = b / GMP_NUMB_BITS; /* lowest limb from high part of rp[] */
      save = rp[bn];
      /* we have to save rp[bn] up to rp[nl-1], i.e. 1 or 2 limbs */
      nl = 1 + (kk + b - 1) / GMP_NUMB_BITS - (kk / GMP_NUMB_BITS);
      if (nl - 1 > bn)
	save2 = rp[bn + 1];
      MPN_RSHIFT (cy, rp, up + kk / GMP_NUMB_BITS, nl, kk % GMP_NUMB_BITS);
      rp[bn] &= ((mp_limb_t) 1 << (b % GMP_NUMB_BITS)) - 1;
      rp[bn] |= save;
      if (nl - 1 > bn)
	rp[bn + 1] = save2; /* the low b bits go in rp[0..bn] only */

      if (i == ni)
	wn = mpn_pow_1 (wp, sp, sn, k - 1, tp);
      cy = mpn_mul_1 (pp, wp, wn, k);
      pn = wn;
      pp[pn] = cy;
      pn += cy != 0;

      if (rn < pn)
	{
	  qn = 0;
	}
      else
	{
	  mp_size_t cut;
	  qn = rn - pn;
	  cut = (pn < qn + 2) ? 0 : pn - (qn + 2);
	  /*	  printf ("%5ld %5ld %5ld %5ld\n", rn, pn, qn, cut); */
	  mpn_tdiv_qr (qp, tp, 0, rp + cut, rn - cut, pp + cut, pn - cut);
	  qn += qp[qn] != 0;
	}

      bn = (b - 1) / GMP_NUMB_BITS + 1; /* number of limbs used by b bits,
					   when least significant bit is
					   aligned to least limb */

      /* q should be smaller than 2^b */
      if (qn > bn || (qn == bn && (b % GMP_NUMB_BITS != 0) &&
		      qp[qn - 1] >= ((mp_limb_t) 1 << (b % GMP_NUMB_BITS))))
	{
	  qn = b / GMP_NUMB_BITS + 1; /* b+1 bits */
	  MPN_ZERO (qp, qn);
	  qp[qn - 1] = (mp_limb_t) 1 << (b % GMP_NUMB_BITS);
	  mpn_sub_1 (qp, qp, qn, (mp_limb_t) 1);
	  qn -= qp[qn - 1] == 0;
	}

      MPN_LSHIFT (cy, sp + b / GMP_NUMB_BITS, sp, sn, b % GMP_NUMB_BITS);
      sn = sn + b / GMP_NUMB_BITS;
      sp[sn] = cy;
      sn += cy != 0;

      ASSERT_ALWAYS (bn >= qn);

      MPN_ZERO (qp + qn, bn - qn);
      /* Combine sB and q to form sB + q.  */
      save = sp[b / GMP_NUMB_BITS];
      MPN_COPY (sp, qp, bn);
      if (b % GMP_NUMB_BITS) /* if b is a multiple of GMP_NUMB_BITS,
				there is no overlap */
	mpn_incr_u (sp + b / GMP_NUMB_BITS, save);

      kk -= (k - 1) * b;
      MPN_RSHIFT (cy, tp, up + kk / GMP_NUMB_BITS, un - kk / GMP_NUMB_BITS, kk % GMP_NUMB_BITS);
      tn = un - kk / GMP_NUMB_BITS;
      tn -= tp[tn - 1] == 0;

      for (c = 0;; c++)
	{
	  /* Compute T.  */
	  if (i == 1)
	    {
	      /* Last iteration */
	      pn = mpn_pow_1 (pp, sp, sn, k, qp);
	    }
	  else
	    {
	      wn = mpn_pow_1 (wp, sp, sn, k - 1, qp);
	      mpn_mul (pp, wp, wn, sp, sn);
	      pn = wn + sn;
	      pn -= pp[pn - 1] == 0;
	    }

	  if (pn > tn || (pn == tn && mpn_cmp (pp, tp, pn) > 0))
	    mpn_decr_u (sp, 1);
	  else
	    break;
	}
      if (c < MAXOFF)
	off[c]++;
#ifdef DEBUG
      if (c > 1)
	{
	  printf ("b =%5lu  B =%5ld  k =%6ld  c =%5d  un =%7ld\n", b, B, k, c, un);
	}
#endif
      B += b;

      ASSERT_ALWAYS (tn >= pn);
      mpn_sub (rp, tp, tn, pp, pn);
      rn = tn;
      MPN_NORMALIZE (rp, rn);
    }
  MPN_COPY (rootp, sp, sn);
  if (remp != NULL)
    MPN_COPY (remp, rp, rn);
  return rn;
}


