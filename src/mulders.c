/* Mulders' MulHigh function (short product)

Copyright 2005, 2006, 2007, 2008, 2009, 2010 Free Software Foundation, Inc.
Contributed by the Arenaire and Caramel projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

#ifndef MUL_FFT_THRESHOLD
#define MUL_FFT_THRESHOLD 8448
#endif

/* Don't use MPFR_MULHIGH_SIZE since it is handled by tuneup */
#ifdef MPFR_MULHIGH_TAB_SIZE
static short mulhigh_ktab[MPFR_MULHIGH_TAB_SIZE];
#else
static short mulhigh_ktab[] = {MPFR_MULHIGH_TAB};
#define MPFR_MULHIGH_TAB_SIZE \
  ((mp_size_t) (sizeof(mulhigh_ktab) / sizeof(mulhigh_ktab[0])))
#endif

/* Put in  rp[n..2n-1] an approximation of the n high limbs
   of {up, n} * {vp, n}. The error is less than n-1 ulps of rp[n]. */
static void
mpfr_mulhigh_n_basecase (mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n)
{
  mp_size_t i;

  rp += n - 1;
  umul_ppmm (rp[1], rp[0], up[n-1], vp[0]); /* we neglect up[0..n-2]*vp[0],
                                               which is less than B^n */
  for (i = 1 ; i < n ; i++)
    /* here, we neglect up[0..n-i-2] * vp[i], which is less than B^n too */
    rp[i + 1] = mpn_addmul_1 (rp, up + (n - i - 1), i + 1, vp[i]);
  /* in total, we neglect less than n*B^n, i.e., n ulps of rp[n]. */
}

/* Put in  rp[n..2n-1] an approximation of the n high limbs
   of {np, n} * {mp, n}. The error is less than n-1 ulps of rp[n]. */
void
mpfr_mulhigh_n (mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
  mp_size_t k;

  MPFR_ASSERTN (MPFR_MULHIGH_TAB_SIZE >= 8); /* so that 3*(n/4) > n/2 */
  k = MPFR_LIKELY (n < MPFR_MULHIGH_TAB_SIZE) ? mulhigh_ktab[n] : 3*(n/4);
  MPFR_ASSERTD (k == -1 || k == 0 || (k > n/2 && k < n));
  if (k < 0)
    mpn_mul_basecase (rp, np, n, mp, n); /* result is exact, no error */
  else if (k == 0)
    mpfr_mulhigh_n_basecase (rp, np, mp, n); /* basecase error < n-1 ulps */
  else if (n > MUL_FFT_THRESHOLD)
    mpn_mul_n (rp, np, mp, n); /* result is exact, no error */
  else
    {
      mp_size_t l = n - k;
      mp_limb_t cy;

      mpn_mul_n (rp + 2 * l, np + l, mp + l, k); /* fills rp[2l..2n-1] */
      mpfr_mulhigh_n (rp, np + k, mp, l);        /* fills rp[l-1..2l-1] */
      cy = mpn_add_n (rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
      mpfr_mulhigh_n (rp, np, mp + k, l);        /* fills rp[l-1..2l-1] */
      cy += mpn_add_n (rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
      mpn_add_1 (rp + n + l, rp + n + l, k, cy); /* propagate carry */
      /* the neglected terms are in the two recursive calls to mpfr_mulhigh_n,
         where in each case by induction the error is at most l-1 ulps, plus
         the two overlapping products {np, l} * {mp, k} and {np, k} * {mp, l},
         which are altogether bounded by B^n, thus 1 ulp each, thus
         the total error is at most 2l ulps. Since k > n/2, l < n/2 which
         gives an error < n-1 ulps. */
    }
}

#ifdef MPFR_SQRHIGH_TAB_SIZE
static short sqrhigh_ktab[MPFR_SQRHIGH_TAB_SIZE];
#else
static short sqrhigh_ktab[] = {MPFR_SQRHIGH_TAB};
#define MPFR_SQRHIGH_TAB_SIZE (sizeof(sqrhigh_ktab) / sizeof(sqrhigh_ktab[0]))
#endif

/* Put in  rp[n..2n-1] an approximation of the n high limbs
   of {np, n}^2. The error is less than n ulps of rp[n]. */
void
mpfr_sqrhigh_n (mp_ptr rp, mp_srcptr np, mp_size_t n)
{
  mp_size_t k;

  MPFR_ASSERTN (MPFR_SQRHIGH_TAB_SIZE > 2); /* ensures k < n */
  k = MPFR_LIKELY (n < MPFR_SQRHIGH_TAB_SIZE) ? sqrhigh_ktab[n]
    : n/2 + n/50 + 1; /* +1 ensures that k > n/2 */
  MPFR_ASSERTD (k == -1 || k == 0 || (k > n/2 && k < n));
  if (k < 0)
    /* we can't use mpn_sqr_basecase here, since it requires
       n <= SQR_KARATSUBA_THRESHOLD, where SQR_KARATSUBA_THRESHOLD
       is not exported by GMP */
    mpn_sqr_n (rp, np, n);
  else if (k == 0)
    mpfr_mulhigh_n_basecase (rp, np, np, n);
  else
    {
      mp_size_t l = n - k;
      mp_limb_t cy;

      mpn_sqr_n (rp + 2 * l, np + l, k);          /* fills rp[2l..2n-1] */
      mpfr_mulhigh_n (rp, np, np + k, l);         /* fills rp[l-1..2l-1] */
      /* {rp+n-1,l+1} += 2 * {rp+l-1,l+1} */
      cy = mpn_lshift (rp + l - 1, rp + l - 1, l + 1, 1);
      cy += mpn_add_n (rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
      mpn_add_1 (rp + n + l, rp + n + l, k, cy); /* propagate carry */
    }
}

#ifdef MPFR_DIVHIGH_TAB_SIZE
static short divhigh_ktab[MPFR_DIVHIGH_TAB_SIZE];
#else
static short divhigh_ktab[] = {MPFR_DIVHIGH_TAB};
#define MPFR_DIVHIGH_TAB_SIZE (sizeof(divhigh_ktab) / sizeof(divhigh_ktab[0]))
#endif

/* put in {qp, n} an approximation of N={np, 2*n} divided by D={dp, n},
   with the most significant limb of the quotient as return value (0 or 1).
   Assumes the most significant bit of D is set. Clobbers N.
   THIS IS PRELIMINARY CODE, DO NOT USE IT.
*/
mp_limb_t
mpfr_divhigh_n (mp_ptr qp, mp_ptr np, mp_ptr dp, mp_size_t n)
{
  mp_size_t k, l;
  mp_limb_t qh, cy;
  mp_ptr tp;
  MPFR_TMP_DECL(marker);

  k = divhigh_ktab[n];
  MPFR_ASSERTD ((n+1)/2 <= k && k <= n); /* we should have n/2 <= k <= n */

  /* for k=0, we use a full division (mpn_divrem) */

  if (k == n)
    return mpn_divrem (qp, 0, np, 2 * n, dp, n);

  MPFR_TMP_MARK (marker);
  l = n - k;
  /* first divide the most significant 2k limbs from N by the most significant
     k limbs of D */
  qh = mpn_divrem (qp + l, 0, np + 2 * l, 2 * k, dp + l, k); /* exact */

  /* it remains {np,2l+k} = {np,n+l} as remainder */

  /* now we have to subtract high(Q1)*D0 where Q1={qp+l,k} and D0={dp,l} */
  tp = (mp_ptr) MPFR_TMP_ALLOC (2 * l * sizeof (mp_limb_t));
  mpfr_mulhigh_n (tp, qp + k, dp, l);
  /* we are only interested in the upper l limbs from {tp,2l} */
  cy = mpn_sub_n (np + n, np + n, tp + l, l);
  while (cy > 0) /* Q1 was too large: subtract 1 to Q1 and add D to np+l */
    {
      qh -= mpn_sub_1 (qp + l, qp + l, k, MPFR_LIMB_ONE);
      cy -= mpn_add_n (np + l, np + l, dp, n);
    }

  /* now it remains {np,n+l} to divide by D */
  cy = mpfr_divhigh_n (qp, np + k, dp + k, l);
  qh += mpn_add_1 (qp + l, qp + l, k, cy);
  MPFR_TMP_FREE(marker);

  return qh;
}
