/* Mulders' short product, square and division.

Copyright 2005-2025 Free Software Foundation, Inc.
Contributed by the Pascaline and Caramba projects, INRIA.

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
along with the GNU MPFR Library; see the file COPYING.LESSER.
If not, see <https://www.gnu.org/licenses/>. */

/* References:
   [1] Short Division of Long Integers, David Harvey and Paul Zimmermann,
       Proceedings of the 20th Symposium on Computer Arithmetic (ARITH-20),
       July 25-27, 2011, pages 7-14.
   [2] Quadratic Short Division, Juraj Sukop and Paul Zimmermann,
       preprint, https://inria.hal.science/hal-04557431, 2024.
*/

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* Don't use MPFR_MULHIGH_SIZE since it is handled by tuneup */
#ifdef MPFR_MULHIGH_TAB_SIZE
static short mulhigh_ktab[MPFR_MULHIGH_TAB_SIZE];
#else
static short mulhigh_ktab[] = {MPFR_MULHIGH_TAB};
#define MPFR_MULHIGH_TAB_SIZE (numberof_const (mulhigh_ktab))
#endif

/* Put in  rp[n..2n-1] an approximation of the n high limbs
   of {up, n} * {vp, n}. The error is less than n ulps of rp[n] (and the
   approximation is always less or equal to the truncated full product).
   Assume 2n limbs are allocated at rp.

   Implements Algorithm ShortMulNaive from [1].
*/
static void
mpfr_mulhigh_n_basecase (mpfr_limb_ptr rp, mpfr_limb_srcptr up,
                         mpfr_limb_srcptr vp, mp_size_t n)
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
   of {np, n} * {mp, n}. The error is less than n ulps of rp[n] (and the
   approximation is always less or equal to the truncated full product).

   Implements Algorithm ShortMul from [1].
*/
void
mpfr_mulhigh_n (mpfr_limb_ptr rp, mpfr_limb_srcptr np, mpfr_limb_srcptr mp,
                mp_size_t n)
{
  mp_size_t k;

  MPFR_STAT_STATIC_ASSERT (MPFR_MULHIGH_TAB_SIZE >= 8); /* so that 3*(n/4) > n/2 */
  k = MPFR_LIKELY (n < MPFR_MULHIGH_TAB_SIZE) ? mulhigh_ktab[n] : 3*(n/4);
  /* Algorithm ShortMul from [1] requires k >= (n+3)/2, which translates
     into k >= (n+4)/2 in the C language. */
  MPFR_ASSERTD (k == -1 || k == 0 || (k >= (n+4)/2 && k < n));
  if (k < 0)
    mpn_mul_basecase (rp, np, n, mp, n); /* result is exact, no error */
  else if (k == 0)
    mpfr_mulhigh_n_basecase (rp, np, mp, n); /* basecase error < n ulps */
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
    }
}

#ifdef MPFR_SQRHIGH_TAB_SIZE
static short sqrhigh_ktab[MPFR_SQRHIGH_TAB_SIZE];
#else
static short sqrhigh_ktab[] = {MPFR_SQRHIGH_TAB};
#define MPFR_SQRHIGH_TAB_SIZE (numberof_const (sqrhigh_ktab))
#endif

/* Put in  rp[n..2n-1] an approximation of the n high limbs
   of {np, n}^2. The error is less than n ulps of rp[n]. */
void
mpfr_sqrhigh_n (mpfr_limb_ptr rp, mpfr_limb_srcptr np, mp_size_t n)
{
  mp_size_t k;

  MPFR_STAT_STATIC_ASSERT (MPFR_SQRHIGH_TAB_SIZE > 2); /* ensures k < n */
  k = MPFR_LIKELY (n < MPFR_SQRHIGH_TAB_SIZE) ? sqrhigh_ktab[n]
    : (n+4)/2; /* ensures that k >= (n+3)/2 */
  MPFR_ASSERTD (k == -1 || k == 0 || (k >= (n+4)/2 && k < n));
  if (k < 0)
    /* we can't use mpn_sqr_basecase here, since it requires
       n <= SQR_KARATSUBA_THRESHOLD, where SQR_KARATSUBA_THRESHOLD
       is not exported by GMP */
    mpn_sqr (rp, np, n);
  else if (k == 0)
    mpfr_mulhigh_n_basecase (rp, np, np, n);
  else
    {
      mp_size_t l = n - k;
      mp_limb_t cy;

      mpn_sqr (rp + 2 * l, np + l, k);            /* fills rp[2l..2n-1] */
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
#define MPFR_DIVHIGH_TAB_SIZE (numberof_const (divhigh_ktab))
#endif

/* Put in Q={qp, n} an approximation of N={np, 2*n} divided by D={dp, n},
   with the most significant limb of the quotient as return value (0 or 1).
   Assumes the most significant bit of D is set. Clobbers N.

   The approximate quotient Q satisfies - 2(n-1) < N/D - Q <= 4.

   This implements Algorithm BasecaseShortDiv from [2], with a 3/2
   quotient selection (Remark 5 from [2]).

   Assumes n >= 2.

   We define B = 2^mp_bits_per_limb.
*/
static mp_limb_t
mpfr_divhigh_n_basecase (mpfr_limb_ptr qp, mpfr_limb_ptr np,
                         mpfr_limb_srcptr dp, mp_size_t n)
{
  mp_limb_t qh, d1, d0, q2, q1, q0;
  mpfr_pi1_t dinv2;

  MPFR_ASSERTD(n >= 2);

  np += n;

  if ((qh = (mpn_cmp (np, dp, n) >= 0)))
    mpn_sub_n (np, np, dp, n);

  /* now {np, n} is less than D={dp, n}, which implies np[n-1] <= dp[n-1] */

  d1 = dp[n - 1];

  /* we assumed n >= 2 */
  d0 = dp[n - 2];
  invert_pi1 (dinv2, d1, d0);
  /* dinv2.inv32 = floor ((B^3 - 1) / (d0 + d1 B)) - B */
  while (n > 0)
    {
      /* In Remark 5 of [2], we should divide by the upper 2 limbs of H_j,
         which coincide with the upper 2 limbs of the divisor, except for
         the last iteration where H_j has only one non-zero limb. */
      if (MPFR_UNLIKELY(np[n-1] > d1 || (np[n-1] == d1 && np[n-2] >= d0)))
        q2 = MPFR_LIMB_MAX;
      else if (n > 1)
        udiv_qr_3by2 (q2, q1, q0, np[n - 1], np[n - 2], np[n - 3],
                      d1, d0, dinv2.inv32);
      else
        { /* Case n=1: since dinv2 was computed for d0 = dp[n - 2],
                we can't use it for d0 = 0 */
          if (MPFR_UNLIKELY(np[n-1] == d1))
            q2 = MPFR_LIMB_MAX;
          else /* now np[n-1] < d1 */
            udiv_qrnnd (q2, q1, np[n - 1], np[n - 2], d1);
        }

      q0 = mpn_submul_1 (np - 1, dp, n, q2);
      /* q0 is the upper limb of the product {dp, n} * q2, minus the
         potential borrow-out from the subtraction {np-1, n} - {dp, n} * q2,
         thus it is the value we have to subtract (in theory) to np[n-1].
         We want q0 = np[n-1] to make np[n-1] vanish, so that we can
         proceed to the next limb np[n-2].
         If q0 > np[n-1], the partial quotient q2 was too large.
         If q0 < np[n-1], the partial quotient q2 was too small. */
      if (MPFR_UNLIKELY(q0 > np[n - 1])) /* q2 was too large */
        {
          q0 -= mpn_add_n (np - 1, np - 1, dp, n);
          q2 --;
        }
      if (MPFR_UNLIKELY(q0 < np[n - 1])) /* q2 was too small */
        {
          /* This implements the "early exit" of Algorithm BasecaseShortDiv
             from [2] (step 10). */
          while (n != 0)
            qp[--n] = MPFR_LIMB_MAX;
          break;
        }
      MPFR_ASSERTN(q0 == np[n - 1]);
      qp[--n] = q2;
      dp ++;
    }

  return qh;
}

/* Put in {qp, n} an approximation of N={np, 2*n} divided by D={dp, n},
   with the most significant limb of the quotient as return value (0 or 1).
   Assumes the most significant bit of D is set. Clobbers N.

   This implements the ShortDiv algorithm from reference [1].

   Assumes n >= 2 (which should be fulfilled also in the recursive calls).
*/
mp_limb_t
mpfr_divhigh_n (mpfr_limb_ptr qp, mpfr_limb_ptr np, mpfr_limb_ptr dp,
                mp_size_t n)
{
  mp_size_t k, l;
  mp_limb_t qh, cy;
  mpfr_limb_ptr tp;
  MPFR_TMP_DECL(marker);

  MPFR_LOG_FUNC
    (("n=%Pd", (mpfr_prec_t) n),
     ("k=%Pd qh=%Mu", (mpfr_prec_t) k, k == 0 ? MPFR_LIMB_ZERO : qh));

  MPFR_STAT_STATIC_ASSERT (MPFR_DIVHIGH_TAB_SIZE >= 15); /* so that 2*(n/3) >= (n+4)/2 */
  MPFR_ASSERTD(n >= 2);
  k = MPFR_LIKELY (n < MPFR_DIVHIGH_TAB_SIZE) ? divhigh_ktab[n] : 2*(n/3);

  if (k == 0)
    {
#if defined(WANT_GMP_INTERNALS) && defined(HAVE___GMPN_SBPI1_DIVAPPR_Q)
      mpfr_pi1_t dinv2;
      invert_pi1 (dinv2, dp[n - 1], dp[n - 2]);
      if (n > 2) /* sbpi1_divappr_q wants n > 2 */
        return __gmpn_sbpi1_divappr_q (qp, np, n + n, dp, n, dinv2.inv32);
      else
        return mpfr_divhigh_n_basecase (qp, np, dp, n);
#else /* use our own code for base-case short division */
      return mpfr_divhigh_n_basecase (qp, np, dp, n);
#endif
    }

  /* Check the bounds from [1]. In addition, we forbid k=n-1, which would
     give l=1 in the recursive call. It follows n >= 5. */
  MPFR_ASSERTD ((n+4)/2 <= k && k < n-1);

  MPFR_TMP_MARK (marker);
  l = n - k;
  /* first divide the most significant 2k limbs from N by the most significant
     k limbs of D */
  qh = mpn_divrem (qp + l, 0, np + 2 * l, 2 * k, dp + l, k); /* exact */

  /* it remains {np,2l+k} = {np,n+l} as remainder */

  /* now we have to subtract high(Q1)*D0 where Q1=qh*B^k+{qp+l,k} and
     D0={dp,l} */
  tp = MPFR_TMP_LIMBS_ALLOC (2 * l);
  mpfr_mulhigh_n (tp, qp + k, dp, l);
  /* we are only interested in the upper l limbs from {tp,2l} */
  cy = mpn_sub_n (np + n, np + n, tp + l, l);
  if (qh)
    cy += mpn_sub_n (np + n, np + n, dp, l);
  while (cy > 0) /* Q1 was too large: subtract 1 from Q1 and add D to np+l */
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
