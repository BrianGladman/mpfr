/* mpfr_sqrt -- square root of a floating-point number

Copyright 1999-2016 Free Software Foundation, Inc.
Contributed by the AriC and Caramba projects, INRIA.

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

#if !defined(MPFR_GENERIC_ABI) && GMP_NUMB_BITS == 64

/* The Taylor coefficient of order 0 of sqrt(i/2^8+x) is
   U[i-64][0]/2^64 + U[i-64][1]/2^128, then the Taylor coefficient of order j is
   (up to sign) U[i-64][j+1]/2^(64-8*j).
   The maximal number of bits is:
   j=1:64 j=2:56 j=3:49 j=4:43 j=5:36 j=6:30 j=7:23
   The sign is implicit: u[j] < 0 for j even except j=0.
   The maximal error is < .927e-21 (attained for i=64). */

#include "sqrt_tab.h"

/* Return an approximation of sqrt(2^64*n), with 2^62 <= n < 2^64,
   and error < 1 ulp (in unknown direction).
   We use a Taylor polynomial of degree 7. */
static mp_limb_t
mpn_sqrtrem2_approx (mp_limb_t n)
{
  int i = n >> 56;
  mp_limb_t x, h, l;
  const mp_limb_t *u;

  x = n << 8;
  u = U[i - 64];
  umul_ppmm (h, l, u[8], x);
  /* the truncation error on h is at most 1 here */
  umul_ppmm (h, l, u[7] - h, x);
  /* the truncation error on h is at most 2 */
  umul_ppmm (h, l, u[6] - h, x);
  /* the truncation error on h is at most 3 */
  umul_ppmm (h, l, u[5] - h, x);
  /* the truncation error on h is at most 4 */
  umul_ppmm (h, l, u[4] - h, x);
  /* the truncation error on h is at most 5 */
  umul_ppmm (h, l, u[3] - h, x);
  /* the truncation error on h is at most 6 */
  umul_ppmm (h, l, u[2] - h, x >> 8); /* here we shift by 8 since u[0] has weight
                                         1/2^64 and u[2] has weight 1/2^72, the
                                         truncation error on h+l/2^64 is <= 6/2^8 */
  add_ssaaaa (h, l, h, l, u[0], u[1]);
  /* Since the above addition is exact, the truncation error on h + l/2^64
     is still at most 6/2^8. Together with the mathematical error < .927e-21*2^64,
     the total error on h + l/2^64 is < 0.0406 */
  return h + (l >> 63); /* round to nearest */
}

/* put in rh,rl the upper 2 limbs of the product xh,xl * yh,yl,
   with error less than 3 ulps */
#define umul_ppmm2(rh,rl,xh,xl,yh,yl)    \
  do                                     \
    {                                    \
      mp_limb_t _h, _l;                  \
      umul_ppmm (rh, rl, xh, yh);        \
      umul_ppmm (_h, _l, xh, yl);        \
      rl += _h;                          \
      rh += (rl < _h);                   \
      umul_ppmm (_h, _l, xl, yh);        \
      rl += _h;                          \
      rh += (rl < _h);                   \
    }                                    \
  while (0)

/* Put in rp[1]*2^64+rp[0] an approximation of sqrt(2^128*n),
   with 2^126 <= n := np[1]*2^64 + np[0] < 2^128.
   The error on {rp, 2} is less than 43 ulps (in unknown direction).
*/
static void
mpn_sqrtrem4_approx (mpfr_limb_ptr rp, mpfr_limb_srcptr np)
{
  int i = np[1] >> 56;
  mp_limb_t x, r1, r0, h, l, t;
  const mp_limb_t *u;
  const mp_limb_t magic = 0xda9fbe76c8b43800; /* ceil(0.854*2^64) */

  x = np[1] << 8 | (np[0] >> 56);
  u = V[i - 64];
  umul_ppmm (h, l, u[8], x);
  /* the truncation error on h is at most 1 here */
  umul_ppmm (h, l, u[7] - h, x);
  /* the truncation error on h is at most 2 */
  umul_ppmm (h, l, u[6] - h, x);
  /* the truncation error on h is at most 3 */
  umul_ppmm (h, l, u[5] - h, x);
  /* the truncation error on h is at most 4 */
  umul_ppmm (h, l, u[4] - h, x);
  /* the truncation error on h is at most 5 */
  umul_ppmm (h, l, u[3] - h, x);
  /* the truncation error on h is at most 6 */
  umul_ppmm (h, l, u[2] - h, x >> 6); /* here we shift by 6 since u[0] has weight
                                         1/2^64 and u[2] has weight 1/2^70, the
                                         truncation error on h+l/2^64 is <= 6/2^6 */
  sub_ddmmss (h, l, u[0], u[1], h, l);
  /* Since the mathematical error is < 0.412e-19*2^64, the total error on
     h + l/2^64 is less than 0.854; magic = ceil(0.854*2^64). */
  x = h - (l < magic && h != 0);

  /* now 2^64 + x is an approximation of 2^96/sqrt(np[1]),
     with 2^64 + x < 2^96/sqrt(np[1]) */

  umul_ppmm (r1, l, np[1], x);
  r1 += np[1];

  /* now r1 is an approximation of sqrt(2^64*np[1]), with r1 < sqrt(2^64*np[1]) */

  /* make sure r1 >= 2^63 */
  if (r1 < MPFR_LIMB_HIGHBIT)
    r1 = MPFR_LIMB_HIGHBIT;

  umul_ppmm (h, l, r1, r1);
  sub_ddmmss (h, l, np[1], np[0], h, l);
  l = (h << 63) | (l >> 1);
  h = h >> 1;

  /* now add (2^64+x) * (h*2^64+l) / 2^64 to r0 */

  umul_ppmm (r0, t, x, l); /* x * l */
  r0 += l;
  r1 += h + (r0 < l); /* now we have added 2^64 * (h*2^64+l) */
  while (h--)
    {
      r0 += x;
      r1 += (r0 < x); /* add x */
    }

  if (r1 & MPFR_LIMB_HIGHBIT)
    {
      rp[0] = r0;
      rp[1] = r1;
    }
  else /* overflow occurred in r1 */
    {
      rp[0] = ~MPFR_LIMB_ZERO;
      rp[1] = ~MPFR_LIMB_ZERO;
    }
}

/* Special code for prec(r), prec(u) < GMP_NUMB_BITS. We cannot have
   prec(u) = GMP_NUMB_BITS here, since when the exponent of u is odd,
   we need to shift u by one bit to the right without losing any bit.
   Assumes GMP_NUMB_BITS = 64. */
static int
mpfr_sqrt1 (mpfr_ptr r, mpfr_srcptr u, mpfr_rnd_t rnd_mode)
{
  mpfr_prec_t p = MPFR_GET_PREC(r);
  mpfr_prec_t exp_u = MPFR_EXP(u), exp_r, sh = GMP_NUMB_BITS - p;
  mp_limb_t u0, r0, rb, sb, mask = MPFR_LIMB_MASK(sh);
  mpfr_limb_ptr rp = MPFR_MANT(r);

  MPFR_STAT_STATIC_ASSERT (GMP_NUMB_BITS == 64);

  /* first make the exponent even */
  u0 = MPFR_MANT(u)[0];
  if (((unsigned int) exp_u & 1) != 0)
    {
      u0 >>= 1;
      exp_u ++;
    }
  MPFR_ASSERTD (((unsigned int) exp_u & 1) == 0);
  exp_r = exp_u / 2;

  /* then compute the integer square root of u0*2^GMP_NUMB_BITS */
  r0 = mpn_sqrtrem2_approx (u0);
  sb = 1; /* when we can round correctly with the approximation, the sticky bit
             is non-zero */

  /* Since the exact square root is in [r0 - 0.5406, r0 + 0.5406], we can round
     correctly except when the last sh-1 bits of r0 are 000...000. */
  if (MPFR_UNLIKELY((r0 & (mask >> 1)) == 0))
    {
      umul_ppmm (rb, sb, r0, r0);
      /* for the exact square root, we should have 0 <= (u0-rb)*2^64 - sb <= 2*r0 */
      if (rb > u0 || (rb == u0 && sb > 0)) /* r0 is too large */
        {
          r0 --;
          umul_ppmm (rb, sb, r0, r0);
        }
      /* if u0 <= rb + 1, then (u0-rb)*2^64 - sb <= 2^64 <= 2*r0
         if u0 >= rb + 3, then (u0-rb)*2^64 - sb > 2*2*64 > 2*r0 */
      else if (u0 > rb + 2 || (u0 == rb + 2 && -sb > 2 * r0))
        {
          r0 ++;
          umul_ppmm (rb, sb, r0, r0);
        }
      sub_ddmmss (rb, sb, u0, 0, rb, sb);
      /* now we should have rb*2^64 + sb <= 2*r0 */
      MPFR_ASSERTN(rb == 0 || (rb == 1 && sb <= 2 * r0));
      sb = rb | sb;
    }

  rb = r0 & (MPFR_LIMB_ONE << (sh - 1));
  sb |= (r0 & mask) ^ rb;
  rp[0] = r0 & ~mask;

  /* rounding */

  /* Note: if 1 and 2 are in [emin,emax], no overflow nor underflow
     is possible */
  if (MPFR_UNLIKELY (exp_r > __gmpfr_emax))
    return mpfr_overflow (r, rnd_mode, 1);

  /* See comments in mpfr_div_1 */
  if (MPFR_UNLIKELY (exp_r < __gmpfr_emin))
    {
      if (rnd_mode == MPFR_RNDN)
        {
          if ((exp_r == __gmpfr_emin - 1) && (rp[0] == ~mask) && rb)
            goto rounding; /* no underflow */
          if (exp_r < __gmpfr_emin - 1 || (rp[0] == MPFR_LIMB_HIGHBIT && sb == 0))
            rnd_mode = MPFR_RNDZ;
        }
      else if (!MPFR_IS_LIKE_RNDZ(rnd_mode, 0))
        {
          if ((exp_r == __gmpfr_emin - 1) && (rp[0] == ~mask) && (rb | sb))
            goto rounding; /* no underflow */
        }
      return mpfr_underflow (r, rnd_mode, 1);
    }

 rounding:
  MPFR_EXP (r) = exp_r;
  if (rb == 0 && sb == 0)
    {
      MPFR_ASSERTD(exp_r >= __gmpfr_emin);
      MPFR_ASSERTD(exp_r <= __gmpfr_emax);
      return 0; /* idem than MPFR_RET(0) but faster */
    }
  else if (rnd_mode == MPFR_RNDN)
    {
      if (rb == 0 || (rb && sb == 0 &&
                      (rp[0] & (MPFR_LIMB_ONE << sh)) == 0))
        goto truncate;
      else
        goto add_one_ulp;
    }
  else if (MPFR_IS_LIKE_RNDZ(rnd_mode, 0))
    {
    truncate:
      MPFR_ASSERTD(exp_r >= __gmpfr_emin);
      MPFR_ASSERTD(exp_r <= __gmpfr_emax);
      MPFR_RET(-1);
    }
  else /* round away from zero */
    {
    add_one_ulp:
      rp[0] += MPFR_LIMB_ONE << sh;
      if (rp[0] == 0)
        {
          rp[0] = MPFR_LIMB_HIGHBIT;
          if (MPFR_UNLIKELY(exp_r + 1 > __gmpfr_emax))
            return mpfr_overflow (r, rnd_mode, 1);
          MPFR_ASSERTD(exp_r + 1 <= __gmpfr_emax);
          MPFR_ASSERTD(exp_r + 1 >= __gmpfr_emin);
          MPFR_SET_EXP (r, exp_r + 1);
        }
      MPFR_RET(1);
    }
}

/* Special code for GMP_NUMB_BITS < prec(r) < 2*GMP_NUMB_BITS,
   and GMP_NUMB_BITS < prec(u) <= 2*GMP_NUMB_BITS.
   Assumes GMP_NUMB_BITS=64. */
static int
mpfr_sqrt2 (mpfr_ptr r, mpfr_srcptr u, mpfr_rnd_t rnd_mode)
{
  mpfr_prec_t p = MPFR_GET_PREC(r);
  mpfr_limb_ptr up = MPFR_MANT(u), rp = MPFR_MANT(r);
  mp_limb_t np[4], rb, sb, mask;
  mpfr_prec_t exp_u = MPFR_EXP(u), exp_r, sh = 2 * GMP_NUMB_BITS - p;

  MPFR_STAT_STATIC_ASSERT (GMP_NUMB_BITS == 64);

  if (((unsigned int) exp_u & 1) != 0)
    {
      np[3] = up[1] >> 1;
      np[2] = (up[1] << (GMP_NUMB_BITS - 1)) | (up[0] >> 1);
      np[1] = up[0] << (GMP_NUMB_BITS - 1);
      exp_u ++;
    }
  else
    {
      np[3] = up[1];
      np[2] = up[0];
      np[1] = 0;
    }
  exp_r = exp_u / 2;

  mask = MPFR_LIMB_MASK(sh);

  mpn_sqrtrem4_approx (rp, np + 2);
  /* the error is less than 43 ulps on rp[0] */
  if (((rp[0] + 42) & (mask >> 1)) > 84)
    sb = 1;
  else
    {
      mp_limb_t tp[4], cy;

      np[0] = 0;
      mpn_mul_n (tp, rp, rp, 2);
      while (mpn_cmp (tp, np, 4) > 0) /* approximation is too large */
        {
          mpn_sub_1 (rp, rp, 2, 1);
          /* subtract 2*{rp,2}+1 to {tp,4} */
          cy = mpn_submul_1 (tp, rp, 2, 2);
          mpn_sub_1 (tp + 2, tp + 2, 2, cy);
          mpn_sub_1 (tp, tp, 4, 1);
        }
      /* now {tp, 4} <= {np, 4} */
      mpn_sub_n (tp, np, tp, 4);
      /* now we want {tp, 4} <= 2 * {rp, 2}, which implies tp[2] <= 1 */
      MPFR_ASSERTD(tp[3] == 0);
      while (tp[2] > 1 ||
             (tp[2] == 1 && tp[1] > ((rp[1] << 1) | (rp[0] >> 63))) ||
             (tp[2] == 1 && tp[1] == ((rp[1] << 1) | (rp[0] >> 63)) &&
              tp[0] > (rp[0] << 1)))
        {
          /* subtract 2*{rp,2}+1 to {tp,3} (we know tp[3] = 0) */
          tp[2] -= mpn_submul_1 (tp, rp, 2, 2);
          mpn_sub_1 (tp, tp, 3, 1);
          mpn_add_1 (rp, rp, 2, 1);
        }
      sb = tp[2] | tp[0] | tp[1];
    }

  rb = rp[0] & (MPFR_LIMB_ONE << (sh - 1));
  sb |= (rp[0] & mask) ^ rb;
  rp[0] = rp[0] & ~mask;

  /* rounding */
  if (MPFR_UNLIKELY (exp_r > __gmpfr_emax))
    return mpfr_overflow (r, rnd_mode, 1);

  /* See comments in mpfr_div_1 */
  if (MPFR_UNLIKELY (exp_r < __gmpfr_emin))
    {
      if (rnd_mode == MPFR_RNDN)
        {
          if (exp_r == __gmpfr_emin - 1 && (rp[1] == MPFR_LIMB_MAX &&
                                            rp[0] == ~mask) && rb)
            goto rounding; /* no underflow */
          if (exp_r < __gmpfr_emin - 1 || (rp[1] == MPFR_LIMB_HIGHBIT &&
                                           rp[0] == MPFR_LIMB_ZERO && sb == 0))
            rnd_mode = MPFR_RNDZ;
        }
      else if (!MPFR_IS_LIKE_RNDZ(rnd_mode, 0))
        {
          if (exp_r == __gmpfr_emin - 1 && (rp[1] == MPFR_LIMB_MAX &&
                                            rp[0] == ~mask) && (rb | sb))
            goto rounding; /* no underflow */
        }
      return mpfr_underflow (r, rnd_mode, 1);
    }

 rounding:
  MPFR_EXP (r) = exp_r;
  if (rb == 0 && sb == 0)
    {
      MPFR_ASSERTD(exp_r >= __gmpfr_emin);
      MPFR_ASSERTD(exp_r <= __gmpfr_emax);
      return 0; /* idem than MPFR_RET(0) but faster */
    }
  else if (rnd_mode == MPFR_RNDN)
    {
      if (rb == 0 || (rb && sb == 0 &&
                      (rp[0] & (MPFR_LIMB_ONE << sh)) == 0))
        goto truncate;
      else
        goto add_one_ulp;
    }
  else if (MPFR_IS_LIKE_RNDZ(rnd_mode, 0))
    {
    truncate:
      MPFR_ASSERTD(exp_r >= __gmpfr_emin);
      MPFR_ASSERTD(exp_r <= __gmpfr_emax);
      MPFR_RET(-1);
    }
  else /* round away from zero */
    {
    add_one_ulp:
      rp[0] += MPFR_LIMB_ONE << sh;
      rp[1] += rp[0] == 0;
      if (rp[1] == 0)
        {
          rp[1] = MPFR_LIMB_HIGHBIT;
          if (MPFR_UNLIKELY(exp_r + 1 > __gmpfr_emax))
            return mpfr_overflow (r, rnd_mode, 1);
          MPFR_ASSERTD(exp_r + 1 <= __gmpfr_emax);
          MPFR_ASSERTD(exp_r + 1 >= __gmpfr_emin);
          MPFR_SET_EXP (r, exp_r + 1);
        }
      MPFR_RET(1);
    }
}
#endif /* !defined(MPFR_GENERIC_ABI) && GMP_NUMB_BITS == 64 */

int
mpfr_sqrt (mpfr_ptr r, mpfr_srcptr u, mpfr_rnd_t rnd_mode)
{
  mp_size_t rsize; /* number of limbs of r (plus 1 if exact limb multiple) */
  mp_size_t rrsize;
  mp_size_t usize; /* number of limbs of u */
  mp_size_t tsize; /* number of limbs of the sqrtrem remainder */
  mp_size_t k;
  mp_size_t l;
  mpfr_limb_ptr rp, rp0;
  mpfr_limb_ptr up;
  mpfr_limb_ptr sp;
  mp_limb_t sticky0; /* truncated part of input */
  mp_limb_t sticky1; /* truncated part of rp[0] */
  mp_limb_t sticky;
  int odd_exp;
  int sh; /* number of extra bits in rp[0] */
  int inexact; /* return ternary flag */
  mpfr_exp_t expr;
  MPFR_TMP_DECL(marker);

  MPFR_LOG_FUNC
    (("x[%Pu]=%.*Rg rnd=%d", mpfr_get_prec (u), mpfr_log_prec, u, rnd_mode),
     ("y[%Pu]=%.*Rg inexact=%d",
      mpfr_get_prec (r), mpfr_log_prec, r, inexact));

  if (MPFR_UNLIKELY(MPFR_IS_SINGULAR(u)))
    {
      if (MPFR_IS_NAN(u))
        {
          MPFR_SET_NAN(r);
          MPFR_RET_NAN;
        }
      else if (MPFR_IS_ZERO(u))
        {
          /* 0+ or 0- */
          MPFR_SET_SAME_SIGN(r, u);
          MPFR_SET_ZERO(r);
          MPFR_RET(0); /* zero is exact */
        }
      else
        {
          MPFR_ASSERTD(MPFR_IS_INF(u));
          /* sqrt(-Inf) = NAN */
          if (MPFR_IS_NEG(u))
            {
              MPFR_SET_NAN(r);
              MPFR_RET_NAN;
            }
          MPFR_SET_POS(r);
          MPFR_SET_INF(r);
          MPFR_RET(0);
        }
    }
  if (MPFR_UNLIKELY(MPFR_IS_NEG(u)))
    {
      MPFR_SET_NAN(r);
      MPFR_RET_NAN;
    }
  MPFR_SET_POS(r);

  /* See the note at the beginning of this file about __GNUC__. */
#if !defined(MPFR_GENERIC_ABI) && GMP_NUMB_BITS == 64
  if (MPFR_GET_PREC (r) < GMP_NUMB_BITS && MPFR_GET_PREC (u) < GMP_NUMB_BITS)
    return mpfr_sqrt1 (r, u, rnd_mode);

  if (GMP_NUMB_BITS < MPFR_GET_PREC (r) && MPFR_GET_PREC (r) < 2*GMP_NUMB_BITS
      && MPFR_LIMB_SIZE(u) == 2)
    return mpfr_sqrt2 (r, u, rnd_mode);
#endif

  MPFR_TMP_MARK (marker);
  MPFR_UNSIGNED_MINUS_MODULO (sh, MPFR_GET_PREC (r));
  if (sh == 0 && rnd_mode == MPFR_RNDN)
    sh = GMP_NUMB_BITS; /* ugly case */
  rsize = MPFR_LIMB_SIZE(r) + (sh == GMP_NUMB_BITS);
  /* rsize is the number of limbs of r + 1 if exact limb multiple and rounding
     to nearest, this is the number of wanted limbs for the square root */
  rrsize = rsize + rsize;
  usize = MPFR_LIMB_SIZE(u); /* number of limbs of u */
  rp0 = MPFR_MANT(r);
  rp = (sh < GMP_NUMB_BITS) ? rp0 : MPFR_TMP_LIMBS_ALLOC (rsize);
  up = MPFR_MANT(u);
  sticky0 = MPFR_LIMB_ZERO; /* truncated part of input */
  sticky1 = MPFR_LIMB_ZERO; /* truncated part of rp[0] */
  odd_exp = (unsigned int) MPFR_GET_EXP (u) & 1;
  inexact = -1; /* return ternary flag */

  sp = MPFR_TMP_LIMBS_ALLOC (rrsize);

  /* copy the most significant limbs of u to {sp, rrsize} */
  if (MPFR_LIKELY(usize <= rrsize)) /* in case r and u have the same precision,
                                       we have indeed rrsize = 2 * usize */
    {
      k = rrsize - usize;
      if (MPFR_LIKELY(k))
        MPN_ZERO (sp, k);
      if (odd_exp)
        {
          if (MPFR_LIKELY(k))
            sp[k - 1] = mpn_rshift (sp + k, up, usize, 1);
          else
            sticky0 = mpn_rshift (sp, up, usize, 1);
        }
      else
        MPN_COPY (sp + rrsize - usize, up, usize);
    }
  else /* usize > rrsize: truncate the input */
    {
      k = usize - rrsize;
      if (odd_exp)
        sticky0 = mpn_rshift (sp, up + k, rrsize, 1);
      else
        MPN_COPY (sp, up + k, rrsize);
      l = k;
      while (sticky0 == MPFR_LIMB_ZERO && l != 0)
        sticky0 = up[--l];
    }

  /* sticky0 is non-zero iff the truncated part of the input is non-zero */

  tsize = mpn_sqrtrem (rp, NULL, sp, rrsize);

  /* a return value of zero in mpn_sqrtrem indicates a perfect square */
  sticky = sticky0 || tsize != 0;

  /* truncate low bits of rp[0] */
  sticky1 = rp[0] & ((sh < GMP_NUMB_BITS) ? MPFR_LIMB_MASK(sh)
                     : MPFR_LIMB_MAX);
  rp[0] -= sticky1;

  sticky = sticky || sticky1;

  expr = (MPFR_GET_EXP(u) + odd_exp) / 2;  /* exact */

  if (rnd_mode == MPFR_RNDZ || rnd_mode == MPFR_RNDD || sticky == MPFR_LIMB_ZERO)
    {
      inexact = (sticky == MPFR_LIMB_ZERO) ? 0 : -1;
      goto truncate;
    }
  else if (rnd_mode == MPFR_RNDN)
    {
      /* if sh < GMP_NUMB_BITS, the round bit is bit (sh-1) of sticky1
                  and the sticky bit is formed by the low sh-1 bits from
                  sticky1, together with the sqrtrem remainder and sticky0. */
      if (sh < GMP_NUMB_BITS)
        {
          if (sticky1 & (MPFR_LIMB_ONE << (sh - 1)))
            { /* round bit is set */
              if (sticky1 == (MPFR_LIMB_ONE << (sh - 1)) && tsize == 0
                  && sticky0 == 0)
                goto even_rule;
              else
                goto add_one_ulp;
            }
          else /* round bit is zero */
            goto truncate; /* with the default inexact=-1 */
        }
      else /* sh = GMP_NUMB_BITS: the round bit is the most significant bit
              of rp[0], and the remaining GMP_NUMB_BITS-1 bits contribute to
              the sticky bit */
        {
          if (sticky1 & MPFR_LIMB_HIGHBIT)
            { /* round bit is set */
              if (sticky1 == MPFR_LIMB_HIGHBIT && tsize == 0 && sticky0 == 0)
                goto even_rule;
              else
                goto add_one_ulp;
            }
          else /* round bit is zero */
            goto truncate; /* with the default inexact=-1 */
        }
    }
  else /* rnd_mode=GMP_RDNU, necessarily sticky <> 0, thus add 1 ulp */
    goto add_one_ulp;

 even_rule: /* has to set inexact */
  if (sh < GMP_NUMB_BITS)
    inexact = (rp[0] & (MPFR_LIMB_ONE << sh)) ? 1 : -1;
  else
    inexact = (rp[1] & MPFR_LIMB_ONE) ? 1 : -1;
  if (inexact == -1)
    goto truncate;
  /* else go through add_one_ulp */

 add_one_ulp:
  inexact = 1; /* always here */
  if (sh == GMP_NUMB_BITS)
    {
      rp ++;
      rsize --;
      sh = 0;
    }
  /* now rsize = MPFR_LIMB_SIZE(r) */
  if (mpn_add_1 (rp0, rp, rsize, MPFR_LIMB_ONE << sh))
    {
      expr ++;
      rp0[rsize - 1] = MPFR_LIMB_HIGHBIT;
    }
  goto end;

 truncate: /* inexact = 0 or -1 */
  if (sh == GMP_NUMB_BITS)
    MPN_COPY (rp0, rp + 1, rsize - 1);

 end:
  /* Do not use MPFR_SET_EXP because the range has not been checked yet. */
  MPFR_ASSERTN (expr >= MPFR_EMIN_MIN && expr <= MPFR_EMAX_MAX);
  MPFR_EXP (r) = expr;
  MPFR_TMP_FREE(marker);

  return mpfr_check_range (r, inexact, rnd_mode);
}
