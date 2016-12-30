/* mpfr_sqr -- Floating square

Copyright 2004-2016 Free Software Foundation, Inc.
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

/* Special code for prec(a) < GMP_NUMB_BITS and prec(b) <= GMP_NUMB_BITS.
   Note: this function was copied from mpfr_mul_1 in file mul.c, thus any change
   here should be done also in mpfr_mul_1. */
static int
mpfr_sqr_1 (mpfr_ptr a, mpfr_srcptr b, mpfr_rnd_t rnd_mode, mpfr_prec_t p)
{
  mp_limb_t a0;
  mpfr_limb_ptr ap = MPFR_MANT(a);
  mpfr_limb_ptr bp = MPFR_MANT(b);
  mpfr_exp_t ax;
  mpfr_prec_t sh = GMP_NUMB_BITS - p;
  mp_limb_t rb, sb, mask = MPFR_LIMB_MASK(sh);

  /* When prec(b) <= GMP_NUMB_BITS / 2, we could replace umul_ppmm
     by a limb multiplication as follows, but we assume umul_ppmm is as fast
     as a limb multiplication on modern processors:
      a0 = (bp[0] >> (GMP_NUMB_BITS / 2)) * (bp[0] >> (GMP_NUMB_BITS / 2));
      sb = 0;
  */
  ax = MPFR_GET_EXP(b) << 1;
  umul_ppmm (a0, sb, bp[0], bp[0]);
  if (a0 < MPFR_LIMB_HIGHBIT)
    {
      ax --;
      a0 = (a0 << 1) | (sb >> (GMP_NUMB_BITS - 1));
      sb = sb << 1;
    }
  rb = a0 & (MPFR_LIMB_ONE << (sh - 1));
  sb |= (a0 & mask) ^ rb;
  ap[0] = a0 & ~mask;

  MPFR_SIGN(a) = MPFR_SIGN_POS;

  /* rounding */
  if (MPFR_UNLIKELY(ax > __gmpfr_emax))
    return mpfr_overflow (a, rnd_mode, MPFR_SIGN(a));

  /* Warning: underflow should be checked *after* rounding, thus when rounding
     away and when a > 0.111...111*2^(emin-1), or when rounding to nearest and
     a >= 0.111...111[1]*2^(emin-1), there is no underflow. */
  if (MPFR_UNLIKELY(ax < __gmpfr_emin))
    {
      if ((ax == __gmpfr_emin - 1) && (ap[0] == ~mask) &&
          ((rnd_mode == MPFR_RNDN && rb) ||
           (!MPFR_IS_LIKE_RNDZ(rnd_mode, MPFR_IS_NEG (a)) && (rb | sb))))
        goto rounding; /* no underflow */
      /* For RNDN, mpfr_underflow always rounds away, thus for |a| <= 2^(emin-2)
         we have to change to RNDZ. This corresponds to:
         (a) either ax < emin - 1
         (b) or ax = emin - 1 and ap[0] = 1000....000 and rb = sb = 0 */
      if (rnd_mode == MPFR_RNDN &&
          (ax < __gmpfr_emin - 1 || (ap[0] == MPFR_LIMB_HIGHBIT && (rb | sb) == 0)))
        rnd_mode = MPFR_RNDZ;
      return mpfr_underflow (a, rnd_mode, MPFR_SIGN(a));
    }

 rounding:
  MPFR_EXP (a) = ax; /* Don't use MPFR_SET_EXP since ax might be < __gmpfr_emin
                        in the cases "goto rounding" above. */
  if (rb == 0 && sb == 0)
    {
      MPFR_ASSERTD(ax >= __gmpfr_emin);
      return 0; /* idem than MPFR_RET(0) but faster */
    }
  else if (rnd_mode == MPFR_RNDN)
    {
      if (rb == 0 || (sb == 0 && (ap[0] & (MPFR_LIMB_ONE << sh)) == 0))
        goto truncate;
      else
        goto add_one_ulp;
    }
  else if (MPFR_IS_LIKE_RNDZ(rnd_mode, MPFR_IS_NEG(a)))
    {
    truncate:
      MPFR_ASSERTD(ax >= __gmpfr_emin);
      MPFR_RET(-MPFR_SIGN(a));
    }
  else /* round away from zero */
    {
    add_one_ulp:
      ap[0] += MPFR_LIMB_ONE << sh;
      if (ap[0] == 0)
        {
          ap[0] = MPFR_LIMB_HIGHBIT;
          if (MPFR_UNLIKELY(ax + 1 > __gmpfr_emax))
            return mpfr_overflow (a, rnd_mode, MPFR_SIGN(a));
          MPFR_ASSERTD(ax + 1 <= __gmpfr_emax);
          MPFR_ASSERTD(ax + 1 >= __gmpfr_emin);
          MPFR_SET_EXP (a, ax + 1);
        }
      MPFR_RET(MPFR_SIGN(a));
    }
}

/* Special code for GMP_NUMB_BITS < prec(a) < 2*GMP_NUMB_BITS and
   GMP_NUMB_BITS < prec(b) <= 2*GMP_NUMB_BITS.
   Note: this function was copied from mpfr_mul_2 in file mul.c, thus any change
   here should be done also in mpfr_mul_2. */
static int
mpfr_sqr_2 (mpfr_ptr a, mpfr_srcptr b, mpfr_rnd_t rnd_mode, mpfr_prec_t p)
{
  mp_limb_t h, l, u, v;
  mpfr_limb_ptr ap = MPFR_MANT(a);
  mpfr_exp_t ax = MPFR_GET_EXP(b) << 1;
  mpfr_prec_t sh = 2 * GMP_NUMB_BITS - p;
  mp_limb_t rb, sb, sb2, mask = MPFR_LIMB_MASK(sh);
  mp_limb_t *bp = MPFR_MANT(b);

  /* we store the 4-limb product in h=ap[1], l=ap[0], sb=ap[-1], sb2=ap[-2] */
  umul_ppmm (h, l, bp[1], bp[1]);
  umul_ppmm (sb, sb2, bp[0], bp[0]);
  umul_ppmm (u, v, bp[1], bp[0]);
  add_ssaaaa (l, sb, l, sb, u, v);
  /* warning: (l < u) is incorrect to detect a carry out of add_ssaaaa, since we
     might have u = 111...111, a carry coming from sb+v, thus l = u */
  h += (l < u) || (l == u && sb < v);
  add_ssaaaa (l, sb, l, sb, u, v);
  h += (l < u) || (l == u && sb < v);
  if (h < MPFR_LIMB_HIGHBIT)
    {
      ax --;
      h = (h << 1) | (l >> (GMP_NUMB_BITS - 1));
      l = (l << 1) | (sb >> (GMP_NUMB_BITS - 1));
      sb = sb << 1;
      /* no need to shift sb2 since we only want to know if it is zero or not */
    }
  ap[1] = h;
  rb = l & (MPFR_LIMB_ONE << (sh - 1));
  sb |= ((l & mask) ^ rb) | sb2;
  ap[0] = l & ~mask;

  MPFR_SIGN(a) = MPFR_SIGN_POS;

  /* rounding */
  if (MPFR_UNLIKELY(ax > __gmpfr_emax))
    return mpfr_overflow (a, rnd_mode, MPFR_SIGN(a));

  /* Warning: underflow should be checked *after* rounding, thus when rounding
     away and when a > 0.111...111*2^(emin-1), or when rounding to nearest and
     a >= 0.111...111[1]*2^(emin-1), there is no underflow. */
  if (MPFR_UNLIKELY(ax < __gmpfr_emin))
    {
      if ((ax == __gmpfr_emin - 1) &&
          (ap[1] == MPFR_LIMB_MAX) &&
          (ap[0] == ~mask) &&
          ((rnd_mode == MPFR_RNDN && rb) ||
           (!MPFR_IS_LIKE_RNDZ(rnd_mode, MPFR_IS_NEG (a)) && (rb | sb))))
        goto rounding; /* no underflow */
      /* for RNDN, mpfr_underflow always rounds away, thus for |a| <= 2^(emin-2)
         we have to change to RNDZ */
      if (rnd_mode == MPFR_RNDN &&
          (ax < __gmpfr_emin - 1 ||
           (ap[1] == MPFR_LIMB_HIGHBIT && ap[0] == 0 && (rb | sb) == 0)))
            rnd_mode = MPFR_RNDZ;
      return mpfr_underflow (a, rnd_mode, MPFR_SIGN(a));
    }

 rounding:
  MPFR_EXP (a) = ax; /* Don't use MPFR_SET_EXP since ax might be < __gmpfr_emin
                        in the cases "goto rounding" above. */
  if (rb == 0 && sb == 0)
    {
      MPFR_ASSERTD(ax >= __gmpfr_emin);
      return 0; /* idem than MPFR_RET(0) but faster */
    }
  else if (rnd_mode == MPFR_RNDN)
    {
      if (rb == 0 || (sb == 0 && (ap[0] & (MPFR_LIMB_ONE << sh)) == 0))
        goto truncate;
      else
        goto add_one_ulp;
    }
  else if (MPFR_IS_LIKE_RNDZ(rnd_mode, MPFR_IS_NEG(a)))
    {
    truncate:
      MPFR_ASSERTD(ax >= __gmpfr_emin);
      MPFR_RET(-MPFR_SIGN(a));
    }
  else /* round away from zero */
    {
    add_one_ulp:
      ap[0] += MPFR_LIMB_ONE << sh;
      ap[1] += (ap[0] == 0);
      if (ap[1] == 0)
        {
          ap[1] = MPFR_LIMB_HIGHBIT;
          if (MPFR_UNLIKELY(ax + 1 > __gmpfr_emax))
            return mpfr_overflow (a, rnd_mode, MPFR_SIGN(a));
          MPFR_ASSERTD(ax + 1 <= __gmpfr_emax);
          MPFR_ASSERTD(ax + 1 >= __gmpfr_emin);
          MPFR_SET_EXP (a, ax + 1);
        }
      MPFR_RET(MPFR_SIGN(a));
    }
}

int
mpfr_sqr (mpfr_ptr a, mpfr_srcptr b, mpfr_rnd_t rnd_mode)
{
  int cc, inexact;
  mpfr_exp_t ax;
  mp_limb_t *tmp;
  mp_limb_t b1;
  mpfr_prec_t bq;
  mp_size_t bn, tn;
  MPFR_TMP_DECL(marker);

  MPFR_LOG_FUNC
    (("x[%Pu]=%.*Rg rnd=%d", mpfr_get_prec (b), mpfr_log_prec, b, rnd_mode),
     ("y[%Pu]=%.*Rg inexact=%d",
      mpfr_get_prec (a), mpfr_log_prec, a, inexact));

  /* deal with special cases */
  if (MPFR_UNLIKELY(MPFR_IS_SINGULAR(b)))
    {
      if (MPFR_IS_NAN(b))
        {
          MPFR_SET_NAN(a);
          MPFR_RET_NAN;
        }
      MPFR_SET_POS (a);
      if (MPFR_IS_INF(b))
        MPFR_SET_INF(a);
      else
        ( MPFR_ASSERTD(MPFR_IS_ZERO(b)), MPFR_SET_ZERO(a) );
      MPFR_RET(0);
    }
  bq = MPFR_PREC(b);

  if (MPFR_GET_PREC(a) < GMP_NUMB_BITS && bq <= GMP_NUMB_BITS)
    return mpfr_sqr_1 (a, b, rnd_mode, MPFR_GET_PREC(a));

  if (GMP_NUMB_BITS < MPFR_GET_PREC(a) && MPFR_GET_PREC(a) < 2 * GMP_NUMB_BITS
      && GMP_NUMB_BITS < bq && bq <= 2 * GMP_NUMB_BITS)
    return mpfr_sqr_2 (a, b, rnd_mode, MPFR_GET_PREC(a));

  ax = 2 * MPFR_GET_EXP (b);
  MPFR_ASSERTN (2 * (mpfr_uprec_t) bq <= MPFR_PREC_MAX);

  bn = MPFR_LIMB_SIZE (b); /* number of limbs of b */
  tn = MPFR_PREC2LIMBS (2 * bq); /* number of limbs of square,
                                    2*bn or 2*bn-1 */

  if (MPFR_UNLIKELY(bn > MPFR_SQR_THRESHOLD))
    return mpfr_mul (a, b, b, rnd_mode);

  MPFR_TMP_MARK(marker);
  tmp = MPFR_TMP_LIMBS_ALLOC (2 * bn);

  /* Multiplies the mantissa in temporary allocated space */
  mpn_sqr_n (tmp, MPFR_MANT(b), bn);
  b1 = tmp[2 * bn - 1];

  /* now tmp[0]..tmp[2*bn-1] contains the product of both mantissa,
     with tmp[2*bn-1]>=2^(GMP_NUMB_BITS-2) */
  b1 >>= GMP_NUMB_BITS - 1; /* msb from the product */

  /* if the mantissas of b and c are uniformly distributed in ]1/2, 1],
     then their product is in ]1/4, 1/2] with probability 2*ln(2)-1 ~ 0.386
     and in [1/2, 1] with probability 2-2*ln(2) ~ 0.614 */
  tmp += 2 * bn - tn; /* +0 or +1 */
  if (MPFR_UNLIKELY(b1 == 0))
    mpn_lshift (tmp, tmp, tn, 1); /* tn <= k, so no stack corruption */

  cc = mpfr_round_raw (MPFR_MANT (a), tmp, 2 * bq, 0,
                       MPFR_PREC (a), rnd_mode, &inexact);
  /* cc = 1 ==> result is a power of two */
  if (MPFR_UNLIKELY(cc))
    MPFR_MANT(a)[MPFR_LIMB_SIZE(a)-1] = MPFR_LIMB_HIGHBIT;

  MPFR_TMP_FREE(marker);
  {
    mpfr_exp_t ax2 = ax + (mpfr_exp_t) (b1 - 1 + cc);
    if (MPFR_UNLIKELY( ax2 > __gmpfr_emax))
      return mpfr_overflow (a, rnd_mode, MPFR_SIGN_POS);
    if (MPFR_UNLIKELY( ax2 < __gmpfr_emin))
      {
        /* In the rounding to the nearest mode, if the exponent of the exact
           result (i.e. before rounding, i.e. without taking cc into account)
           is < __gmpfr_emin - 1 or the exact result is a power of 2 (i.e. if
           both arguments are powers of 2), then round to zero. */
        if (rnd_mode == MPFR_RNDN &&
            (ax + (mpfr_exp_t) b1 < __gmpfr_emin || mpfr_powerof2_raw (b)))
          rnd_mode = MPFR_RNDZ;
        return mpfr_underflow (a, rnd_mode, MPFR_SIGN_POS);
      }
    MPFR_SET_EXP (a, ax2);
    MPFR_SET_POS (a);
  }
  MPFR_RET (inexact);
}
