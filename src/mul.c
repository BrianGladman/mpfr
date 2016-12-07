/* mpfr_mul -- multiply two floating-point numbers

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


/********* BEGINNING CHECK *************/

/* Check if we have to check the result of mpfr_mul.
   TODO: Find a better (and faster?) check than using old implementation */
#if MPFR_WANT_ASSERT >= 2

int mpfr_mul2 (mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mpfr_rnd_t rnd_mode);
static int
mpfr_mul3 (mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mpfr_rnd_t rnd_mode)
{
  /* Old implementation */
  int sign_product, cc, inexact;
  mpfr_exp_t ax;
  mp_limb_t *tmp;
  mp_limb_t b1;
  mpfr_prec_t bq, cq;
  mp_size_t bn, cn, tn, k;
  MPFR_TMP_DECL(marker);

  /* deal with special cases */
  if (MPFR_ARE_SINGULAR(b,c))
    {
      if (MPFR_IS_NAN(b) || MPFR_IS_NAN(c))
        {
          MPFR_SET_NAN(a);
          MPFR_RET_NAN;
        }
      sign_product = MPFR_MULT_SIGN( MPFR_SIGN(b) , MPFR_SIGN(c) );
      if (MPFR_IS_INF(b))
        {
          if (MPFR_IS_INF(c) || MPFR_NOTZERO(c))
            {
              MPFR_SET_SIGN(a, sign_product);
              MPFR_SET_INF(a);
              MPFR_RET(0); /* exact */
            }
          else
            {
              MPFR_SET_NAN(a);
              MPFR_RET_NAN;
            }
        }
      else if (MPFR_IS_INF(c))
        {
          if (MPFR_NOTZERO(b))
            {
              MPFR_SET_SIGN(a, sign_product);
              MPFR_SET_INF(a);
              MPFR_RET(0); /* exact */
            }
          else
            {
              MPFR_SET_NAN(a);
              MPFR_RET_NAN;
            }
        }
      else
        {
          MPFR_ASSERTD(MPFR_IS_ZERO(b) || MPFR_IS_ZERO(c));
          MPFR_SET_SIGN(a, sign_product);
          MPFR_SET_ZERO(a);
          MPFR_RET(0); /* 0 * 0 is exact */
        }
    }
  sign_product = MPFR_MULT_SIGN( MPFR_SIGN(b) , MPFR_SIGN(c) );

  ax = MPFR_GET_EXP (b) + MPFR_GET_EXP (c);

  bq = MPFR_PREC (b);
  cq = MPFR_PREC (c);

  MPFR_ASSERTN ((mpfr_uprec_t) bq + cq <= MPFR_PREC_MAX);

  bn = MPFR_PREC2LIMBS (bq); /* number of limbs of b */
  cn = MPFR_PREC2LIMBS (cq); /* number of limbs of c */
  k = bn + cn; /* effective nb of limbs used by b*c (= tn or tn+1) below */
  tn = MPFR_PREC2LIMBS (bq + cq);
  /* <= k, thus no int overflow */
  MPFR_ASSERTD(tn <= k);

  /* Check for no size_t overflow*/
  MPFR_ASSERTD((size_t) k <= ((size_t) -1) / MPFR_BYTES_PER_MP_LIMB);
  MPFR_TMP_MARK(marker);
  tmp = MPFR_TMP_LIMBS_ALLOC (k);

  /* multiplies two mantissa in temporary allocated space */
  b1 = (MPFR_LIKELY(bn >= cn)) ?
    mpn_mul (tmp, MPFR_MANT(b), bn, MPFR_MANT(c), cn)
    : mpn_mul (tmp, MPFR_MANT(c), cn, MPFR_MANT(b), bn);

  /* now tmp[0]..tmp[k-1] contains the product of both mantissa,
     with tmp[k-1]>=2^(GMP_NUMB_BITS-2) */
  b1 >>= GMP_NUMB_BITS - 1; /* msb from the product */

  /* if the mantissas of b and c are uniformly distributed in ]1/2, 1],
     then their product is in ]1/4, 1/2] with probability 2*ln(2)-1 ~ 0.386
     and in [1/2, 1] with probability 2-2*ln(2) ~ 0.614 */
  tmp += k - tn;
  if (MPFR_UNLIKELY(b1 == 0))
    mpn_lshift (tmp, tmp, tn, 1); /* tn <= k, so no stack corruption */
  cc = mpfr_round_raw (MPFR_MANT (a), tmp, bq + cq,
                       MPFR_IS_NEG_SIGN(sign_product),
                       MPFR_PREC (a), rnd_mode, &inexact);

  /* cc = 1 ==> result is a power of two */
  if (MPFR_UNLIKELY(cc))
    MPFR_MANT(a)[MPFR_LIMB_SIZE(a)-1] = MPFR_LIMB_HIGHBIT;

  MPFR_TMP_FREE(marker);

  {
    mpfr_exp_t ax2 = ax + (mpfr_exp_t) (b1 - 1 + cc);
    if (MPFR_UNLIKELY( ax2 > __gmpfr_emax))
      return mpfr_overflow (a, rnd_mode, sign_product);
    if (MPFR_UNLIKELY( ax2 < __gmpfr_emin))
      {
        /* In the rounding to the nearest mode, if the exponent of the exact
           result (i.e. before rounding, i.e. without taking cc into account)
           is < __gmpfr_emin - 1 or the exact result is a power of 2 (i.e. if
           both arguments are powers of 2) in absolute value, then round to
           zero. */
        if (rnd_mode == MPFR_RNDN &&
            (ax + (mpfr_exp_t) b1 < __gmpfr_emin ||
             (mpfr_powerof2_raw (b) && mpfr_powerof2_raw (c))))
          rnd_mode = MPFR_RNDZ;
        return mpfr_underflow (a, rnd_mode, sign_product);
      }
    MPFR_SET_EXP (a, ax2);
    MPFR_SET_SIGN(a, sign_product);
  }
  MPFR_RET (inexact);
}

int
mpfr_mul (mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mpfr_rnd_t rnd_mode)
{
  mpfr_t ta, tb, tc;
  int inexact1, inexact2;

  if (rnd_mode == MPFR_RNDF)
    return mpfr_mul2 (a, b, c, rnd_mode);

  mpfr_init2 (ta, MPFR_PREC (a));
  mpfr_init2 (tb, MPFR_PREC (b));
  mpfr_init2 (tc, MPFR_PREC (c));
  MPFR_ASSERTN (mpfr_set (tb, b, MPFR_RNDN) == 0);
  MPFR_ASSERTN (mpfr_set (tc, c, MPFR_RNDN) == 0);

  inexact2 = mpfr_mul3 (ta, tb, tc, rnd_mode);
  inexact1 = mpfr_mul2 (a, b, c, rnd_mode);
  if (MPFR_IS_NAN (ta) && MPFR_IS_NAN (a))
    {
      /* Getting both NaN is OK. */
    }
  else if (! mpfr_equal_p (ta, a) || ! SAME_SIGN (inexact1, inexact2))
    {
      fprintf (stderr, "mpfr_mul return different values for %s\n"
               "Prec_a = %lu, Prec_b = %lu, Prec_c = %lu\nb = ",
               mpfr_print_rnd_mode (rnd_mode),
               MPFR_PREC (a), MPFR_PREC (b), MPFR_PREC (c));
      mpfr_fprint_binary (stderr, b);
      fprintf (stderr, "\nc = ");
      mpfr_fprint_binary (stderr, c);
      fprintf (stderr, "\nOldMul: ");
      mpfr_fprint_binary (stderr, ta);
      fprintf (stderr, "\nNewMul: ");
      mpfr_fprint_binary (stderr, a);
      fprintf (stderr, "\nNewInexact = %d | OldInexact = %d\n",
               inexact1, inexact2);
      MPFR_ASSERTN(0);
    }

  mpfr_clears (ta, tb, tc, (mpfr_ptr) 0);
  return inexact1;
}

# define mpfr_mul mpfr_mul2

#endif  /* MPFR_WANT_ASSERT >= 2 */

/****** END OF CHECK *******/

/* Multiply 2 mpfr_t */

/* special code for prec(a) < GMP_NUMB_BITS and
   prec(b), prec(c) <= GMP_NUMB_BITS */
static int
mpfr_mul_1 (mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mpfr_rnd_t rnd_mode,
              mpfr_prec_t p)
{
  mp_limb_t h;
  mpfr_limb_ptr ap = MPFR_MANT(a);
  mpfr_exp_t ax = MPFR_GET_EXP(b) + MPFR_GET_EXP(c);
  mpfr_prec_t sh = GMP_NUMB_BITS - p;
  mp_limb_t rb, sb, mask = MPFR_LIMB_MASK(sh);

  /* When prec(b), prec(c) <= GMP_NUMB_BITS / 2, we could replace umul_ppmm
     by a limb multiplication as follows, but we assume umul_ppmm is as fast
     as a limb multiplication on modern processors:
      h = (MPFR_MANT(b)[0] >> (GMP_NUMB_BITS / 2))
        * (MPFR_MANT(c)[0] >> (GMP_NUMB_BITS / 2));
      sb = 0;
  */
  umul_ppmm (h, sb, MPFR_MANT(b)[0], MPFR_MANT(c)[0]);
  if (h < MPFR_LIMB_HIGHBIT)
    {
      ax --;
      h = (h << 1) | (sb >> (GMP_NUMB_BITS - 1));
      sb = sb << 1;
    }
  rb = h & (MPFR_LIMB_ONE << (sh - 1));
  sb |= (h & mask) ^ rb;
  ap[0] = h & ~mask;

  MPFR_SIGN(a) = MPFR_MULT_SIGN (MPFR_SIGN (b), MPFR_SIGN (c));

  /* rounding */
  if (MPFR_UNLIKELY(ax > __gmpfr_emax))
    return mpfr_overflow (a, rnd_mode, MPFR_SIGN(a));

  /* Warning: underflow should be checked *after* rounding, thus when rounding
     away and when a > 0.111...111*2^(emin-1), or when rounding to nearest and
     a >= 0.111...111[1]*2^(emin-1), there is no underflow. */
  if (MPFR_UNLIKELY(ax < __gmpfr_emin))
    {
      /* for RNDN, mpfr_underflow always rounds away, thus for |a| <= 2^(emin-2)
         we have to change to RNDZ */
      if (rnd_mode == MPFR_RNDN)
        {
          if ((ax == __gmpfr_emin - 1) && (ap[0] == ~mask) && rb)
            goto rounding; /* no underflow */
          if (ax < __gmpfr_emin - 1 || (ap[0] == MPFR_LIMB_HIGHBIT && sb == 0))
            rnd_mode = MPFR_RNDZ;
        }
      else if (!MPFR_IS_LIKE_RNDZ(rnd_mode, MPFR_IS_NEG (a)))
        {
          if ((ax == __gmpfr_emin - 1) && (ap[0] == ~mask) && (rb | sb))
            goto rounding; /* no underflow */
        }
      return mpfr_underflow (a, rnd_mode, MPFR_SIGN(a));
    }

 rounding:
  MPFR_EXP (a) = ax; /* Don't use MPFR_SET_EXP since ax might be < __gmpfr_emin
                        in the cases "goto rounding" above. */
  if ((rb == 0 && sb == 0) || (rnd_mode == MPFR_RNDF))
    {
      MPFR_ASSERTD(ax >= __gmpfr_emin);
      return 0; /* idem than MPFR_RET(0) but faster */
    }
  else if (rnd_mode == MPFR_RNDN)
    {
      if (rb == 0 || (rb && sb == 0 &&
                      (ap[0] & (MPFR_LIMB_ONE << sh)) == 0))
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

/* special code for GMP_NUMB_BITS < prec(a) < 2*GMP_NUMB_BITS and
   GMP_NUMB_BITS < prec(b), prec(c) <= 2*GMP_NUMB_BITS */
static int
mpfr_mul_2 (mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mpfr_rnd_t rnd_mode,
              mpfr_prec_t p)
{
  mp_limb_t h, l, u, v;
  mpfr_limb_ptr ap = MPFR_MANT(a);
  mpfr_exp_t ax = MPFR_GET_EXP(b) + MPFR_GET_EXP(c);
  mpfr_prec_t sh = 2 * GMP_NUMB_BITS - p;
  mp_limb_t rb, sb, sb2, mask = MPFR_LIMB_MASK(sh);
  mp_limb_t *bp = MPFR_MANT(b), *cp = MPFR_MANT(c);

  /* we store the 4-limb product in h=ap[1], l=ap[0], sb=ap[-1], sb2=ap[-2] */
  umul_ppmm (h, l, bp[1], cp[1]);
  umul_ppmm (sb, sb2, bp[0], cp[0]);
  umul_ppmm (u, v, bp[1], cp[0]);
  add_ssaaaa (l, sb, l, sb, u, v);
  /* warning: (l < u) is incorrect to detect a carry out of add_ssaaaa, since we
     might have u = 111...111, a carry coming from l+v, thus l = u */
  h += (l < u) || (l == u && sb < v);
  umul_ppmm (u, v, bp[0], cp[1]);
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

  MPFR_SIGN(a) = MPFR_MULT_SIGN (MPFR_SIGN (b), MPFR_SIGN (c));

  /* rounding */
  if (MPFR_UNLIKELY(ax > __gmpfr_emax))
    return mpfr_overflow (a, rnd_mode, MPFR_SIGN(a));

  /* Warning: underflow should be checked *after* rounding, thus when rounding
     away and when a > 0.111...111*2^(emin-1), or when rounding to nearest and
     a >= 0.111...111[1]*2^(emin-1), there is no underflow. */
  if (MPFR_UNLIKELY(ax < __gmpfr_emin))
    {
      /* for RNDN, mpfr_underflow always rounds away, thus for |a| <= 2^(emin-2)
         we have to change to RNDZ */
      if (rnd_mode == MPFR_RNDN)
        {
          if ((ax == __gmpfr_emin - 1) && (~ap[1] == 0) && (ap[0] == ~mask) && rb)
            goto rounding; /* no underflow */
          if (ax < __gmpfr_emin - 1 ||
              (ap[1] == MPFR_LIMB_HIGHBIT && ap[0] == 0 && sb == 0))
            rnd_mode = MPFR_RNDZ;
        }
      else if (!MPFR_IS_LIKE_RNDZ(rnd_mode, MPFR_IS_NEG (a)))
        {
          if ((ax == __gmpfr_emin - 1) && (~ap[1] == 0) && (ap[0] == ~mask)
              && (rb | sb))
            goto rounding; /* no underflow */
        }
      return mpfr_underflow (a, rnd_mode, MPFR_SIGN(a));
    }

 rounding:
  MPFR_EXP (a) = ax; /* Don't use MPFR_SET_EXP since ax might be < __gmpfr_emin
                        in the cases "goto rounding" above. */
  if ((rb == 0 && sb == 0) || (rnd_mode == MPFR_RNDF))
    {
      MPFR_ASSERTD(ax >= __gmpfr_emin);
      return 0; /* idem than MPFR_RET(0) but faster */
    }
  else if (rnd_mode == MPFR_RNDN)
    {
      if (rb == 0 || (rb && sb == 0 &&
                      (ap[0] & (MPFR_LIMB_ONE << sh)) == 0))
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

/* Note: mpfr_sqr will call mpfr_mul if bn > MPFR_SQR_THRESHOLD,
   in order to use Mulders' mulhigh, which is handled only here
   to avoid partial code duplication. There is some overhead due
   to the additional tests, but slowdown should not be noticeable
   as this code is not executed in very small precisions. */

MPFR_HOT_FUNCTION_ATTR int
mpfr_mul (mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mpfr_rnd_t rnd_mode)
{
  int sign, inexact;
  mpfr_exp_t ax, ax2;
  mp_limb_t *tmp;
  mp_limb_t b1;
  mpfr_prec_t bq, cq;
  mp_size_t bn, cn, tn, k, threshold;
  MPFR_TMP_DECL (marker);

  MPFR_LOG_FUNC
    (("b[%Pu]=%.*Rg c[%Pu]=%.*Rg rnd=%d",
      mpfr_get_prec (b), mpfr_log_prec, b,
      mpfr_get_prec (c), mpfr_log_prec, c, rnd_mode),
     ("a[%Pu]=%.*Rg inexact=%d",
      mpfr_get_prec (a), mpfr_log_prec, a, inexact));

  /* deal with special cases */
  if (MPFR_ARE_SINGULAR (b, c))
    {
      if (MPFR_IS_NAN (b) || MPFR_IS_NAN (c))
        {
          MPFR_SET_NAN (a);
          MPFR_RET_NAN;
        }
      sign = MPFR_MULT_SIGN (MPFR_SIGN (b), MPFR_SIGN (c));
      if (MPFR_IS_INF (b))
        {
          if (!MPFR_IS_ZERO (c))
            {
              MPFR_SET_SIGN (a, sign);
              MPFR_SET_INF (a);
              MPFR_RET (0);
            }
          else
            {
              MPFR_SET_NAN (a);
              MPFR_RET_NAN;
            }
        }
      else if (MPFR_IS_INF (c))
        {
          if (!MPFR_IS_ZERO (b))
            {
              MPFR_SET_SIGN (a, sign);
              MPFR_SET_INF (a);
              MPFR_RET(0);
            }
          else
            {
              MPFR_SET_NAN (a);
              MPFR_RET_NAN;
            }
        }
      else
        {
          MPFR_ASSERTD (MPFR_IS_ZERO(b) || MPFR_IS_ZERO(c));
          MPFR_SET_SIGN (a, sign);
          MPFR_SET_ZERO (a);
          MPFR_RET (0);
        }
    }

  bq = MPFR_GET_PREC (b);
  cq = MPFR_GET_PREC (c);
  if (MPFR_GET_PREC(a) < GMP_NUMB_BITS &&
      bq <= GMP_NUMB_BITS && cq <= GMP_NUMB_BITS)
    return mpfr_mul_1 (a, b, c, rnd_mode, MPFR_GET_PREC(a));

  if (GMP_NUMB_BITS < MPFR_GET_PREC(a) && MPFR_GET_PREC(a) < 2 * GMP_NUMB_BITS
      && GMP_NUMB_BITS < bq && bq <= 2 * GMP_NUMB_BITS
      && GMP_NUMB_BITS < cq && cq <= 2 * GMP_NUMB_BITS)
    return mpfr_mul_2 (a, b, c, rnd_mode, MPFR_GET_PREC(a));

  sign = MPFR_MULT_SIGN (MPFR_SIGN (b), MPFR_SIGN (c));

  ax = MPFR_GET_EXP (b) + MPFR_GET_EXP (c);
  /* Note: the exponent of the exact result will be e = bx + cx + ec with
     ec in {-1,0,1} and the following assumes that e is representable. */

  /* FIXME: Useful since we do an exponent check after?
   * It is useful iff the precision is big, there is an overflow
   * and we are doing further mults...*/
#ifdef HUGE
  if (MPFR_UNLIKELY (ax > __gmpfr_emax + 1))
    return mpfr_overflow (a, rnd_mode, sign);
  if (MPFR_UNLIKELY (ax < __gmpfr_emin - 2))
    return mpfr_underflow (a, rnd_mode == MPFR_RNDN ? MPFR_RNDZ : rnd_mode,
                           sign);
#endif

  MPFR_ASSERTN ((mpfr_uprec_t) bq + cq <= MPFR_PREC_MAX);

  bn = MPFR_PREC2LIMBS (bq); /* number of limbs of b */
  cn = MPFR_PREC2LIMBS (cq); /* number of limbs of c */
  k = bn + cn; /* effective nb of limbs used by b*c (= tn or tn+1) below */
  tn = MPFR_PREC2LIMBS (bq + cq);
  MPFR_ASSERTD (tn <= k); /* tn <= k, thus no int overflow */

  /* Check for no size_t overflow. */
  MPFR_ASSERTD ((size_t) k <= ((size_t) -1) / MPFR_BYTES_PER_MP_LIMB);
  MPFR_TMP_MARK (marker);
  tmp = MPFR_TMP_LIMBS_ALLOC (k);

  /* multiplies two mantissa in temporary allocated space */
  if (MPFR_UNLIKELY (bn < cn))
    {
      mpfr_srcptr z = b;
      mp_size_t zn  = bn;
      b = c;
      bn = cn;
      c = z;
      cn = zn;
    }
  MPFR_ASSERTD (bn >= cn);
  if (bn <= 2)
    {
      /* The 3 cases perform the same first operation. */
      umul_ppmm (tmp[1], tmp[0], MPFR_MANT (b)[0], MPFR_MANT (c)[0]);
      if (bn == 1)
        {
          /* 1 limb * 1 limb */
          b1 = tmp[1];
        }
      else if (MPFR_UNLIKELY (cn == 1))
        {
          /* 2 limbs * 1 limb */
          mp_limb_t t;
          umul_ppmm (tmp[2], t, MPFR_MANT (b)[1], MPFR_MANT (c)[0]);
          add_ssaaaa (tmp[2], tmp[1], tmp[2], tmp[1], 0, t);
          b1 = tmp[2];
        }
      else
        {
          /* 2 limbs * 2 limbs */
          mp_limb_t t1, t2, t3;
          /* First 2 limbs * 1 limb */
          umul_ppmm (tmp[2], t1, MPFR_MANT (b)[1], MPFR_MANT (c)[0]);
          add_ssaaaa (tmp[2], tmp[1], tmp[2], tmp[1], 0, t1);
          /* Second, the other 2 limbs * 1 limb product */
          umul_ppmm (t1, t2, MPFR_MANT (b)[0], MPFR_MANT (c)[1]);
          umul_ppmm (tmp[3], t3, MPFR_MANT (b)[1], MPFR_MANT (c)[1]);
          add_ssaaaa (tmp[3], t1, tmp[3], t1, 0, t3);
          /* Sum those two partial products */
          add_ssaaaa (tmp[2], tmp[1], tmp[2], tmp[1], t1, t2);
          tmp[3] += (tmp[2] < t1);
          b1 = tmp[3];
        }
      b1 >>= (GMP_NUMB_BITS - 1);
      tmp += k - tn;
      if (MPFR_UNLIKELY (b1 == 0))
        mpn_lshift (tmp, tmp, tn, 1); /* tn <= k, so no stack corruption */
    }
  else
    /* Mulders' mulhigh. This code can also be used via mpfr_sqr,
       hence the tests b != c. */
    if (MPFR_UNLIKELY (bn > (threshold = b != c ?
                             MPFR_MUL_THRESHOLD : MPFR_SQR_THRESHOLD)))
      {
        mp_limb_t *bp, *cp;
        mp_size_t n;
        mpfr_prec_t p;

        /* First check if we can reduce the precision of b or c:
           exact values are a nightmare for the short product trick */
        bp = MPFR_MANT (b);
        cp = MPFR_MANT (c);
        MPFR_STAT_STATIC_ASSERT (MPFR_MUL_THRESHOLD >= 1 &&
                                 MPFR_SQR_THRESHOLD >= 1);
        if (MPFR_UNLIKELY ((bp[0] == 0 && bp[1] == 0) ||
                           (cp[0] == 0 && cp[1] == 0)))
          {
            mpfr_t b_tmp, c_tmp;

            MPFR_TMP_FREE (marker);
            /* Check for b */
            while (*bp == 0)
              {
                bp++;
                bn--;
                MPFR_ASSERTD (bn > 0);
              } /* This must end since the most significant limb is != 0 */

            /* Check for c too: if b == c, this will do nothing */
            while (*cp == 0)
              {
                cp++;
                cn--;
                MPFR_ASSERTD (cn > 0);
              } /* This must end since the most significant limb is != 0 */

            /* It is not the fastest way, but it is safer. */
            MPFR_SET_SAME_SIGN (b_tmp, b);
            MPFR_SET_EXP (b_tmp, MPFR_GET_EXP (b));
            MPFR_PREC (b_tmp) = bn * GMP_NUMB_BITS;
            MPFR_MANT (b_tmp) = bp;

            if (b != c)
              {
                MPFR_SET_SAME_SIGN (c_tmp, c);
                MPFR_SET_EXP (c_tmp, MPFR_GET_EXP (c));
                MPFR_PREC (c_tmp) = cn * GMP_NUMB_BITS;
                MPFR_MANT (c_tmp) = cp;

                /* Call again mpfr_mul with the fixed arguments */
                return mpfr_mul (a, b_tmp, c_tmp, rnd_mode);
              }
            else
              /* Call mpfr_mul instead of mpfr_sqr as the precision
                 is probably still high enough. */
              return mpfr_mul (a, b_tmp, b_tmp, rnd_mode);
          }

        /* Compute estimated precision of mulhigh.
           We could use `+ (n < cn) + (n < bn)' instead of `+ 2',
           but does it worth it? */
        n = MPFR_LIMB_SIZE (a) + 1;
        n = MIN (n, cn);
        MPFR_ASSERTD (n >= 1 && 2*n <= k && n <= cn && n <= bn);
        p = n * GMP_NUMB_BITS - MPFR_INT_CEIL_LOG2 (n + 2);
        bp += bn - n;
        cp += cn - n;

        /* Check if MulHigh can produce a roundable result.
           We may lose 1 bit due to RNDN, 1 due to final shift. */
        if (MPFR_UNLIKELY (MPFR_GET_PREC (a) > p - 5))
          {
            if (MPFR_UNLIKELY (MPFR_GET_PREC (a) > p - 5 + GMP_NUMB_BITS
                               || bn <= threshold + 1))
              {
                /* MulHigh can't produce a roundable result. */
                MPFR_LOG_MSG (("mpfr_mulhigh can't be used (%lu VS %lu)\n",
                               MPFR_GET_PREC (a), p));
                goto full_multiply;
              }
            /* Add one extra limb to mantissa of b and c. */
            if (bn > n)
              bp --;
            else
              {
                bp = MPFR_TMP_LIMBS_ALLOC (n + 1);
                bp[0] = 0;
                MPN_COPY (bp + 1, MPFR_MANT (b) + bn - n, n);
              }
            if (b != c)
              {
                if (cn > n)
                  cp --; /* FIXME: Could this happen? */
                else
                  {
                    cp = MPFR_TMP_LIMBS_ALLOC (n + 1);
                    cp[0] = 0;
                    MPN_COPY (cp + 1, MPFR_MANT (c) + cn - n, n);
                  }
              }
            /* We will compute with one extra limb */
            n++;
            /* ceil(log2(n+2)) takes into account the lost bits due to
               Mulders' short product */
            p = n * GMP_NUMB_BITS - MPFR_INT_CEIL_LOG2 (n + 2);
            /* Due to some nasty reasons we can have only 4 bits */
            MPFR_ASSERTD (MPFR_GET_PREC (a) <= p - 4);

            if (MPFR_LIKELY (k < 2*n))
              {
                tmp = MPFR_TMP_LIMBS_ALLOC (2 * n);
                tmp += 2*n-k; /* `tmp' still points to an area of `k' limbs */
              }
          }
        MPFR_LOG_MSG (("Use mpfr_mulhigh (%lu VS %lu)\n", MPFR_GET_PREC (a), p));
        /* Compute an approximation of the product of b and c */
        if (b != c)
          mpfr_mulhigh_n (tmp + k - 2 * n, bp, cp, n);
        else
          mpfr_sqrhigh_n (tmp + k - 2 * n, bp, n);
        /* now tmp[0]..tmp[k-1] contains the product of both mantissa,
           with tmp[k-1]>=2^(GMP_NUMB_BITS-2) */
        /* [VL] FIXME: This cannot be true: mpfr_mulhigh_n only
           depends on pointers and n. As k can be arbitrarily larger,
           the result cannot depend on k. And indeed, with GMP compiled
           with --enable-alloca=debug, valgrind was complaining, at
           least because MPFR_RNDRAW at the end tried to compute the
           sticky bit even when not necessary; this problem is fixed,
           but there's at least something wrong with the comment above. */
        b1 = tmp[k-1] >> (GMP_NUMB_BITS - 1); /* msb from the product */

        /* If the mantissas of b and c are uniformly distributed in (1/2, 1],
           then their product is in (1/4, 1/2] with probability 2*ln(2)-1
           ~ 0.386 and in [1/2, 1] with probability 2-2*ln(2) ~ 0.614 */
        if (MPFR_UNLIKELY (b1 == 0))
          /* Warning: the mpfr_mulhigh_n call above only surely affects
             tmp[k-n-1..k-1], thus we shift only those limbs */
          mpn_lshift (tmp + k - n - 1, tmp + k - n - 1, n + 1, 1);
        tmp += k - tn;
        MPFR_ASSERTD (MPFR_LIMB_MSB (tmp[tn-1]) != 0);

        /* if the most significant bit b1 is zero, we have only p-1 correct
           bits */
        if (rnd_mode == MPFR_RNDF && p + b1 - 1 >= MPFR_GET_PREC(a))
          {
            /* we can round */
          }
        else if (MPFR_UNLIKELY (!mpfr_round_p (tmp, tn, p + b1 - 1,
                            MPFR_GET_PREC(a) + (rnd_mode == MPFR_RNDN))))
          {
            tmp -= k - tn; /* tmp may have changed, FIX IT!!!!! */
            goto full_multiply;
          }
      }
    else
      {
      full_multiply:
        MPFR_LOG_MSG (("Use mpn_mul\n", 0));
        b1 = mpn_mul (tmp, MPFR_MANT (b), bn, MPFR_MANT (c), cn);

        /* now tmp[0]..tmp[k-1] contains the product of both mantissa,
           with tmp[k-1]>=2^(GMP_NUMB_BITS-2) */
        b1 >>= GMP_NUMB_BITS - 1; /* msb from the product */

        /* if the mantissas of b and c are uniformly distributed in (1/2, 1],
           then their product is in (1/4, 1/2] with probability 2*ln(2)-1
           ~ 0.386 and in [1/2, 1] with probability 2-2*ln(2) ~ 0.614 */
        tmp += k - tn;
        if (MPFR_UNLIKELY (b1 == 0))
          mpn_lshift (tmp, tmp, tn, 1); /* tn <= k, so no stack corruption */
      }

  ax2 = ax + (mpfr_exp_t) (b1 - 1);
  MPFR_RNDRAW (inexact, a, tmp, bq+cq, rnd_mode, sign, ax2++);
  MPFR_TMP_FREE (marker);
  MPFR_EXP (a) = ax2; /* Can't use MPFR_SET_EXP: Expo may be out of range */
  MPFR_SET_SIGN (a, sign);
  if (MPFR_UNLIKELY (ax2 > __gmpfr_emax))
    return mpfr_overflow (a, rnd_mode, sign);
  if (MPFR_UNLIKELY (ax2 < __gmpfr_emin))
    {
      /* In the rounding to the nearest mode, if the exponent of the exact
         result (i.e. before rounding, i.e. without taking cc into account)
         is < __gmpfr_emin - 1 or the exact result is a power of 2 (i.e. if
         both arguments are powers of 2), then round to zero. */
      if (rnd_mode == MPFR_RNDN
          && (ax + (mpfr_exp_t) b1 < __gmpfr_emin
              || (mpfr_powerof2_raw (b) && mpfr_powerof2_raw (c))))
        rnd_mode = MPFR_RNDZ;
      return mpfr_underflow (a, rnd_mode, sign);
    }
  MPFR_RET (inexact);
}
