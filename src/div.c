/* mpfr_div -- divide two floating-point numbers

Copyright 1999, 2001-2016 Free Software Foundation, Inc.
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

/* References:
   [1] Short Division of Long Integers, David Harvey and Paul Zimmermann,
       Proceedings of the 20th Symposium on Computer Arithmetic (ARITH-20),
       July 25-27, 2011, pages 7-14.
*/

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

#ifdef DEBUG2
#define mpfr_mpn_print(ap,n) mpfr_mpn_print3 (ap,n,MPFR_LIMB_ZERO)
static void
mpfr_mpn_print3 (mpfr_limb_ptr ap, mp_size_t n, mp_limb_t cy)
{
  mp_size_t i;
  for (i = 0; i < n; i++)
    printf ("+%lu*2^%lu", (unsigned long) ap[i], (unsigned long)
            (GMP_NUMB_BITS * i));
  if (cy)
    printf ("+2^%lu", (unsigned long) (GMP_NUMB_BITS * n));
  printf ("\n");
}
#endif

/* check if {ap, an} is zero */
static int
mpfr_mpn_cmpzero (mpfr_limb_ptr ap, mp_size_t an)
{
  MPFR_ASSERTD (an >= 0);
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
mpfr_mpn_cmp_aux (mpfr_limb_ptr ap, mp_size_t an,
                  mpfr_limb_ptr bp, mp_size_t bn, int extra)
{
  int cmp = 0;
  mp_size_t k;
  mp_limb_t bb;

  MPFR_ASSERTD (an >= 0);
  MPFR_ASSERTD (bn >= 0);
  MPFR_ASSERTD (extra == 0 || extra == 1);

  if (an >= bn)
    {
      k = an - bn;
      while (cmp == 0 && bn > 0)
        {
          bn --;
          bb = (extra) ? ((bp[bn+1] << (GMP_NUMB_BITS - 1)) | (bp[bn] >> 1))
            : bp[bn];
          cmp = (ap[k + bn] > bb) ? 1 : ((ap[k + bn] < bb) ? -1 : 0);
        }
      bb = (extra) ? bp[0] << (GMP_NUMB_BITS - 1) : MPFR_LIMB_ZERO;
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
          bb = (extra) ? ((bp[k+an+1] << (GMP_NUMB_BITS - 1)) | (bp[k+an] >> 1))
            : bp[k+an];
          if (ap[an] > bb)
            cmp = 1;
          else if (ap[an] < bb)
            cmp = -1;
        }
      while (cmp == 0 && k > 0)
        {
          k--;
          bb = (extra) ? ((bp[k+1] << (GMP_NUMB_BITS - 1)) | (bp[k] >> 1))
            : bp[k];
          cmp = (bb != MPFR_LIMB_ZERO) ? -1 : 0;
        }
      if (cmp == 0 && extra && (bp[0] & MPFR_LIMB_ONE))
        cmp = -1;
    }
  return cmp;
}

/* {ap, n} <- {ap, n} - {bp, n} >> extra - cy, with cy = 0 or 1.
   Return borrow out.
*/
static mp_limb_t
mpfr_mpn_sub_aux (mpfr_limb_ptr ap, mpfr_limb_ptr bp, mp_size_t n,
                  mp_limb_t cy, int extra)
{
  mp_limb_t bb, rp;

  MPFR_ASSERTD (cy <= 1);
  MPFR_ASSERTD (n >= 0);

  while (n--)
    {
      bb = (extra) ? ((bp[1] << (GMP_NUMB_BITS-1)) | (bp[0] >> 1)) : bp[0];
      rp = ap[0] - bb - cy;
      cy = (ap[0] < bb) || (cy && ~rp == MPFR_LIMB_ZERO) ?
        MPFR_LIMB_ONE : MPFR_LIMB_ZERO;
      ap[0] = rp;
      ap ++;
      bp ++;
    }
  MPFR_ASSERTD (cy <= 1);
  return cy;
}

/* For large precision, mpz_tdiv_q (which computes only quotient)
   is faster than mpn_divrem (which computes also the remainder).
   Unfortunately as of GMP 6.0.0 the corresponding mpn_div_q function
   is not in the public interface, thus we call mpz_tdiv_q.

   If this function succeeds in computing the correct rounding, return 1,
   and put the ternary value in inex.

   Otherwise return 0 (and inex is undefined).
*/
static int
mpfr_div_with_mpz_tdiv_q (mpfr_ptr q, mpfr_srcptr u, mpfr_srcptr v,
                          mpfr_rnd_t rnd_mode, int *inex)
{
  mpz_t qm, um, vm;
  mpfr_exp_t ue, ve;
  mpfr_prec_t qp = MPFR_PREC(q), wp = qp + GMP_NUMB_BITS;
  mp_size_t up, vp, k;
  int ok;

  mpz_init (qm);
  mpz_init (um);
  mpz_init (vm);

  ue = mpfr_get_z_2exp (um, u); /* u = um * 2^ue */
  ve = mpfr_get_z_2exp (vm, v); /* v = vm * 2^ve */

  vp = mpz_sizeinbase (vm, 2);
  if (vp > wp)
    {
      k = vp - wp; /* truncate k bits of vm */
      mpz_tdiv_q_2exp (vm, vm, k);
      ve += k;
      vp -= k;
    }

  /* we want about qp + GMP_NUMB_BITS bits of the quotient, thus um should
     have qp + GMP_NUMB_BITS more bits than vm */

  up = mpz_sizeinbase (um, 2);
  if (up > vp + wp)
    {
      k = up - (vp + wp); /* truncate k bits of um */
      mpz_tdiv_q_2exp (um, um, k);
      ue += k;
      up -= k;
    }
  else if (up < vp + wp) /* we need more bits */
    {
      k = (vp + wp) - up;
      mpz_mul_2exp (um, um, k);
      ue -= k;
      up += k;
    }

  /* now um has exactly wp more bits than vp */
  mpz_tdiv_q (qm, um, vm);
  /* qm has either wp or wp+1 bits, and we have:
     (a) um = u/2^ue*(1-tu) with tu=0 if no truncation of um,
                            and 0 <= tu < 2^(1-wp) otherwise;
     (b) vm = v/2^ve*(1-tv) with tv=0 if no truncation of vm,
                             and 0 <= tv < 2^(1-wp) otherwise;
     (c) um/vm - 1 < qm <= um/vm, thus qm = um/vm*(1-tq) with
         0 <= tw < 2^(1-wp) since um/vm >= 2^(wp-1)
     Altogether we have:
     q = u/v*2^(ve-ue)*(1-tu)/(1-tv)*(1-tq)
     Thus:
     u/v*2^(ve-ue)*(1-2^(2-wp)) < q < u/v*2^(ve-ue)*(1+2^(2-wp)).
     If q has wp bits, the error is less than 2^(wp-1)*2^(2-wp) <= 2.
     If q has wp+1 bits, the error is less than 2^wp*2^(2-wp) <= 4.
  */

  k = mpz_sizeinbase (qm, 2) - wp; /* 0 or 1 */
  /* Assume qm has wp bits (i.e. k=0) and a directed rounding: if the first
     set bit after position 1 has position less than GMP_NUMB_BITS, then
     subtracting 2 to qm will not change the bits beyond the GMP_NUMB_BITS
     low ones, thus we get correct rounding.
     For k=1, we need to start at position 2, and the first set bit has to be
     in posiiton less than GMP_NUMB_BITS+1.
     For rounding to nearest, the first set bit has to be in position less
     than GMP_NUMB_BITS-1 for k=0 (or less than GMP_NUMB_BITS for k=1).
  */
  if (mpz_scan1 (qm, k + 1) < GMP_NUMB_BITS + k - (rnd_mode == MPFR_RNDN) &&
      mpz_scan0 (qm, k + 1) < GMP_NUMB_BITS + k - (rnd_mode == MPFR_RNDN))
    {
      MPFR_SAVE_EXPO_DECL (expo);
      ok = 1;
      MPFR_SAVE_EXPO_MARK (expo);
      *inex = mpfr_set_z (q, qm, rnd_mode);
      MPFR_SAVE_EXPO_FREE (expo);
      /* if we got an underflow or overflow, the result is not valid */
      if (MPFR_IS_SINGULAR(q) || MPFR_EXP(q) == MPFR_EXT_EMIN ||
          MPFR_EXP(q) == MPFR_EXT_EMAX)
        ok =  0;
      MPFR_EXP(q) += ue - ve;
      *inex = mpfr_check_range (q, *inex, rnd_mode);
    }
  else
    ok = 0;

  mpz_clear (qm);
  mpz_clear (um);
  mpz_clear (vm);

  return ok;
}

/* special code for p=PREC(q) < GMP_NUMB_BITS, PREC(u), PREC(v) <= GMP_NUMB_BITS */
static int
mpfr_divsp1 (mpfr_ptr q, mpfr_srcptr u, mpfr_srcptr v, mpfr_rnd_t rnd_mode)
{
  mpfr_prec_t p = MPFR_GET_PREC(q);
  mpfr_limb_ptr up = MPFR_MANT(u);
  mpfr_limb_ptr vp = MPFR_MANT(v);
  mpfr_limb_ptr qp = MPFR_MANT(q);
  mpfr_exp_t qx = MPFR_GET_EXP(u) - MPFR_GET_EXP(v);
  mpfr_prec_t sh = GMP_NUMB_BITS - p;
  mp_limb_t h, rb, sb, mask = MPFR_LIMB_MASK(sh);

  if (up[0] >= vp[0])
    {
      udiv_qrnnd (h, sb, up[0] - vp[0], 0, vp[0]);
      /* Noting W = 2^GMP_NUMB_BITS, we have up[0]*W = (W + h) * vp[0] + sb,
         thus up[0]/vp[0] = 1 + h/W + sb/vp[0]/W, with 0 <= sb < vp[0]. */
      qx ++;
      sb |= h & 1;
      h = MPFR_LIMB_HIGHBIT | (h >> 1);
      rb = h & (MPFR_LIMB_ONE << (sh - 1));
      sb |= (h & mask) ^ rb;
      qp[0] = h & ~mask;
    }
  else
    {
      udiv_qrnnd (h, sb, up[0], 0, vp[0]);
      rb = h & (MPFR_LIMB_ONE << (sh - 1));
      sb |= (h & mask) ^ rb;
      qp[0] = h & ~mask;
    }
  
  MPFR_SIGN(q) = MPFR_MULT_SIGN (MPFR_SIGN (u), MPFR_SIGN (v));

  /* rounding */
  if (qx > __gmpfr_emax)
    return mpfr_overflow (q, rnd_mode, MPFR_SIGN(q));

  /* Warning: underflow should be checked *after* rounding, thus when rounding
     away and when q > 0.111...111*2^(emin-1), or when rounding to nearest and
     q >= 0.111...111[1]*2^(emin-1), there is no underflow. */
  if (qx < __gmpfr_emin)
    {
      /* for RNDN, mpfr_underflow always rounds away, thus for |q|<=2^(emin-2)
         we have to change to RNDZ */
      if (rnd_mode == MPFR_RNDN)
        {
          if ((qx == __gmpfr_emin - 1) && (qp[0] == ~mask) && rb)
            goto rounding; /* no underflow */
          if (qx < __gmpfr_emin - 1 || (qp[0] == MPFR_LIMB_HIGHBIT && sb == 0))
            rnd_mode = MPFR_RNDZ;
        }
      else if (!MPFR_IS_LIKE_RNDZ(rnd_mode, MPFR_IS_NEG (q)))
        {
          if ((qx == __gmpfr_emin - 1) && (qp[0] == ~mask) && (rb | sb))
            goto rounding; /* no underflow */
        }
      return mpfr_underflow (q, rnd_mode, MPFR_SIGN(q));
    }

 rounding:
  MPFR_EXP (q) = qx; /* Don't use MPFR_SET_EXP since qx might be < __gmpfr_emin
                        in the cases "goto rounding" above. */
  if (rb == 0 && sb == 0)
    {
      MPFR_ASSERTD(qx >= __gmpfr_emin);
      return 0; /* idem than MPFR_RET(0) but faster */
    }
  else if (rnd_mode == MPFR_RNDN)
    {
      if (rb == 0 || (rb && sb == 0 &&
                      (qp[0] & (MPFR_LIMB_ONE << sh)) == 0))
        goto truncate;
      else
        goto add_one_ulp;
    }
  else if (MPFR_IS_LIKE_RNDZ(rnd_mode, MPFR_IS_NEG(q)))
    {
    truncate:
      MPFR_ASSERTD(qx >= __gmpfr_emin);
      MPFR_RET(-MPFR_SIGN(q));
    }
  else /* round away from zero */
    {
    add_one_ulp:
      qp[0] += MPFR_LIMB_ONE << sh;
      if (qp[0] == 0)
        {
          qp[0] = MPFR_LIMB_HIGHBIT;
          if (MPFR_UNLIKELY(qx + 1 > __gmpfr_emax))
            return mpfr_overflow (q, rnd_mode, MPFR_SIGN(q));
          MPFR_ASSERTD(qx + 1 <= __gmpfr_emax);
          MPFR_ASSERTD(qx + 1 >= __gmpfr_emin);
          MPFR_SET_EXP (q, qx + 1);
        }
      MPFR_RET(MPFR_SIGN(q));
    }
}

MPFR_HOT_FUNCTION_ATTR int
mpfr_div (mpfr_ptr q, mpfr_srcptr u, mpfr_srcptr v, mpfr_rnd_t rnd_mode)
{
  mp_size_t q0size, usize, vsize;
  mp_size_t qsize; /* number of limbs wanted for the computed quotient */
  mp_size_t qqsize;
  mp_size_t k;
  mpfr_limb_ptr q0p, qp;
  mpfr_limb_ptr up, vp;
  mpfr_limb_ptr ap;
  mpfr_limb_ptr bp;
  mp_limb_t qh;
  mp_limb_t sticky_u, sticky_v;
  mp_limb_t low_u;
  mp_limb_t sticky;
  mp_limb_t sticky3;
  mp_limb_t round_bit;
  mpfr_exp_t qexp;
  int sign_quotient;
  int extra_bit;
  int sh, sh2;
  int inex;
  int like_rndz;
  MPFR_TMP_DECL(marker);

  MPFR_LOG_FUNC (
    ("u[%Pu]=%.*Rg v[%Pu]=%.*Rg rnd=%d",
     mpfr_get_prec(u), mpfr_log_prec, u,
     mpfr_get_prec (v),mpfr_log_prec, v, rnd_mode),
    ("q[%Pu]=%.*Rg inexact=%d", mpfr_get_prec(q), mpfr_log_prec, q, inex));

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
              MPFR_ASSERTD (! MPFR_IS_INF (u));
              MPFR_SET_INF(q);
              MPFR_SET_DIVBY0 ();
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

  if (MPFR_GET_PREC(q) < GMP_NUMB_BITS && MPFR_GET_PREC(u) <= GMP_NUMB_BITS
      && MPFR_GET_PREC(v) <= GMP_NUMB_BITS)
    return mpfr_divsp1 (q, u, v, rnd_mode);

  q0size = MPFR_LIMB_SIZE(q); /* number of limbs of destination */
  usize = MPFR_LIMB_SIZE(u);
  vsize = MPFR_LIMB_SIZE(v);
  q0p = MPFR_MANT(q);
  up = MPFR_MANT(u);
  vp = MPFR_MANT(v);
  sticky_u = MPFR_LIMB_ZERO;
  sticky_v = MPFR_LIMB_ZERO;
  round_bit = MPFR_LIMB_ZERO;

  /**************************************************************************
   *                                                                        *
   *              End of the part concerning special values.                *
   *                                                                        *
   **************************************************************************/

  /* when the divisor has one limb, we can use mpfr_div_ui, which should be
     faster, assuming there is no intermediate overflow or underflow.
     The divisor interpreted as an integer satisfies
     2^(GMP_NUMB_BITS-1) <= vm < 2^GMP_NUMB_BITS, thus the quotient
     satisfies 2^(EXP(u)-1-GMP_NUMB_BITS) < u/vm < 2^(EXP(u)-GMP_NUMB_BITS+1)
     and its exponent is either EXP(u)-GMP_NUMB_BITS or one more. */
  if (vsize <= 1 && __gmpfr_emin <= MPFR_EXP(u) - GMP_NUMB_BITS
      && MPFR_EXP(u) - GMP_NUMB_BITS + 1 <= __gmpfr_emax
      && vp[0] <= ULONG_MAX)
    {
      mpfr_exp_t exp_v = MPFR_EXP(v); /* save it in case q=v */
      if (MPFR_IS_POS (v))
        inex = mpfr_div_ui (q, u, vp[0], rnd_mode);
      else
        {
          inex = -mpfr_div_ui (q, u, vp[0], MPFR_INVERT_RND(rnd_mode));
          MPFR_CHANGE_SIGN(q);
        }
      /* q did not under/overflow */
      MPFR_EXP(q) -= exp_v;
      /* The following test is needed, otherwise the next addition
         on the exponent may overflow, e.g. when dividing the
         largest finite MPFR number by the smallest positive one. */
      if (MPFR_UNLIKELY (MPFR_EXP(q) > __gmpfr_emax - GMP_NUMB_BITS))
        return mpfr_overflow (q, rnd_mode, MPFR_SIGN(q));
      MPFR_EXP(q) += GMP_NUMB_BITS;
      return mpfr_check_range (q, inex, rnd_mode);
    }

  /* for large precisions, try using truncated division first */
  if (q0size >= 32 && mpfr_div_with_mpz_tdiv_q (q, u, v, rnd_mode, &inex))
    return inex;

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
      mp_size_t l;

      k = usize - 1;
      l = vsize - 1;
      while (k != 0 && l != 0 && up[--k] == vp[--l]);
      /* now k=0 or l=0 or up[k] != vp[l] */
      if (up[k] != vp[l])
        extra_bit = (up[k] > vp[l]);
      /* now up[k] = vp[l], thus either k=0 or l=0 */
      else if (l == 0) /* no more divisor limb */
        extra_bit = 1;
      else /* k=0: no more dividend limb */
        extra_bit = mpfr_mpn_cmpzero (vp, l) == 0;
    }
#ifdef DEBUG
  printf ("extra_bit=%d\n", extra_bit);
#endif

  /* set exponent */
  qexp = MPFR_GET_EXP (u) - MPFR_GET_EXP (v) + extra_bit;

  /* sh is the number of zero bits in the low limb of the quotient */
  MPFR_UNSIGNED_MINUS_MODULO(sh, MPFR_PREC(q));

  like_rndz = rnd_mode == MPFR_RNDZ ||
    rnd_mode == (sign_quotient < 0 ? MPFR_RNDU : MPFR_RNDD);

  /**************************************************************************
   *                                                                        *
   *       We first try Mulders' short division (for large operands)        *
   *                                                                        *
   **************************************************************************/

  if (MPFR_UNLIKELY(q0size >= MPFR_DIV_THRESHOLD &&
                    vsize >= MPFR_DIV_THRESHOLD))
    {
      mp_size_t n = q0size + 1; /* we will perform a short (2n)/n division */
      mpfr_limb_ptr ap, bp, qp;
      mpfr_prec_t p;

      /* since Mulders' short division clobbers the dividend, we have to
         copy it */
      ap = MPFR_TMP_LIMBS_ALLOC (n + n);
      if (usize >= n + n) /* truncate the dividend */
        MPN_COPY(ap, up + usize - (n + n), n + n);
      else                /* zero-pad the dividend */
        {
          MPN_COPY(ap + (n + n) - usize, up, usize);
          MPN_ZERO(ap, (n + n) - usize);
        }

      if (vsize >= n) /* truncate the divisor */
        bp = vp + vsize - n;
      else            /* zero-pad the divisor */
        {
          bp = MPFR_TMP_LIMBS_ALLOC (n);
          MPN_COPY(bp + n - vsize, vp, vsize);
          MPN_ZERO(bp, n - vsize);
        }

      qp = MPFR_TMP_LIMBS_ALLOC (n);
      qh = mpfr_divhigh_n (qp, ap, bp, n);
      MPFR_ASSERTD (qh == 0 || qh == 1);
      /* in all cases, the error is at most (2n+2) ulps on qh*B^n+{qp,n},
         cf algorithms.tex */

      p = n * GMP_NUMB_BITS - MPFR_INT_CEIL_LOG2 (2 * n + 2);
      /* If rnd=RNDN, we need to be able to round with a directed rounding
         and one more bit. */
      if (qh == 1)
        {
          mpn_rshift (qp, qp, n, 1);
          qp[n - 1] |= MPFR_LIMB_HIGHBIT;
        }
      if (MPFR_LIKELY (mpfr_round_p (qp, n, p,
                                     MPFR_PREC(q) + (rnd_mode == MPFR_RNDN))))
        {
          /* we can round correctly whatever the rounding mode */
          MPN_COPY (q0p, qp + 1, q0size);
          q0p[0] &= ~MPFR_LIMB_MASK(sh); /* put to zero low sh bits */

          if (rnd_mode == MPFR_RNDN) /* round to nearest */
            {
              /* we know we can round, thus we are never in the even rule case:
                 if the round bit is 0, we truncate
                 if the round bit is 1, we add 1 */
              if (sh > 0)
                round_bit = (qp[1] >> (sh - 1)) & 1;
              else
                round_bit = qp[0] >> (GMP_NUMB_BITS - 1);
              /* TODO: add value coverage tests in tdiv to check that
                 we reach this part with different values of qh and
                 round_bit (4 cases). */
              if (round_bit == 0)
                {
                  inex = -1;
                  goto truncate;
                }
              else /* round_bit = 1 */
                goto add_one_ulp;
            }
          else if (! like_rndz) /* round away */
            goto add_one_ulp;
          else /* round to zero: nothing to do */
            {
              inex = -1;
              goto truncate;
            }
        }
    }

  /**************************************************************************
   *                                                                        *
   *     Mulders' short division failed: we revert to integer division      *
   *                                                                        *
   **************************************************************************/

  if (MPFR_UNLIKELY(rnd_mode == MPFR_RNDN && sh == 0))
    { /* we compute the quotient with one more limb, in order to get
         the round bit in the quotient, and the remainder only contains
         sticky bits */
      qsize = q0size + 1;
      /* need to allocate memory for the quotient */
      qp = MPFR_TMP_LIMBS_ALLOC (qsize);
    }
  else
    {
      qsize = q0size;
      qp = q0p; /* directly put the quotient in the destination */
    }
  qqsize = qsize + qsize;

  /* prepare the dividend */
  ap = MPFR_TMP_LIMBS_ALLOC (qqsize);
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
      sticky_u = sticky_u || mpfr_mpn_cmpzero (up, k);
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
          bp = MPFR_TMP_LIMBS_ALLOC (qsize);
          MPN_COPY(bp, vp, vsize);
        }
      sticky_v = sticky_v || mpfr_mpn_cmpzero (vp, k);
      k = 0;
    }
  else /* vsize < qsize: small divisor case */
    {
      bp = vp;
      k = qsize - vsize;
    }

  /**************************************************************************
   *                                                                        *
   *  Here we perform the real division of {ap+k,qqsize-k} by {bp,qsize-k}  *
   *                                                                        *
   **************************************************************************/

  /* if Mulders' short division failed, we revert to division with remainder */
  qh = mpn_divrem (qp, 0, ap + k, qqsize - k, bp, qsize - k);
  /* warning: qh may be 1 if u1 == v1, but u < v */
#ifdef DEBUG2
  printf ("q="); mpfr_mpn_print (qp, qsize);
  printf ("r="); mpfr_mpn_print (ap, qsize);
#endif

  k = qsize;
  sticky_u = sticky_u || mpfr_mpn_cmpzero (ap, k);

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
  else /* qsize = q0size + 1: only happens when rnd_mode=MPFR_RNDN and sh=0 */
    {
      MPN_COPY (q0p, qp + 1, q0size);
      sticky3 = qp[0];
      sh2 = GMP_NUMB_BITS;
    }
  qp[0] ^= sticky3;
  /* sticky3 contains the truncated bits from the quotient,
     including the round bit, and 1 <= sh2 <= GMP_NUMB_BITS
     is the number of bits in sticky3 */
  inex = (sticky != MPFR_LIMB_ZERO) || (sticky3 != MPFR_LIMB_ZERO);
#ifdef DEBUG
  printf ("sticky=%lu sticky3=%lu inex=%d\n",
          (unsigned long) sticky, (unsigned long) sticky3, inex);
#endif

  /* to round, we distinguish two cases:
     (a) vsize <= qsize: we used the full divisor
     (b) vsize > qsize: the divisor was truncated
  */

#ifdef DEBUG
  printf ("vsize=%lu qsize=%lu\n",
          (unsigned long) vsize, (unsigned long) qsize);
#endif
  if (MPFR_LIKELY(vsize <= qsize)) /* use the full divisor */
    {
      if (MPFR_LIKELY(rnd_mode == MPFR_RNDN))
        {
          round_bit = sticky3 & (MPFR_LIMB_ONE << (sh2 - 1));
          sticky = (sticky3 ^ round_bit) | sticky_u;
        }
      else if (like_rndz || inex == 0)
        sticky = (inex == 0) ? MPFR_LIMB_ZERO : MPFR_LIMB_ONE;
      else  /* round away from zero */
        sticky = MPFR_LIMB_ONE;
      goto case_1;
    }
  else /* vsize > qsize: need to truncate the divisor */
    {
      if (inex == 0)
        goto truncate;
      else
        {
          /* We know the estimated quotient is an upper bound of the exact
             quotient (with rounding toward zero), with a difference of at
             most 2 in qp[0].
             Thus we can round except when sticky3 is 000...000 or 000...001
             for directed rounding, and 100...000 or 100...001 for rounding
             to nearest. (For rounding to nearest, we cannot determine the
             inexact flag for 000...000 or 000...001.)
          */
          mp_limb_t sticky3orig = sticky3;
          if (rnd_mode == MPFR_RNDN)
            {
              round_bit = sticky3 & (MPFR_LIMB_ONE << (sh2 - 1));
              sticky3   = sticky3 ^ round_bit;
#ifdef DEBUG
              printf ("rb=%lu sb=%lu\n",
                      (unsigned long) round_bit, (unsigned long) sticky3);
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
              mpfr_limb_ptr sp;
              int cmp_s_r;
              mp_limb_t qh2;

              sp = MPFR_TMP_LIMBS_ALLOC (vsize);
              k = vsize - qsize;
              /* sp <- {qp, qsize} * {vp, vsize-qsize} */
              qp[0] ^= sticky3orig; /* restore original quotient */
              if (qsize >= k)
                mpn_mul (sp, qp, qsize, vp, k);
              else
                mpn_mul (sp, vp, k, qp, qsize);
              if (qh)
                qh2 = mpn_add_n (sp + qsize, sp + qsize, vp, k);
              else
                qh2 = MPFR_LIMB_ZERO;
              qp[0] ^= sticky3orig; /* restore truncated quotient */

              /* compare qh2 + {sp, k + qsize} to {ap, qsize} + low(u) */
              cmp_s_r = (qh2 != 0) ? 1 : mpn_cmp (sp + k, ap, qsize);
              if (cmp_s_r == 0) /* compare {sp, k} and low(u) */
                {
                  cmp_s_r = (usize >= qqsize) ?
                    mpfr_mpn_cmp_aux (sp, k, up, usize - qqsize, extra_bit) :
                    mpfr_mpn_cmpzero (sp, k);
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
              else /* cmp_s_r > 0, quotient is < q1: to determine if it is
                      in [q1-2,q1-1] or in [q1-1,q1], we need to subtract
                      the low part u0 of the dividend u0 from q*v0 */
                {
                  mp_limb_t cy = MPFR_LIMB_ZERO;

                  /* subtract low(u)>>extra_bit if non-zero */
                  if (qh2 != 0) /* whatever the value of {up, m + k}, it
                                   will be smaller than qh2 + {sp, k} */
                    cmp_s_r = 1;
                  else
                    {
                      if (low_u != MPFR_LIMB_ZERO)
                        {
                          mp_size_t m;
                          l = usize - qqsize; /* number of low limbs in u */
                          m = (l > k) ? l - k : 0;
                          cy = (extra_bit) ?
                            (up[m] & MPFR_LIMB_ONE) : MPFR_LIMB_ZERO;
                          if (l >= k) /* u0 has more limbs than s:
                                         first look if {up, m} is not zero,
                                         and compare {sp, k} and {up + m, k} */
                            {
                              cy = cy || mpfr_mpn_cmpzero (up, m);
                              low_u = cy;
                              cy = mpfr_mpn_sub_aux (sp, up + m, k,
                                                     cy, extra_bit);
                            }
                          else /* l < k: s has more limbs than u0 */
                            {
                              low_u = MPFR_LIMB_ZERO;
                              if (cy != MPFR_LIMB_ZERO)
                                cy = mpn_sub_1 (sp + k - l - 1, sp + k - l - 1,
                                                1, MPFR_LIMB_HIGHBIT);
                              cy = mpfr_mpn_sub_aux (sp + k - l, up, l,
                                                     cy, extra_bit);
                            }
                        }
                      MPFR_ASSERTD (cy <= 1);
                      cy = mpn_sub_1 (sp + k, sp + k, qsize, cy);
                      /* subtract r */
                      cy += mpn_sub_n (sp + k, sp + k, ap, qsize);
                      MPFR_ASSERTD (cy <= 1);
                      /* now compare {sp, ssize} to v */
                      cmp_s_r = mpn_cmp (sp, vp, vsize);
                      if (cmp_s_r == 0 && low_u != MPFR_LIMB_ZERO)
                        cmp_s_r = 1; /* since in fact we subtracted
                                        less than 1 */
                    }
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
                          if (rnd_mode == MPFR_RNDN ||
                              (! like_rndz && inex != 0))
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
                      /* if rnd=MPFR_RNDN, the result is q1 when
                         q1-2 >= q1-2^(sh-1), i.e. sh >= 2,
                         otherwise (sh=1) it is q1-2 */
                      if (rnd_mode == MPFR_RNDN) /* sh > 0 */
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
                      else if (like_rndz)
                        {
                          /* the result is down(q1-2), i.e. subtract one
                             ulp if sh > 0, and two ulps if sh=0 */
                          inex = -1;
                          if (sh > 0)
                            goto sub_one_ulp;
                          else
                            goto sub_two_ulp;
                        }
                      /* if round away from zero, the result is up(q1-1),
                         which is q1 unless sh = 0, where it is q1-1 */
                      else
                        {
                          inex = 1;
                          if (sh > 0)
                            goto truncate_check_qh;
                          else /* sh = 0 */
                            goto sub_one_ulp;
                        }
                    }
                }
            }
        }
    }

 case_1: /* quotient is in [q1, q1+1),
            round_bit is the round_bit (0 for directed rounding),
            sticky the sticky bit */
  if (like_rndz || (round_bit == MPFR_LIMB_ZERO && sticky == MPFR_LIMB_ZERO))
    {
      inex = round_bit == MPFR_LIMB_ZERO && sticky == MPFR_LIMB_ZERO ? 0 : -1;
      goto truncate;
    }
  else if (rnd_mode == MPFR_RNDN) /* sticky <> 0 or round <> 0 */
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
  else /* round away from zero, sticky <> 0 */
    goto add_one_ulp; /* with inex=1 */

 sub_two_ulp:
  /* we cannot subtract MPFR_LIMB_MPFR_LIMB_ONE << (sh+1) since this is
     undefined for sh = GMP_NUMB_BITS */
  qh -= mpn_sub_1 (q0p, q0p, q0size, MPFR_LIMB_ONE << sh);
  /* go through */

 sub_one_ulp:
  qh -= mpn_sub_1 (q0p, q0p, q0size, MPFR_LIMB_ONE << sh);
  /* go through truncate_check_qh */

 truncate_check_qh:
  if (qh)
    {
      if (MPFR_LIKELY (qexp < MPFR_EXP_MAX))
        qexp ++;
      /* else qexp is now incorrect, but one will still get an overflow */
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
      if (MPFR_LIKELY (qexp < MPFR_EXP_MAX))
        qexp ++;
      /* else qexp is now incorrect, but one will still get an overflow */
      q0p[q0size - 1] = MPFR_LIMB_HIGHBIT;
    }

 truncate: /* inex already set */

  MPFR_TMP_FREE(marker);

  /* check for underflow/overflow */
  if (MPFR_UNLIKELY(qexp > __gmpfr_emax))
    return mpfr_overflow (q, rnd_mode, sign_quotient);
  else if (MPFR_UNLIKELY(qexp < __gmpfr_emin))
    {
      if (rnd_mode == MPFR_RNDN && ((qexp < __gmpfr_emin - 1) ||
                                   (inex >= 0 && mpfr_powerof2_raw (q))))
        rnd_mode = MPFR_RNDZ;
      return mpfr_underflow (q, rnd_mode, sign_quotient);
    }
  MPFR_SET_EXP(q, qexp);

  inex *= sign_quotient;
  MPFR_RET (inex);
}
