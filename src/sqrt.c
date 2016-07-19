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

/* the following routine is copied from GMP 6.1.1, file mpn/generic/sqrtrem.c */

/* invsqrttab[x] ~ sqrt(2^25/x)-2^8 for 2^7 <= x < 2^9 */
static const unsigned char invsqrttab[384] = /* The common 0x100 was removed */
{
  0xff,0xfd,0xfb,0xf9,0xf7,0xf5,0xf3,0xf2, /* sqrt(1/80)..sqrt(1/87) */
  0xf0,0xee,0xec,0xea,0xe9,0xe7,0xe5,0xe4, /* sqrt(1/88)..sqrt(1/8f) */
  0xe2,0xe0,0xdf,0xdd,0xdb,0xda,0xd8,0xd7, /* sqrt(1/90)..sqrt(1/97) */
  0xd5,0xd4,0xd2,0xd1,0xcf,0xce,0xcc,0xcb, /* sqrt(1/98)..sqrt(1/9f) */
  0xc9,0xc8,0xc6,0xc5,0xc4,0xc2,0xc1,0xc0, /* sqrt(1/a0)..sqrt(1/a7) */
  0xbe,0xbd,0xbc,0xba,0xb9,0xb8,0xb7,0xb5, /* sqrt(1/a8)..sqrt(1/af) */
  0xb4,0xb3,0xb2,0xb0,0xaf,0xae,0xad,0xac, /* sqrt(1/b0)..sqrt(1/b7) */
  0xaa,0xa9,0xa8,0xa7,0xa6,0xa5,0xa4,0xa3, /* sqrt(1/b8)..sqrt(1/bf) */
  0xa2,0xa0,0x9f,0x9e,0x9d,0x9c,0x9b,0x9a, /* sqrt(1/c0)..sqrt(1/c7) */
  0x99,0x98,0x97,0x96,0x95,0x94,0x93,0x92, /* sqrt(1/c8)..sqrt(1/cf) */
  0x91,0x90,0x8f,0x8e,0x8d,0x8c,0x8c,0x8b, /* sqrt(1/d0)..sqrt(1/d7) */
  0x8a,0x89,0x88,0x87,0x86,0x85,0x84,0x83, /* sqrt(1/d8)..sqrt(1/df) */
  0x83,0x82,0x81,0x80,0x7f,0x7e,0x7e,0x7d, /* sqrt(1/e0)..sqrt(1/e7) */
  0x7c,0x7b,0x7a,0x79,0x79,0x78,0x77,0x76, /* sqrt(1/e8)..sqrt(1/ef) */
  0x76,0x75,0x74,0x73,0x72,0x72,0x71,0x70, /* sqrt(1/f0)..sqrt(1/f7) */
  0x6f,0x6f,0x6e,0x6d,0x6d,0x6c,0x6b,0x6a, /* sqrt(1/f8)..sqrt(1/ff) */
  0x6a,0x69,0x68,0x68,0x67,0x66,0x66,0x65, /* sqrt(1/100)..sqrt(1/107) */
  0x64,0x64,0x63,0x62,0x62,0x61,0x60,0x60, /* sqrt(1/108)..sqrt(1/10f) */
  0x5f,0x5e,0x5e,0x5d,0x5c,0x5c,0x5b,0x5a, /* sqrt(1/110)..sqrt(1/117) */
  0x5a,0x59,0x59,0x58,0x57,0x57,0x56,0x56, /* sqrt(1/118)..sqrt(1/11f) */
  0x55,0x54,0x54,0x53,0x53,0x52,0x52,0x51, /* sqrt(1/120)..sqrt(1/127) */
  0x50,0x50,0x4f,0x4f,0x4e,0x4e,0x4d,0x4d, /* sqrt(1/128)..sqrt(1/12f) */
  0x4c,0x4b,0x4b,0x4a,0x4a,0x49,0x49,0x48, /* sqrt(1/130)..sqrt(1/137) */
  0x48,0x47,0x47,0x46,0x46,0x45,0x45,0x44, /* sqrt(1/138)..sqrt(1/13f) */
  0x44,0x43,0x43,0x42,0x42,0x41,0x41,0x40, /* sqrt(1/140)..sqrt(1/147) */
  0x40,0x3f,0x3f,0x3e,0x3e,0x3d,0x3d,0x3c, /* sqrt(1/148)..sqrt(1/14f) */
  0x3c,0x3b,0x3b,0x3a,0x3a,0x39,0x39,0x39, /* sqrt(1/150)..sqrt(1/157) */
  0x38,0x38,0x37,0x37,0x36,0x36,0x35,0x35, /* sqrt(1/158)..sqrt(1/15f) */
  0x35,0x34,0x34,0x33,0x33,0x32,0x32,0x32, /* sqrt(1/160)..sqrt(1/167) */
  0x31,0x31,0x30,0x30,0x2f,0x2f,0x2f,0x2e, /* sqrt(1/168)..sqrt(1/16f) */
  0x2e,0x2d,0x2d,0x2d,0x2c,0x2c,0x2b,0x2b, /* sqrt(1/170)..sqrt(1/177) */
  0x2b,0x2a,0x2a,0x29,0x29,0x29,0x28,0x28, /* sqrt(1/178)..sqrt(1/17f) */
  0x27,0x27,0x27,0x26,0x26,0x26,0x25,0x25, /* sqrt(1/180)..sqrt(1/187) */
  0x24,0x24,0x24,0x23,0x23,0x23,0x22,0x22, /* sqrt(1/188)..sqrt(1/18f) */
  0x21,0x21,0x21,0x20,0x20,0x20,0x1f,0x1f, /* sqrt(1/190)..sqrt(1/197) */
  0x1f,0x1e,0x1e,0x1e,0x1d,0x1d,0x1d,0x1c, /* sqrt(1/198)..sqrt(1/19f) */
  0x1c,0x1b,0x1b,0x1b,0x1a,0x1a,0x1a,0x19, /* sqrt(1/1a0)..sqrt(1/1a7) */
  0x19,0x19,0x18,0x18,0x18,0x18,0x17,0x17, /* sqrt(1/1a8)..sqrt(1/1af) */
  0x17,0x16,0x16,0x16,0x15,0x15,0x15,0x14, /* sqrt(1/1b0)..sqrt(1/1b7) */
  0x14,0x14,0x13,0x13,0x13,0x12,0x12,0x12, /* sqrt(1/1b8)..sqrt(1/1bf) */
  0x12,0x11,0x11,0x11,0x10,0x10,0x10,0x0f, /* sqrt(1/1c0)..sqrt(1/1c7) */
  0x0f,0x0f,0x0f,0x0e,0x0e,0x0e,0x0d,0x0d, /* sqrt(1/1c8)..sqrt(1/1cf) */
  0x0d,0x0c,0x0c,0x0c,0x0c,0x0b,0x0b,0x0b, /* sqrt(1/1d0)..sqrt(1/1d7) */
  0x0a,0x0a,0x0a,0x0a,0x09,0x09,0x09,0x09, /* sqrt(1/1d8)..sqrt(1/1df) */
  0x08,0x08,0x08,0x07,0x07,0x07,0x07,0x06, /* sqrt(1/1e0)..sqrt(1/1e7) */
  0x06,0x06,0x06,0x05,0x05,0x05,0x04,0x04, /* sqrt(1/1e8)..sqrt(1/1ef) */
  0x04,0x04,0x03,0x03,0x03,0x03,0x02,0x02, /* sqrt(1/1f0)..sqrt(1/1f7) */
  0x02,0x02,0x01,0x01,0x01,0x01,0x00,0x00  /* sqrt(1/1f8)..sqrt(1/1ff) */
};

#if GMP_NUMB_BITS > 32
#define MAGIC 0x10000000000	/* 0xffe7debbfc < MAGIC < 0x232b1850f410 */
#else
#define MAGIC 0x100000		/* 0xfee6f < MAGIC < 0x29cbc8 */
#endif

static mp_limb_t
mpn_sqrtrem1 (mp_ptr rp, mp_limb_t a0)
{
#if GMP_NUMB_BITS > 32
  mp_limb_t a1;
#endif
  mp_limb_t x0, t2, t, x2;
  unsigned abits;

  /* Use Newton iterations for approximating 1/sqrt(a) instead of sqrt(a),
     since we can do the former without division.  As part of the last
     iteration convert from 1/sqrt(a) to sqrt(a).  */

  abits = a0 >> (GMP_LIMB_BITS - 1 - 8);	/* extract bits for table lookup */
  x0 = 0x100 | invsqrttab[abits - 0x80];	/* initial 1/sqrt(a) */

  /* x0 is now an 8 bits approximation of 1/sqrt(a0) */

#if GMP_NUMB_BITS > 32
  a1 = a0 >> (GMP_LIMB_BITS - 1 - 32);
  t = (mp_limb_signed_t) ((0x2000000000000) - 0x30000 - a1 * x0 * x0) >> 16;
  x0 = (x0 << 16) + ((mp_limb_signed_t) (x0 * t) >> (16+2));

  /* x0 is now a 16 bits approximation of 1/sqrt(a0) */

  t2 = x0 * (a0 >> (32-8));
  t = t2 >> 25;
  t = ((mp_limb_signed_t) ((a0 << 14) - t * t - MAGIC) >> (32-8));
  x0 = t2 + ((mp_limb_signed_t) (x0 * t) >> 15);
  x0 >>= 32;
#else
  t2 = x0 * (a0 >> (16-8));
  t = t2 >> 13;
  t = ((mp_limb_signed_t) ((a0 << 6) - t * t - MAGIC) >> (16-8));
  x0 = t2 + ((mp_limb_signed_t) (x0 * t) >> 7);
  x0 >>= 16;
#endif

  /* x0 is now a full limb approximation of sqrt(a0) */

  x2 = x0 * x0;
  if (x2 + 2*x0 <= a0 - 1)
    {
      x2 += 2*x0 + 1;
      x0++;
    }

  *rp = a0 - x2;
  return x0;
}
/* end of code copied from GMP 6.1.1 */

/* return s such that a*2^GMP_NUMB_BITS = s^2 + r with 0 <= r <= 2*s,
   and *inexact = 0 iff r <> 0. */
static mp_limb_t
mpfr_sqrt_limb2 (mp_limb_t a, mp_limb_t *inexact)
{
  mp_limb_t s, r, q, h, l;

  s = mpn_sqrtrem1 (&r, a);
  /* now a = s^2 + r[1], with 0 <= r[1] <= 2*s, and s < 2^32 */
  q = (r << 31) / s;
  s = (s << 32) + q;
  umul_ppmm (h, l, s, s);
  if (h > a || (h == a && l > 0))
    {
      h -= (l < 2 * s - 1);
      l -= 2 * s - 1;
      s--;
    }
  *inexact = (h != a) || (l != 0);
  return s;
}

/* Special code for prec(r), prec(u) < GMP_NUMB_BITS. */
static int
mpfr_sqrt1 (mpfr_ptr r, mpfr_srcptr u, mpfr_rnd_t rnd_mode)
{
  mpfr_prec_t p = MPFR_GET_PREC(r);
  mpfr_prec_t exp_u = MPFR_EXP(u), exp_r, sh = GMP_NUMB_BITS - p;
  mp_limb_t u0, r0, rb, sb, mask;
  mpfr_limb_ptr rp = MPFR_MANT(r);

  u0 = MPFR_MANT(u)[0];
  if (exp_u & 1)
    {
      u0 >>= 1;
      exp_u ++;
    }
  exp_r = exp_u >> 1;

  if (p < GMP_NUMB_BITS / 2)
    r0 = mpn_sqrtrem1 (&sb, u0) << (GMP_NUMB_BITS / 2);
  else
    r0 = mpfr_sqrt_limb2 (u0, &sb);

  rb = r0 & (MPFR_LIMB_ONE << (sh - 1));
  mask = MPFR_LIMB_MASK(sh);
  sb |= (r0 & mask) ^ rb;
  rp[0] = r0 & ~mask;

  /* rounding */
  if (exp_r > __gmpfr_emax)
    return mpfr_overflow (r, rnd_mode, 1);

  /* See comments in mpfr_divsp1 */
  if (exp_r < __gmpfr_emin)
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

  if (MPFR_GET_PREC (r) < GMP_NUMB_BITS && MPFR_GET_PREC (u) < GMP_NUMB_BITS)
    return mpfr_sqrt1 (r, u, rnd_mode);

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
                     : ~MPFR_LIMB_ZERO);
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
