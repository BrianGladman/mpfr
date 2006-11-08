/* mpfr_get_decimal64 -- convert a multiple precision floating-point number
                         to a IEEE 754r decimal64 float

See http://gcc.gnu.org/ml/gcc/2006-06/msg00691.html
and http://gcc.gnu.org/onlinedocs/gcc/Decimal-Float.html.

Copyright 2006 Free Software Foundation, Inc.

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
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <stdio.h>  /* DEBUG */
#include <string.h> /* for strlen */
#include "mpfr-impl.h"

#if MPFR_WANT_DECIMAL_FLOATS

#ifdef DPD_FORMAT
static _Decimal64
bid_to_dpd (_Decimal64 d)
{
  union ieee_double_extract x;
  union ieee_double_decimal64 y;
  mp_limb_t h, l;
  unsigned long d0, d1, d2, d3, d4, d5; /* declets: 0 <= d0..d5 < 999 */

  y.d64 = d;
  x.d = y.d;
  G = x.s.exp << 2; /* x.s.exp is the biased exponent, in [0, 2048) */
#if BITS_PER_MP_LIMB == 32
  G |= x.s.manh >> 18; /* manh has 20 bits */
#elif BITS_PER_MP_LIMB == 64
  G |= x.s.manl >> 50; /* manl has 52 bits */
#else
#error "wrong value of BITS_PER_MP_LIMB"
#endif
  /* now G has 13 bits: 0 <= G < 8192 */
  Gh = G >> 8;
  if (Gh >= 30) /* NaN or Inf have the same encoding in DPD or BID */
    return d;

  if (Gh >= 24)
    {
      h = 8 | (G & 1); /* 8 + G[w+4] */
      exp = (G >> 1) & 1023;
    }
  else /* 0 <= Gh < 24 */
    {
      h = G & 7;
      exp = G >> 3;
    }
  /* now significand is h*2^50 plus the remaining 50 bits */
#if BITS_PER_MP_LIMB == 32
  h = h << 18 | (x.s.manh & 262143);
  l = x.s.manl;
  /* 2^32*h+l < 10^16 i.e. h <= 2328306 */
  d5 = (296 * h + l) % 1000; /* 2^32 = 296 mod 1000 */
  sub_ddmmss (h, l, h, l, 0, d5);
  /* now h*2^32 + l is an exact multiple of 1000 */
  /* first divide by 8 */
  l = ((h & 7) << 29) | (l >> 3);
  h = h >> 3;
  /* now divide exactly by 125 */
  l = l * 652835029; /* 1/125 mod 2^32 */
  h = h / 125;
  /* now 2^32*h+l < 10^13 i.e. h <= 2328 */
  d4 = (296 * h + l) % 1000;
  sub_ddmmss (h, l, h, l, 0, d4);
  l = ((h & 7) << 29) | (l >> 3);
  h = h >> 3;
  l = l * 652835029;
  h = h / 125;
  /* now 2^32*h+l < 10^10 i.e. h <= 2 */
  d3 = (296 * h + l) % 1000;
  sub_ddmmss (h, l, h, l, 0, d4);
  l = (h << 29) | (l >> 3);
  l = l * 652835029;
  /* now l < 10^7 */
#else
  l = h << 50 | (x.s.manl & 1125899906842623);
  /* l < 10*2^50 theoretically, but l < 10^16 in practice */
  d5 = l % 1000;
  l = l / 1000; /* l < 10^13 */
  d4 = l % 1000;
  l = l / 1000; /* l < 10^10 */
  d3 = l % 1000;
  l = l / 1000; /* l < 10^7 */
#endif
  d2 = l % 1000;
  l = l / 1000;
  d1 = l % 1000;
  d0 = l / 10;
  
  /* now encode the declets */
}
#endif /* DPD_FORMAT */

/* construct a decimal64 NaN */
static _Decimal64
get_decimal64_nan (void)
{
  union ieee_double_extract x;
  union ieee_double_decimal64 y;

  x.s.exp = 1984; /* G[0]..G[4] = 11111: quiet NaN */
  y.d = x.d;
  return y.d64;
}

/* construct the decimal64 Inf with given sign */
static _Decimal64
get_decimal64_inf (int negative)
{
  union ieee_double_extract x;
  union ieee_double_decimal64 y;

  x.s.sig = (negative) ? 1 : 0;
  x.s.exp = 1920; /* G[0]..G[4] = 11110: Inf */
  y.d = x.d;
  return y.d64;
}

/* construct the decimal64 zero with given sign */
static _Decimal64
get_decimal64_zero (int negative)
{
  union ieee_double_decimal64 y;

  /* zero has the same representation in binary64 and decimal64 */
  y.d = negative ? DBL_NEG_ZERO : 0.0;
  return y.d64;
}

/* construct the decimal64 smallest non-zero with given sign */
static _Decimal64
get_decimal64_min (int negative)
{
  union ieee_double_extract x;

  x.s.sig = (negative) ? 1 : 0;
  x.s.exp = 0;
#if BITS_PER_MP_LIMB == 32
  x.s.manh = 0;
#endif
  x.s.manl = 1;
  return x.d;
}

/* construct the decimal64 largest finite number with given sign */
static _Decimal64
get_decimal64_max (int negative)
{
  union ieee_double_extract x;

  x.s.sig = (negative) ? 1 : 0;
  x.s.exp = 1919;
#if BITS_PER_MP_LIMB == 32
  x.s.manh = 1048575; /* 2^20-1 */
#endif
  x.s.manl = ~0;
  return x.d;
}

/* for BID format, the exponent exp is meant with a significand of the
   form [0.]sss...sss */
static _Decimal64
fill_decimal64 (int negative, mp_exp_t exp, char *s)
{
  union ieee_double_extract x;
  union ieee_double_decimal64 y;
  unsigned int i, l;
  mp_limb_t rp[2];
  mp_size_t rn;
  int case_i;

  x.s.sig = (negative == 0) ? 0 : 1;
  l = strlen (s);
  exp -= l;
  exp -= -398; /* biases exponent */
  case_i = (l < 16) || (strcmp (s, "9007199254740992") < 0);
  for (i = 0; i < l; i++)
    s[i] -= '0';
  rn = mpn_set_str (rp, (unsigned char *) s, l, 10);
  if (rn == 1)
    rp[1] = 0;
  if (case_i)
    {  /* s < 2^53: case i) */
      x.s.exp = exp << 1;
#if BITS_PER_MP_LIMB == 32
      x.s.manl = rp[0];           /* 32 bits */
      x.s.manh = rp[1] & 1048575; /* 20 low bits */
      x.s.exp |= rp[1] >> 20;     /* 1 bit */
#else
      x.s.manl = rp[0] & 4503599627370495; /* 52 bits */
      x.s.exp |= rp[0] >> 52;
#endif
    }
  else /* s >= 2^53: case ii) */
    {
      x.s.exp = 1536 | (exp >> 1);
#if BITS_PER_MP_LIMB == 32
      x.s.manl = rp[0];
      x.s.manh = (rp[1] ^ 2097152) | ((exp & 1) << 19);
#else
      x.s.manl = (rp[0] ^ 9007199254740992) | ((exp & 1) << 51);
#endif
    }
  y.d = x.d;
  return y.d64;
}

_Decimal64
mpfr_get_decimal64 (mpfr_srcptr src, mp_rnd_t rnd_mode)
{
  int negative;
  mp_exp_t e;

  /* the encoding of NaN, Inf, zero is the same under DPD or BID */
  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (src)))
    {
      if (MPFR_IS_NAN (src))
	return get_decimal64_nan ();

      negative = MPFR_IS_NEG (src);

      if (MPFR_IS_INF (src))
	return get_decimal64_inf (negative);

      MPFR_ASSERTD (MPFR_IS_ZERO(src));
      return get_decimal64_zero (negative);
    }

  e = MPFR_GET_EXP (src);
  negative = MPFR_IS_NEG (src);

  /* the smallest decimal64 number is 10^(-398),
     with 2^(-1323) < 10^(-398) < 2^(-1322) */
  if (MPFR_UNLIKELY (e < -1323)) /* src <= 2^(-1324) < 1/2*10^(-398) */
    {
      if (rnd_mode == GMP_RNDZ || rnd_mode == GMP_RNDN
	  || (rnd_mode == GMP_RNDD && negative == 0)
	  || (rnd_mode == GMP_RNDU && negative != 0))
	return get_decimal64_zero (negative);
      else /* return the smallest non-zero number */
	return get_decimal64_min (negative);
    }
  /* the largest decimal64 number is just below 10^(385) < 2^1279 */
  else if (MPFR_UNLIKELY (e > 1279)) /* then src >= 2^1279 */
    {
      if (GMP_RNDZ || (rnd_mode == GMP_RNDU && negative != 0)
	  || (rnd_mode == GMP_RNDD && negative == 0))
	return get_decimal64_max (negative);
      else
	return get_decimal64_inf (negative);
    }
  else
    {
      /* we need to store the sign (1), the mantissa (16), and the terminating
	 character, thus we need at least 18 characters in s */
      char s[18];
      mpfr_get_str (s, &e, 10, 16, src, rnd_mode);
      /* the smallest normal number is 1.000...000E-383,
	 which corresponds to s=[0.]1000...000 and e=-382 */
      if (e < -382)
	{
	  /* the smallest subnormal number is 0.000...001E-383 = 1E-398,
	     which corresponds to s=[0.]1000...000 and e=-397 */
	  if (e < -397)
	    {
	      if (rnd_mode == GMP_RNDZ || rnd_mode == GMP_RNDN
		  || (rnd_mode == GMP_RNDD && negative == 0)
		  || (rnd_mode == GMP_RNDU && negative != 0))
		return get_decimal64_zero (negative);
	      else /* return the smallest non-zero number */
		return get_decimal64_min (negative);
	    }
	  else
	    {
	      mp_exp_t e2;
	      /* if e = -397 then 16 - (-382 - e) = 1 */
	      mpfr_get_str (s + negative, &e2, 10, 16 - (-382 - e), src,
			    rnd_mode);
	      /* Warning: we can have e2 = e + 1 here, when rounding to
		 nearest or away from zero. */
	      return fill_decimal64 (negative, e2, s + negative);
	    }
	}
      /* the largest number is 9.999...999E+384,
	 which corresponds to s=[0.]9999...999 and e=385 */
      else if (e > 385)
	{
	  if (GMP_RNDZ || (rnd_mode == GMP_RNDU && negative != 0)
	      || (rnd_mode == GMP_RNDD && negative == 0))
	    return get_decimal64_max (negative);
	  else
	    return get_decimal64_inf (negative);
	}
      else /* -382 <= e <= 385 */
	return fill_decimal64 (negative, e, s + negative);
    }
}

#endif /* MPFR_WANT_DECIMAL_FLOATS */
