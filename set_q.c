/* mpfr_set_q -- set a floating-point number from a multiple-precision rational

Copyright 2000, 2001, 2002, 2004, 2005 Free Software Foundation, Inc.

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

/*
 * Set f to z, choosing the smallest precision for f
 * so that z = f*(2^BPML)*zs*2^(RetVal)
 */
static int
set_z (mpfr_ptr f, mpz_srcptr z, mp_size_t *zs)
{
  mp_limb_t *p;
  mp_size_t s;
  int c;
  mp_prec_t pf;

  MPFR_ASSERTD (mpz_sgn (z) != 0);

  /* Remove useless ending 0 */
  for (p = PTR (z), s = *zs = ABS (SIZ (z)) ; *p == 0; p++, s--)
    MPFR_ASSERTD (s >= 0);

  /* Get working precision */
  count_leading_zeros (c, p[s-1]);
  pf = s * BITS_PER_MP_LIMB - c;
  if (pf < MPFR_PREC_MIN)
    pf = MPFR_PREC_MIN;
  mpfr_init2 (f, pf);

  /* Copy Mantissa */
  if (MPFR_LIKELY (c))
    mpn_lshift (MPFR_MANT (f), p, s, c);
  else
    MPN_COPY (MPFR_MANT (f), p, s);
  
  MPFR_SET_SIGN (f, mpz_sgn (z));
  MPFR_SET_EXP (f, 0);

  return -c;
}

/* set f to the rational q */
int
mpfr_set_q (mpfr_ptr f, mpq_srcptr q, mp_rnd_t rnd)
{
  mpz_srcptr num, den;
  mpfr_t n, d;
  int inexact;
  int cn, cd;
  long shift;
  mp_size_t sn, sd;
  MPFR_SAVE_EXPO_DECL (expo);

  num = mpq_numref (q);
  den = mpq_denref (q);
  /* NAN and INF for mpq are not really documented, but could be found */
  if (MPFR_UNLIKELY (mpz_sgn (num) == 0))
    {
      if (MPFR_UNLIKELY (mpz_sgn (den) == 0))
	{
	  MPFR_SET_NAN (f);
	  MPFR_RET_NAN;
	}
      else
	{
	  MPFR_SET_ZERO (f);
	  MPFR_SET_POS (f);
	  MPFR_RET (0);
	}
    }
  if (MPFR_UNLIKELY (mpz_sgn (den) == 0))
    {
      MPFR_SET_INF (f);
      MPFR_SET_SIGN (f, mpz_sgn (num));
      MPFR_RET (0);
    }

  MPFR_SAVE_EXPO_MARK (expo);

  cn = set_z (n, num, &sn);
  cd = set_z (d, den, &sd);
  
  sn -= sd;
  if (MPFR_UNLIKELY (sn > MPFR_EMAX_MAX / BITS_PER_MP_LIMB))
    {
      inexact = mpfr_overflow (f, rnd, MPFR_SIGN (f));
      goto end;
    }
  if (MPFR_UNLIKELY (sn < MPFR_EMIN_MIN / BITS_PER_MP_LIMB -1))
    {
      if (rnd == GMP_RNDN)
	rnd = GMP_RNDZ;
      inexact = mpfr_underflow (f, rnd, MPFR_SIGN (f));
      goto end;
    }

  inexact = mpfr_div (f, n, d, rnd);
  shift = BITS_PER_MP_LIMB*sn+cn-cd;
  MPFR_ASSERTD (shift == BITS_PER_MP_LIMB*sn+cn-cd);
  cd = mpfr_mul_2si (f, f, shift, rnd);
  MPFR_SAVE_EXPO_FREE (expo);
  if (MPFR_UNLIKELY (cd != 0))
    inexact = cd;
  else
    inexact = mpfr_check_range (f, inexact, rnd);
 end:
  mpfr_clear (d);
  mpfr_clear (n);
  return inexact;
}


