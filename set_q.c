/* mpfr_set_q -- set a floating-point number from a multiple-precision rational

Copyright 2000, 2001, 2002, 2004 Free Software Foundation, Inc.

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

/* set f to the rational q */
int
mpfr_set_q (mpfr_ptr f, mpq_srcptr q, mp_rnd_t rnd)
{
  mpz_srcptr num, den;
  mpfr_t n, d;
  int inexact;
  mp_prec_t pnum, pden;
  MPFR_SAVE_EXPO_DECL (expo);

  num = mpq_numref (q);
  if (MPFR_UNLIKELY (mpz_sgn (num) == 0))
    {
      MPFR_SET_ZERO (f);
      MPFR_SET_POS (f);
      MPFR_RET (0);
    }
  den = mpq_denref (q);
  pnum = mpz_sizeinbase (num, 2);
  pden = mpz_sizeinbase (den, 2);
 
  /* Check for underflow */
  if (MPFR_UNLIKELY ((mp_exp_t) pnum - (mp_exp_t) pden+1 <= __gmpfr_emin - 1))
    {
      if (rnd == GMP_RNDN 
	  && ((mp_exp_t) pnum-pden+1 <= __gmpfr_emin - 2
	      || mpq_cmp_si (q, 1, __gmpfr_emin - 2) <= 0))
	rnd = GMP_RNDZ;
      return mpfr_set_underflow (f, rnd, mpq_sgn (q));
    }

  if (MPFR_UNLIKELY (pnum < MPFR_PREC_MIN))
    pnum = MPFR_PREC_MIN;
  if (MPFR_UNLIKELY (pden < MPFR_PREC_MIN))
    pden = MPFR_PREC_MIN;

  MPFR_SAVE_EXPO_MARK (expo);
  mpfr_init2 (n, pnum);
  inexact = mpfr_set_z (n, num, GMP_RNDZ);
  MPFR_ASSERTN (inexact == 0);
  /* result is exact: overflow can occur but we can't handle it */

  mpfr_init2 (d, pden);
  inexact = mpfr_set_z (d, den, GMP_RNDZ);
  MPFR_ASSERTN (inexact == 0);
  /* result is exact: overflow can occur but we can't handle it */

  inexact = mpfr_div (f, n, d, rnd);
  mpfr_clear (n);
  mpfr_clear (d);
  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (f, inexact, rnd);
}
