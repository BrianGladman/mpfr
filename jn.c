/* mpfr_j0, mpfr_j1, mpfr_jn -- Bessel functions of 1st kind and integer order.
   http://www.opengroup.org/onlinepubs/009695399/functions/j0.html

Copyright 2007 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

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

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* Relations: j(-n,z) = (-1)^n j(n,z) */

int
mpfr_j0 (mpfr_ptr res, mpfr_srcptr z, mp_rnd_t r)
{
  return mpfr_jn_si (res, z, 0, r);
}

int
mpfr_j1 (mpfr_ptr res, mpfr_srcptr z, mp_rnd_t r)
{
  return mpfr_jn_si (res, z, 1, r);
}

int
mpfr_jn_si (mpfr_ptr res, mpfr_srcptr z, long n, mp_rnd_t r)
{
  int inex;
  unsigned long absn;
  mp_prec_t prec, err;
  mp_exp_t exps, expT;
  mpfr_t y, s, t;
  unsigned long k, zz;
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC (("x[%#R]=%R n=%d rnd=%d", z, z, n, r),
                 ("y[%#R]=%R inexact=%d", res, res, inex));

  absn = SAFE_ABS (unsigned long, n);

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (z)))
    {
      if (MPFR_IS_NAN (z))
        {
          MPFR_SET_NAN (res);
          MPFR_RET_NAN;
        }
      /* j(n,z) tends to zero when z goes to +Inf or -Inf, oscillating around
         0. We choose to return +0 in that case. */
      else if (MPFR_IS_INF (z)) /* FIXME: according to j(-n,z) = (-1)^n j(n,z)
                                   we might want to give a sign depending on
                                   z and n */
        return mpfr_set_ui (res, 0, r);
      else /* z=0: j(0,0)=1, j(n odd,+/-0) = +/-0 if n > 0, -/+0 if n < 0,
              j(n even,+/-0) = +0 */
        {
          if (n == 0)
            return mpfr_set_ui (res, 1, r);
          else if (absn & 1) /* n odd */
            return (n > 0) ? mpfr_set (res, z, r) : mpfr_neg (res, z, r);
          else /* n even */
            return mpfr_set_ui (res, 0, r);
        }
    }

  mpfr_init2 (y, 32);

  /* check underflow case: |j(n,z)| <= 1/sqrt(2 Pi n) (ze/2n)^n
     (see algorithms.tex) */
  if (absn > 0)
    {
      /* the following is an upper 32-bit approximation of exp(1)/2 */
      mpfr_set_str_binary (y, "1.0101101111110000101010001011001");
      if (MPFR_SIGN(z) > 0)
	mpfr_mul (y, y, z, GMP_RNDU);
      else
	{
	  mpfr_mul (y, y, z, GMP_RNDD);
	  mpfr_neg (y, y, GMP_RNDU);
	}
      mpfr_div_ui (y, y, absn, GMP_RNDU);
      /* now y is an upper approximation of |ze/2n|: y < 2^EXP(y),
	 thus |j(n,z)| < 1/2*y^n < 2^(n*EXP(y)-1).
	 If n*EXP(y) < __gmpfr_emin then we have an underflow.
	 Warning: absn is an unsigned long. */
      if ((MPFR_EXP(y) < 0 && absn > (unsigned long) (-__gmpfr_emin))
	  || (absn <= (unsigned long) (-MPFR_EMIN_MIN) &&
	      MPFR_EXP(y) < __gmpfr_emin / (mp_exp_t) absn))
	{
	  mpfr_clear (y);
	  return mpfr_underflow (res, (r == GMP_RNDN) ? GMP_RNDZ : r,
			 (n % 2) ? ((n > 0) ? MPFR_SIGN(z) : -MPFR_SIGN(z))
				 : MPFR_SIGN_POS);
	}
    }
  
  mpfr_init (s);
  mpfr_init (t);

  prec = MPFR_PREC (res) + MPFR_INT_CEIL_LOG2 (MPFR_PREC (res)) + 3;

  MPFR_ZIV_INIT (loop, prec);
  for (;;)
    {
      mpfr_set_prec (y, prec);
      mpfr_set_prec (s, prec);
      mpfr_set_prec (t, prec);
      mpfr_pow_ui (t, z, absn, GMP_RNDN); /* z^|n| */
      mpfr_mul (y, z, z, GMP_RNDN);       /* z^2 */
      zz = mpfr_get_ui (y, GMP_RNDU);
      MPFR_ASSERTN (zz < ULONG_MAX);
      mpfr_div_2ui (y, y, 2, GMP_RNDN);   /* z^2/4 */
      mpfr_fac_ui (s, absn, GMP_RNDN);    /* |n|! */
      mpfr_div (t, t, s, GMP_RNDN);
      if (absn > 0)
        mpfr_div_2ui (t, t, absn, GMP_RNDN);
      mpfr_set (s, t, GMP_RNDN);
      exps = MPFR_EXP (s);
      expT = exps;
      for (k = 1; ; k++)
        {
          mpfr_mul (t, t, y, GMP_RNDN);
          mpfr_neg (t, t, GMP_RNDN);
          if (k + absn <= ULONG_MAX / k)
            mpfr_div_ui (t, t, k * (k + absn), GMP_RNDN);
          else
            {
              mpfr_div_ui (t, t, k, GMP_RNDN);
              mpfr_div_ui (t, t, k + absn, GMP_RNDN);
            }
          exps = MPFR_EXP (t);
          if (exps > expT)
            expT = exps;
          mpfr_add (s, s, t, GMP_RNDN);
          exps = MPFR_EXP (s);
          if (exps > expT)
            expT = exps;
          if (MPFR_EXP (t) + (mp_exp_t) prec <= MPFR_EXP (s) &&
              zz / (2 * k) < k + n)
            break;
        }
      /* the error is bounded by (4k^2+21/2k+7) ulp(s)*2^(expT-exps)
         <= (k+2)^2 ulp(s)*2^(2+expT-exps) */
      err = 2 * MPFR_INT_CEIL_LOG2(k + 2) + 2 + expT - MPFR_EXP (s);
      if (MPFR_LIKELY (MPFR_CAN_ROUND (s, prec - err, MPFR_PREC(res), r)))
          break;
      MPFR_ZIV_NEXT (loop, prec);
    }
  MPFR_ZIV_FREE (loop);

  inex = (n >= 0) ? mpfr_set (res, s, r) : mpfr_neg (res, s, r);

  mpfr_clear (y);
  mpfr_clear (s);
  mpfr_clear (t);

  return inex;
}
