/* mpfr_j0, mpfr_j1, mpfr_jn -- Bessel functions of 1st kind, integer order.
   http://www.opengroup.org/onlinepubs/009695399/functions/j0.html

Copyright 2007, 2008, 2009, 2010, 2011 Free Software Foundation, Inc.
Contributed by the Arenaire and Caramel projects, INRIA.

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

/* Relations: j(-n,z) = (-1)^n j(n,z)
              j(n,-z) = (-1)^n j(n,z)
*/

static int mpfr_jn_asympt (mpfr_ptr, long, mpfr_srcptr, mpfr_rnd_t);

int
mpfr_j0 (mpfr_ptr res, mpfr_srcptr z, mpfr_rnd_t r)
{
  return mpfr_jn (res, 0, z, r);
}

int
mpfr_j1 (mpfr_ptr res, mpfr_srcptr z, mpfr_rnd_t r)
{
  return mpfr_jn (res, 1, z, r);
}

/* Estimate k0 such that z^2/4 = k0 * (k0 + n)
   i.e., (sqrt(n^2+z^2)-n)/2 = n/2 * (sqrt(1+(z/n)^2) - 1).
   Return min(2*k0/log(2), ULONG_MAX).
*/
static unsigned long
mpfr_jn_k0 (long n, mpfr_srcptr z)
{
  mpfr_t t, u;
  unsigned long k0;

  mpfr_init2 (t, 32);
  mpfr_init2 (u, 32);
  mpfr_div_si (t, z, n, MPFR_RNDN);
  mpfr_sqr (t, t, MPFR_RNDN);
  mpfr_add_ui (t, t, 1, MPFR_RNDN);
  mpfr_sqrt (t, t, MPFR_RNDN);
  mpfr_sub_ui (t, t, 1, MPFR_RNDN);
  mpfr_mul_si (t, t, n, MPFR_RNDN);
  /* the following is a 32-bit approximation to nearest of log(2) */
  mpfr_set_str_binary (u, "0.10110001011100100001011111111");
  mpfr_div (t, t, u, MPFR_RNDN);
  if (mpfr_fits_ulong_p (t, MPFR_RNDN))
    k0 = mpfr_get_ui (t, MPFR_RNDN);
  else
    k0 = ULONG_MAX;
  mpfr_clear (t);
  mpfr_clear (u);
  return k0;
}

int
mpfr_jn (mpfr_ptr res, long n, mpfr_srcptr z, mpfr_rnd_t r)
{
  int inex;
  unsigned long absn;
  mpfr_prec_t prec, pbound, err;
  mpfr_exp_t exps, expT, diffexp;
  mpfr_t y, s, t, absz;
  unsigned long k, zz, k0;
  MPFR_GROUP_DECL(g);
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC (("x[%#R]=%R n=%d rnd=%d", z, z, n, r),
                 ("y[%#R]=%R", res, res));

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

  MPFR_SAVE_EXPO_MARK (expo);

  /* check for tiny input for j0: j0(z) = 1 - z^2/4 + ..., more precisely
     |j0(z) - 1| <= z^2/4 for -1 <= z <= 1. */
  if (n == 0)
    MPFR_FAST_COMPUTE_IF_SMALL_INPUT (res, __gmpfr_one, -2 * MPFR_GET_EXP (z),
                                      2, 0, r, inex = _inexact; goto end);

  /* idem for j1: j1(z) = z/2 - z^3/16 + ..., more precisely
     |j1(z) - z/2| <= |z^3|/16 for -1 <= z <= 1, with the sign of j1(z) - z/2
     being the opposite of that of z. */
  /* TODO: add a test to trigger an error when
       inex = _inexact; goto end
     is forgotten in MPFR_FAST_COMPUTE_IF_SMALL_INPUT below. */
  if (n == 1)
    /* we first compute 2j1(z) = z - z^3/8 + ..., then divide by 2 using
       the "extra" argument of MPFR_FAST_COMPUTE_IF_SMALL_INPUT. */
    MPFR_FAST_COMPUTE_IF_SMALL_INPUT (res, z, -2 * MPFR_GET_EXP (z), 3,
                                      0, r, mpfr_div_2ui (res, res, 1, r);
                                      inex = _inexact; goto end);

  /* we can use the asymptotic expansion as soon as |z| > p log(2)/2,
     but to get some margin we use it for |z| > p/2 */
  pbound = MPFR_PREC (res) / 2 + 3;
  MPFR_ASSERTN (pbound <= ULONG_MAX);
  MPFR_ALIAS (absz, z, 1, MPFR_EXP (z));
  if (mpfr_cmp_ui (absz, pbound) > 0)
    {
      inex = mpfr_jn_asympt (res, n, z, r);
      if (inex != 0)
        goto end;
    }

  MPFR_GROUP_INIT_3 (g, 32, y, s, t);

  /* check underflow case: |j(n,z)| <= 1/sqrt(2 Pi n) (ze/2n)^n
     (see algorithms.tex) */
  /* FIXME: the code below doesn't detect all the underflow cases. Either
     this should be done, or the generic code should detect underflows. */
  if (absn > 0)
    {
      /* the following is an upper 32-bit approximation to exp(1)/2 */
      mpfr_set_str_binary (y, "1.0101101111110000101010001011001");
      if (MPFR_SIGN(z) > 0)
        mpfr_mul (y, y, z, MPFR_RNDU);
      else
        {
          mpfr_mul (y, y, z, MPFR_RNDD);
          mpfr_neg (y, y, MPFR_RNDU);
        }
      mpfr_div_ui (y, y, absn, MPFR_RNDU);
      /* now y is an upper approximation to |ze/2n|: y < 2^EXP(y),
         thus |j(n,z)| < 1/2*y^n < 2^(n*EXP(y)-1).
         If n*EXP(y) < emin then we have an underflow.
         Note that if emin = MPFR_EMIN_MIN and j = 1, this inequality
         will never be satisfied.
         Warning: absn is an unsigned long. */
      if ((MPFR_GET_EXP (y) < 0 && absn > - expo.saved_emin)
          || (absn <= - MPFR_EMIN_MIN &&
              MPFR_GET_EXP (y) < expo.saved_emin / (mpfr_exp_t) absn))
        {
          MPFR_GROUP_CLEAR (g);
          MPFR_SAVE_EXPO_FREE (expo);
          return mpfr_underflow (res, (r == MPFR_RNDN) ? MPFR_RNDZ : r,
                         (n % 2) ? ((n > 0) ? MPFR_SIGN(z) : -MPFR_SIGN(z))
                                 : MPFR_SIGN_POS);
        }
    }

  /* the logarithm of the ratio between the largest term in the series
     and the first one is roughly bounded by k0, which we add to the
     working precision to take into account this cancellation */
  k0 = mpfr_jn_k0 (absn, z);
  prec = MPFR_PREC (res) + k0 + 2 * MPFR_INT_CEIL_LOG2 (MPFR_PREC (res)) + 3;

  MPFR_ZIV_INIT (loop, prec);
  for (;;)
    {
      MPFR_GROUP_REPREC_3 (g, prec, y, s, t);
      mpfr_pow_ui (t, z, absn, MPFR_RNDN); /* z^|n| */
      mpfr_mul (y, z, z, MPFR_RNDN);       /* z^2 */
      mpfr_clear_erangeflag ();
      zz = mpfr_get_ui (y, MPFR_RNDU);
      /* FIXME: The error analysis is incorrect in case of range error. */
      MPFR_ASSERTN (! mpfr_erangeflag_p ()); /* since mpfr_clear_erangeflag */
      mpfr_div_2ui (y, y, 2, MPFR_RNDN);   /* z^2/4 */
      mpfr_fac_ui (s, absn, MPFR_RNDN);    /* |n|! */
      mpfr_div (t, t, s, MPFR_RNDN);
      if (absn > 0)
        mpfr_div_2ui (t, t, absn, MPFR_RNDN);
      mpfr_set (s, t, MPFR_RNDN);
      /* note: we assume here that the maximal error bound is proportional to
         2^exps, which is true also in the case where s=0 */
      exps = MPFR_IS_ZERO (s) ? MPFR_EMIN_MIN : MPFR_GET_EXP (s);
      expT = exps;
      for (k = 1; ; k++)
        {
          MPFR_LOG_MSG (("loop on k, k = %lu\n", k));
          mpfr_mul (t, t, y, MPFR_RNDN);
          mpfr_neg (t, t, MPFR_RNDN);
          /* Mathematically: absn <= LONG_MAX + 1 <= (ULONG_MAX + 1) / 2,
             and in practice, k is not very large, so that one should have
             k + absn <= ULONG_MAX. */
          MPFR_ASSERTN (absn <= ULONG_MAX - k);
          if (k + absn <= ULONG_MAX / k)
            mpfr_div_ui (t, t, k * (k + absn), MPFR_RNDN);
          else
            {
              mpfr_div_ui (t, t, k, MPFR_RNDN);
              mpfr_div_ui (t, t, k + absn, MPFR_RNDN);
            }
          /* see above note */
          exps = MPFR_IS_ZERO (s) ? MPFR_EMIN_MIN : MPFR_GET_EXP (t);
          if (exps > expT)
            expT = exps;
          mpfr_add (s, s, t, MPFR_RNDN);
          exps = MPFR_IS_ZERO (s) ? MPFR_EMIN_MIN : MPFR_GET_EXP (s);
          if (exps > expT)
            expT = exps;
          /* Above it has been checked that k + absn <= ULONG_MAX. */
          if (MPFR_GET_EXP (t) + (mpfr_exp_t) prec <= exps &&
              zz / (2 * k) < k + absn)
            break;
        }
      /* the error is bounded by (4k^2+21/2k+7) ulp(s)*2^(expT-exps)
         <= (k+2)^2 ulp(s)*2^(2+expT-exps) */
      diffexp = expT - exps;
      err = 2 * MPFR_INT_CEIL_LOG2(k + 2) + 2;
      /* FIXME: Can an overflow occur in the following sum? */
      MPFR_ASSERTN (diffexp >= 0 && err >= 0 &&
                    diffexp <= MPFR_PREC_MAX - err);
      err += diffexp;
      if (MPFR_LIKELY (MPFR_CAN_ROUND (s, prec - err, MPFR_PREC(res), r)))
        break;
      MPFR_ZIV_NEXT (loop, prec);
    }
  MPFR_ZIV_FREE (loop);

  inex = ((n >= 0) || ((n & 1) == 0)) ? mpfr_set (res, s, r)
                                      : mpfr_neg (res, s, r);

  MPFR_GROUP_CLEAR (g);

 end:
  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (res, inex, r);
}

#define MPFR_JN
#include "jyn_asympt.c"
