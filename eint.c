/* mpfr_eint, mpfr_eint1 -- the exponential integral

Copyright 2005 Free Software Foundation, Inc.

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
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA. */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* eint1(x) = -gamma - log(x) - sum((-1)^k*z^k/k/k!, k=1..infinity) for x > 0
            = - eint(-x) for x < 0
   where
   eint (x) = gamma + log(x) + sum(z^k/k/k!, k=1..infinity) for x > 0
   eint (x) is undefined for x < 0.
*/

/* compute in y an approximation of sum(x^k/k/k!, k=1..infinity),
   and return e such that the absolute error is bound by 2^e ulp(y) */
static mp_exp_t
mpfr_eint_aux (mpfr_t y, mpfr_srcptr x)
{
  mp_prec_t w = MPFR_PREC(y), w_minus_e;
  unsigned long k;
  mpz_t m, s, t, u;
  mp_exp_t e;
  mpfr_t eps; /* dynamic (absolute) error bound on t */
  mpfr_t erru, errs;

  mpz_init (s); /* initializes to 0 */
  mpz_init (t);
  mpz_init (u);
  mpz_init (m);
  mpfr_init2 (eps, 31);
  mpfr_init2 (erru, 31);
  mpfr_init2 (errs, 31);
  e = mpfr_get_z_exp (m, x); /* x = m * 2^e */
  if (mpz_sizeinbase (m, 2) > w)
    {
      e += mpz_sizeinbase (m, 2) - w;
      mpz_tdiv_q_2exp (m, m, mpz_sizeinbase (m, 2) - w);
    }
  /* remove trailing zeroes from m: this will speed up much cases where
     x is a small integer divided by a power of 2 */
  k = mpz_scan1 (m, 0);
  mpz_tdiv_q_2exp (m, m, k);
  e += k;
  MPFR_ASSERTN(e < 0 || w >= (mp_prec_t) e);
  w_minus_e = (e < 0) ? w + (-e) : w - (mp_prec_t) e;
  /* initialize t to 2^w */
  mpz_set_ui (t, 1);
  mpz_mul_2exp (t, t, w);
  mpfr_set_ui (eps, 0, GMP_RNDN); /* eps[0] = 0 */
  mpfr_set_ui (errs, 0, GMP_RNDN);
  for (k = 1;; k++)
    {
      /* let eps[k] be the absolute error on t[k]:
         since t[k] = trunc(t[k-1]*m*2^e/k), we have
         eps[k+1] <= 1 + eps[k-1]*m*2^e/k + t[k-1]*m*2^(1-w)*2^e/k
                  =  1 + (eps[k-1] + t[k-1]*2^(1-w))*m*2^e/k
                  = 1 + (eps[k-1]*2^(w-1) + t[k-1])*2^(1-w)*m*2^e/k */
      mpfr_mul_2exp (eps, eps, w - 1, GMP_RNDU);
      mpfr_add_z (eps, eps, t, GMP_RNDU);
      mpfr_div_2exp (eps, eps, w - 1, GMP_RNDU);
      mpfr_mul_2exp (eps, eps, mpz_sizeinbase (m, 2), GMP_RNDU);
      if (e < 0)
        mpfr_div_2exp (eps, eps, -e, GMP_RNDU);
      else
        mpfr_mul_2exp (eps, eps, e, GMP_RNDU);
      mpfr_div_ui (eps, eps, k, GMP_RNDU);
      mpfr_add_ui (eps, eps, 1, GMP_RNDU);
      mpz_mul (t, t, m);
      if (e < 0)
        mpz_tdiv_q_2exp (t, t, -e);
      else
        mpz_mul_2exp (t, t, e);
      mpz_tdiv_q_ui (t, t, k);
      mpz_tdiv_q_ui (u, t, k);
      mpz_add (s, s, u);
      /* the absolute error on u is <= 1 + eps[k]/k */
      mpfr_div_ui (erru, eps, k, GMP_RNDU);
      mpfr_add_ui (erru, erru, 1, GMP_RNDU);
      /* and that on s is the sum of all errors on u */
      mpfr_add (errs, errs, erru, GMP_RNDU);
      /* we are done when t is smaller than errs */
      if (mpz_sizeinbase (t, 2) < MPFR_EXP(errs))
        break;
    }
  /* the truncation error is bounded by (|t|+eps)/k*(|x|/k + |x|^2/k^2 + ...)
     <= (|t|+eps)/k*|x|/(k-|x|) */
  mpz_abs (t, t);
  mpfr_add_z (eps, eps, t, GMP_RNDU);
  mpfr_div_ui (eps, eps, k, GMP_RNDU);
  mpfr_abs (erru, x, GMP_RNDU); /* |x| */
  mpfr_mul (eps, eps, erru, GMP_RNDU);
  mpfr_ui_sub (erru, k, erru, GMP_RNDD);
  if (MPFR_SIGN(erru) < 0)
    {
      /* the truncated series does not converge, return fail */
      e = w;
    }
  else
    {
      mpfr_div (eps, eps, erru, GMP_RNDU);
      mpfr_add (errs, errs, eps, GMP_RNDU);
      mpfr_set_z (y, s, GMP_RNDN);
      mpfr_div_2exp (y, y, w, GMP_RNDN);
      /* errs was an absolute error bound on s. We must convert it to an error
         in terms of ulp(y). Since ulp(y) = 2^(EXP(y)-PREC(y)), we must
         divide the error by 2^(EXP(y)-PREC(y)), but since we divided also
         y by 2^w = 2^PREC(y), we must simply divide by 2^EXP(y). */
      e = MPFR_EXP(errs) - MPFR_EXP(y);
    }
  mpfr_clear (eps);
  mpfr_clear (erru);
  mpfr_clear (errs);
  mpz_clear (s);
  mpz_clear (t);
  mpz_clear (u);
  mpz_clear (m);
  return e;
}

int 
mpfr_eint (mpfr_ptr y, mpfr_srcptr x, mp_rnd_t rnd)
{
  int inex;
  mpfr_t tmp, ump;
  mp_exp_t err, te;
  mp_prec_t prec;
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC (("x[%#R]=%R rnd=%d", x, x, rnd),
                 ("y[%#R]=%R inexact=%d", y, y, inex));

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    {
      /* exp(NaN) = exp(-Inf) = NaN */
      if (MPFR_IS_NAN (x) || (MPFR_IS_INF (x) && MPFR_IS_NEG(x)))
        {
          MPFR_SET_NAN (y);
          MPFR_RET_NAN;
        }
      /* eint(+inf) = +inf */
      else if (MPFR_IS_INF (x))
        {
          MPFR_SET_INF(y);
          MPFR_SET_POS(y);
          MPFR_RET(0);
        }
      else /* eint(+/-0) = -Inf */
        {
          MPFR_SET_INF(y);
          MPFR_SET_NEG(y);
          MPFR_RET(0);
        }
    }

  /* eint(x) = NaN for x < 0 */
  if (MPFR_IS_NEG(x))
    {
      MPFR_SET_NAN (y);
      MPFR_RET_NAN;
    }

  /* Init stuff */
  MPFR_SAVE_EXPO_MARK (expo);  
  prec = MPFR_PREC (y) + 2 * MPFR_INT_CEIL_LOG2 (MPFR_PREC (y)) + 6;
  mpfr_init2 (tmp, prec);
  mpfr_init2 (ump, prec);

  /* eint() has a root 0.37250741078136663446..., so if x is near,
     already take more bits */
  if (MPFR_EXP(x) == -1) /* 1/4 <= x < 1/2 */
    {
      double d;
      d = mpfr_get_d (x, GMP_RNDN) - 0.37250741078136663;
      d = (d == 0.0) ? -53 : __gmpfr_ceil_log2 (d);
      prec += -d;
    }

  MPFR_ZIV_INIT (loop, prec);            /* Initialize the ZivLoop controler */
  for (;;)                               /* Infinite loop */
    {
      err = mpfr_eint_aux (tmp, x); /* error <= 2^err ulp(tmp) */
      te = MPFR_EXP(tmp);
      mpfr_const_euler (ump, GMP_RNDN); /* 0.577 -> EXP(ump)=0 */
      mpfr_add (tmp, tmp, ump, GMP_RNDN);
      /* error <= 1/2 + 1/2*2^(EXP(ump)-EXP(tmp)) + 2^(te-EXP(tmp)+err)
               <= 1/2 + 2^(MAX(EXP(ump), te+err+1) - EXP(tmp))
               <= 2^(MAX(0, 1 + MAX(EXP(ump), te+err+1) - EXP(tmp))) */
      err = MAX(1, te + err + 2) - MPFR_EXP(tmp);
      err = MAX(0, err);
      te = MPFR_EXP(tmp);
      mpfr_log (ump, x, GMP_RNDN);
      mpfr_add (tmp, tmp, ump, GMP_RNDN);
      /* same formula as above, except now EXP(ump) is not 0 */
      err = MAX(0, 1 + MAX(MPFR_EXP(ump), te + err + 1) - MPFR_EXP(tmp));
      err = MPFR_PREC(tmp) - err;
      if (MPFR_LIKELY (MPFR_CAN_ROUND (tmp, err, MPFR_PREC (y), rnd)))
	break;
      MPFR_ZIV_NEXT (loop, prec);        /* Increase used precision */
      mpfr_set_prec (tmp, prec);
      mpfr_set_prec (ump, prec);
    }
  MPFR_ZIV_FREE (loop);                  /* Free the ZivLoop Controler */

  inex = mpfr_set (y, tmp, rnd);    /* Set y to the computed value */
  mpfr_clear (tmp);
  mpfr_clear (ump);

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (y, inex, rnd);
}
