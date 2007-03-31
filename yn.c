/* mpfr_y0, mpfr_y1, mpfr_yn_si -- Bessel functions of 2nd kind, integer order.
   http://www.opengroup.org/onlinepubs/009695399/functions/y0.html

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

int
mpfr_y0 (mpfr_ptr res, mpfr_srcptr z, mp_rnd_t r)
{
  return mpfr_yn_si (res, z, 0, r);
}

int
mpfr_y1 (mpfr_ptr res, mpfr_srcptr z, mp_rnd_t r)
{
  return mpfr_yn_si (res, z, 1, r);
}

/* compute in s an approximation of S1 = sum((n-k)!/k!*y^k,k=0..n)
   return e >= 0 the exponent difference between the maximal value of |s|
   during the for loop and the final value of |s|.
*/
static mp_exp_t
mpfr_yn_s1 (mpfr_ptr s, mpfr_srcptr y, unsigned long n)
{
  unsigned long k;
  mpz_t f;
  mp_exp_t e, emax;
  
  mpz_init_set_ui (f, 1);
  /* we compute n!*S1 = sum(a[k]*y^k,k=0..n) where a[k] = n!*(n-k)!/k!,
     a[0] = (n!)^2, a[1] = n!*(n-1)!, ..., a[n-1] = n, a[n] = 1 */
  mpfr_set_ui (s, 1, GMP_RNDN); /* a[n] */
  emax = MPFR_EXP(s);
  for (k = n; k-- > 0;)
    {
      /* a[k]/a[k+1] = (n-k)!/k!/(n-(k+1))!*(k+1)! = (k+1)*(n-k) */
      mpfr_mul (s, s, y, GMP_RNDN);
      mpz_mul_ui (f, f, n - k);
      mpz_mul_ui (f, f, k + 1);
      /* invariant: f = a[k] */
      mpfr_add_z (s, s, f, GMP_RNDN);
      e = MPFR_EXP(s);
      if (e > emax)
	emax = e;
    }
  /* now we have f = (n!)^2 */
  mpz_sqrt (f, f);
  mpfr_div_z (s, s, f, GMP_RNDN);
  mpz_clear (f);
  return emax - MPFR_EXP(s);
}

/* compute in s an approximation of
   S3 = c*sum((h(k)+h(n+k))*y^k/k!/(n+k)!,k=0..infinity)
   where h(k) = 1 + 1/2 + ... + 1/k
   k=0: h(n)
   k=1: 1+h(n+1)
   k=2: 3/2+h(n+2)
   Returns e such that the error is bounded by 2^e ulp(s).
*/
static mp_exp_t
mpfr_yn_s3 (mpfr_ptr s, mpfr_srcptr y, mpfr_srcptr c, unsigned long n)
{
  unsigned long k, zz;
  mpfr_t t, u;
  mpz_t p, q; /* p/q will store h(k)+h(n+k) */
  mp_exp_t exps, expU;

  zz = mpfr_get_ui (y, GMP_RNDU); /* y = z^2/4 */
  MPFR_ASSERTN (zz < ULONG_MAX - 2);
  zz += 2; /* z^2 <= 2^zz */
  mpz_init_set_ui (p, 0);
  mpz_init_set_ui (q, 1);
  /* initialize p/q to h(n) */
  for (k = 1; k <= n; k++)
    {
      /* p/q + 1/k = (k*p+q)/(q*k) */
      mpz_mul_ui (p, p, k);
      mpz_add (p, p, q);
      mpz_mul_ui (q, q, k);
    }
  mpfr_init2 (t, MPFR_PREC(s));
  mpfr_init2 (u, MPFR_PREC(s));
  mpfr_fac_ui (t, n, GMP_RNDN);
  mpfr_div (t, c, t, GMP_RNDN);    /* c/n! */
  mpfr_mul_z (u, t, p, GMP_RNDN);
  mpfr_div_z (s, u, q, GMP_RNDN);
  exps = MPFR_EXP (s);
  expU = exps;
  for (k = 1; ;k ++)
    {
      /* update t */
      mpfr_mul (t, t, y, GMP_RNDN);
      mpfr_div_ui (t, t, k, GMP_RNDN);
      mpfr_div_ui (t, t, n + k, GMP_RNDN);
      /* update p/q:
	 p/q + 1/k + 1/(n+k) = [p*k*(n+k) + q*(n+k) + q*k]/(q*k*(n+k)) */
      mpz_mul_ui (p, p, k);
      mpz_mul_ui (p, p, n + k);
      mpz_addmul_ui (p, q, n + 2 * k);
      mpz_mul_ui (q, q, k);
      mpz_mul_ui (q, q, n + k);
      mpfr_mul_z (u, t, p, GMP_RNDN);
      mpfr_div_z (u, u, q, GMP_RNDN);
      exps = MPFR_EXP (u);
      if (exps > expU)
	expU = exps;
      mpfr_add (s, s, u, GMP_RNDN);
      exps = MPFR_EXP (s);
      if (exps > expU)
	expU = exps;
      if (MPFR_EXP (u) + (mp_exp_t) MPFR_PREC (u) < MPFR_EXP (s) &&
	  zz / (2 * k) < k + n)
	break;
    }
  mpfr_clear (t);
  mpfr_clear (u);
  mpz_clear (p);
  mpz_clear (q);
  exps = expU - MPFR_EXP (s);
  /* the error is bounded by (6k^2+33/2k+11) 2^exps ulps 
     <= 8*(k+2)^2 2^exps ulps */
  return 3 + 2 * MPFR_INT_CEIL_LOG2(k + 2) + exps;
}

int
mpfr_yn_si (mpfr_ptr res, mpfr_srcptr z, long n, mp_rnd_t r)
{
  int inex;
  unsigned long absn;
  mp_prec_t prec;
  mp_exp_t err1, err2, err3;
  mpfr_t y, s1, s2, s3, zn;
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC (("x[%#R]=%R n=%d rnd=%d", z, z, n, r),
                 ("y[%#R]=%R", res, res));

  absn = SAFE_ABS (unsigned long, n);

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (z)))
    {
      if (MPFR_IS_NAN (z))
        {
          MPFR_SET_NAN (res); /* y(n,NaN) = NaN */
          MPFR_RET_NAN;
        }
      /* y(n,z) tends to zero when z goes to +Inf, oscillating around
         0. We choose to return +0 in that case. */
      else if (MPFR_IS_INF (z))
	{
	  if (MPFR_SIGN(z) > 0)
	    return mpfr_set_ui (res, 0, r);
	  else /* y(n,-Inf) = NaN */
	    {
	      MPFR_SET_NAN (res);
	      MPFR_RET_NAN;
	    }
	}
      else /* y(n,z) tends to -Inf for n >= 0 or n even, to +Inf otherwise,
	      when z goes to zero */
        {
	  MPFR_SET_INF(res);
	  if (n >= 0 || (n & 1) == 0)
	    MPFR_SET_NEG(res);
	  else
	    MPFR_SET_POS(res);
	  MPFR_RET(0);
        }
    }

  /* for z < 0, y(n,z) is imaginary except when j(n,|z|) = 0, which we
     assume does not happen for a rational z. */
  if (MPFR_SIGN(z) < 0)
    {
      MPFR_SET_NAN (res);
      MPFR_RET_NAN;
    }

  /* now z is not singular, and z > 0 */

  mpfr_init (y);
  mpfr_init (s1);
  mpfr_init (s2);
  mpfr_init (s3);

  prec = MPFR_PREC(res) + 2 * MPFR_INT_CEIL_LOG2 (MPFR_PREC (res)) + 13;
  MPFR_ZIV_INIT (loop, prec);
  for (;;)
    {
      mpfr_set_prec (y, prec);
      mpfr_set_prec (s1, prec);
      mpfr_set_prec (s2, prec);
      mpfr_set_prec (s3, prec);

      mpfr_mul (y, z, z, GMP_RNDN);
      mpfr_div_2ui (y, y, 2, GMP_RNDN); /* z^2/4 */

      /* store (z/2)^n temporarily in s2 */
      mpfr_pow_ui (s2, z, absn, GMP_RNDN);
      mpfr_div_2si (s2, s2, absn, GMP_RNDN);

      /* compute S1 * (z/2)^(-n) */
      if (n == 0)
	{
	  mpfr_set_ui (s1, 0, GMP_RNDN);
	  err1 = 0;
	}
      else
	err1 = mpfr_yn_s1 (s1, y, absn - 1);
      mpfr_div (s1, s1, s2, GMP_RNDN); /* (z/2)^(-n) * S1 */
      /* See algorithms.tex: the relative error on s1 is bounded by
	 (3n+3)*2^(e+1-prec). */
      err1 = MPFR_INT_CEIL_LOG2 (3 * absn + 3) + err1 + 1;
      /* rel_err(s1) <= 2^(err1-prec), thus err(s1) <= 2^err1 ulps */

      /* compute (z/2)^n * S3 */
      mpfr_neg (y, y, GMP_RNDN); /* -z^2/4 */
      err3 = mpfr_yn_s3 (s3, y, s2, absn); /* (z/2)^n * S3 */
      /* the error on s3 is bounded by 2^err3 ulps */

      /* add s1+s3 */
      err1 += MPFR_EXP(s1);
      mpfr_add (s1, s1, s3, GMP_RNDN);
      /* the error is bounded by 1/2 + 2^err1*2^(- EXP(s1))
	 + 2^err3*2^(EXP(s3) - EXP(s1)) */
      err3 += MPFR_EXP(s3);
      err1 = (err3 > err1) ? err3 + 1 : err1 + 1;
      err1 -= MPFR_EXP(s1);
      err1 = (err1 >= 0) ? err1 + 1 : 1;
      /* now the error on s1 is bounded by 2^err1*ulp(s1) */

      /* compute S2 */
      mpfr_div_2ui (s2, z, 1, GMP_RNDN); /* z/2 */
      mpfr_log (s2, s2, GMP_RNDN); /* log(z/2) */
      mpfr_const_euler (s3, GMP_RNDN);
      err2 = MPFR_EXP(s2) > MPFR_EXP(s3) ? MPFR_EXP(s2) : MPFR_EXP(s3);
      mpfr_add (s2, s2, s3, GMP_RNDN); /* log(z/2) + gamma */
      err2 -= MPFR_EXP(s2);
      mpfr_mul_2ui (s2, s2, 1, GMP_RNDN); /* 2*(log(z/2) + gamma) */
      mpfr_jn_si (s3, z, absn, GMP_RNDN); /* Jn(z) */
      mpfr_mul (s2, s2, s3, GMP_RNDN); /* 2*(log(z/2) + gamma)*Jn(z) */
      err2 += 4; /* the error on s2 is bounded by 2^err2 ulps, see
		    algorithms.tex */

      /* add all three sums */
      err1 += MPFR_EXP(s1); /* the error on s1 is bounded by 2^err1 */
      err2 += MPFR_EXP(s2); /* the error on s2 is bounded by 2^err2 */
      mpfr_sub (s2, s2, s1, GMP_RNDN); /* s2 - (s1+s3) */
      err2 = (err1 > err2) ? err1 + 1 : err2 + 1;
      err2 -= MPFR_EXP(s2);
      err2 = (err2 >= 0) ? err2 + 1 : 1;
      /* now the error on s2 is bounded by 2^err2*ulp(s2) */
      mpfr_const_pi (y, GMP_RNDN); /* error bounded by 1 ulp */
      mpfr_div (s2, s2, y, GMP_RNDN); /* error bounded by 2^(err2+1)*ulp(s2) */
      err2 ++;

      if (MPFR_LIKELY (MPFR_CAN_ROUND (s2, prec - err2, MPFR_PREC(res), r)))
	break;
      MPFR_ZIV_NEXT (loop, prec);
    }
  MPFR_ZIV_FREE (loop);

  inex = (n >= 0 || (n & 1) == 0)
    ? mpfr_set (res, s2, r)
    : mpfr_neg (res, s2, r);

  mpfr_clear (y);
  mpfr_clear (s1);
  mpfr_clear (s2);
  mpfr_clear (s3);

  return inex;
}
