/* mpfr_erf -- error function of a floating-point number

Copyright 2001, 2003, 2004 Free Software Foundation, Inc.
Contributed by Ludovic Meunier and Paul Zimmermann.

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


#include "mpfr-impl.h"

/* #define DEBUG */

#define EXP1 2.71828182845904523536 /* exp(1) */

int mpfr_erf_0 _MPFR_PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
#if 0
int mpfr_erf_inf _MPFR_PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
int mpfr_erfc_inf _MPFR_PROTO((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
#endif

int
mpfr_erf (mpfr_ptr y, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  double xf;
  int sign_x;
  mp_rnd_t rnd2;
  double n = (double) MPFR_PREC(y);
  int inex;

  sign_x = MPFR_SIGN (x);
  if (MPFR_UNLIKELY( MPFR_IS_SINGULAR(x) ))
    {
      if (MPFR_IS_NAN(x))
	{
	  MPFR_SET_NAN(y);
	  MPFR_RET_NAN;
	}
      else if (MPFR_IS_INF(x)) /* erf(+inf) = +1, erf(-inf) = -1 */
	return mpfr_set_si (y, MPFR_FROM_SIGN_TO_INT(sign_x), GMP_RNDN);
      else /* erf(+0) = +0, erf(-0) = -0 */
        {
          MPFR_ASSERTD(MPFR_IS_ZERO(x));
          return mpfr_set (y, x, GMP_RNDN); /* should keep the sign of x */
        }
    }

  /* now x is neither NaN, Inf nor 0 */

  xf = mpfr_get_d (x, GMP_RNDN);
  xf = xf * xf; /* xf ~ x^2 */

  rnd2 = MPFR_IS_POS_SIGN(sign_x) ? rnd_mode : MPFR_INVERT_RND(rnd_mode);

  /* use expansion at x=0 when e*x^2 <= n (target precision)
     otherwise use asymptotic expansion */

  if (xf > n * LOG2) /* |erf x| = 1 or 1- */
    {
      if (rnd2 == GMP_RNDN || rnd2 == GMP_RNDU)
        {
          if (MPFR_IS_POS_SIGN(sign_x))
            {
              mpfr_set_ui (y, 1, rnd2);
              inex = 1;
            }
          else
            {
              mpfr_set_si (y, -1, rnd2);
              inex = -1;
            }
        }
      else /* round to zero */
        {
          mpfr_setmax (y, 0); /* warning: setmax keeps the old sign of y */
          MPFR_SET_SAME_SIGN(y, x);
          inex = MPFR_IS_POS_SIGN(sign_x) ? -1 : 1;
        }
    }
  else  /* use Taylor */
    {
      inex = mpfr_erf_0 (y, x, rnd_mode);
    }

  return inex;
}

/* return x*2^e */
static
double mul_2exp (double x, mp_exp_t e)
{
  if (e > 0)
    {
      while (e--)
        x *= 2.0;
    }
  else
    {
      while (e++)
        x /= 2.0;
    }

  return x;
}

/* evaluates erf(x) using the expansion at x=0:

   erf(x) = 2/sqrt(Pi) * sum((-1)^k*x^(2k+1)/k!/(2k+1), k=0..infinity)

   Assumes x is neither NaN nor infinite nor zero.
   Assumes also that e*x^2 <= n (target precision).
 */
int
mpfr_erf_0 (mpfr_ptr res, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  mp_prec_t n, m;
  mp_exp_t nuk, sigmak;
  double xf, tauk;
  mpfr_t y, s, t, u;
  unsigned int k;
  long log2tauk;
  int ok;
  int inex;

  n = MPFR_PREC(res); /* target precision */
  xf = mpfr_get_d (x, GMP_RNDN);

  /* initial working precision */
  m = n + (mp_prec_t) (xf * xf / LOG2) + 8;

  mpfr_init2 (y, 2);
  mpfr_init2 (s, 2);
  mpfr_init2 (t, 2);
  mpfr_init2 (u, 2);

  do
    {

      m += __gmpfr_ceil_log2 ((double) n);

      mpfr_set_prec (y, m);
      mpfr_set_prec (s, m);
      mpfr_set_prec (t, m);
      mpfr_set_prec (u, m);

      mpfr_mul (y, x, x, GMP_RNDU); /* err <= 1 ulp */
      mpfr_set_ui (s, 1, GMP_RNDN); /* exact */
      mpfr_set_ui (t, 1, GMP_RNDN); /* exact */
      tauk = 0.0;

      for (k = 1; ; k++)
        {
          mpfr_mul (t, y, t, GMP_RNDU);
          mpfr_div_ui (t, t, k, GMP_RNDU);
          mpfr_div_ui (u, t, 2 * k + 1, GMP_RNDU);
          sigmak = MPFR_EXP(s);
          if (k % 2)
            mpfr_sub (s, s, u, GMP_RNDN);
          else
            mpfr_add (s, s, u, GMP_RNDN);
          sigmak -= MPFR_EXP(s);
          nuk = MPFR_EXP(u) - MPFR_EXP(s);

          if ((nuk < - (mp_exp_t) m) && ((double) k >= xf * xf))
            break;

          /* tauk <- 1/2 + tauk * 2^sigmak + (1+8k)*2^nuk */
          tauk = 0.5 + mul_2exp (tauk, sigmak)
            + mul_2exp (1.0 + 8.0 * (double) k, nuk);
        }

      mpfr_mul (s, x, s, GMP_RNDU);
      MPFR_EXP(s) ++;

      mpfr_const_pi (t, GMP_RNDZ);
      mpfr_sqrt (t, t, GMP_RNDZ);
      mpfr_div (s, s, t, GMP_RNDN);
      tauk = 4.0 * tauk + 11.0; /* final ulp-error on s */
      log2tauk = __gmpfr_ceil_log2 (tauk);

      ok = mpfr_can_round (s, m - log2tauk, GMP_RNDN, GMP_RNDZ,
                           n + (rnd_mode == GMP_RNDN));
    }
  while (ok == 0);

  inex = mpfr_set (res, s, rnd_mode);

  mpfr_clear (y);
  mpfr_clear (t);
  mpfr_clear (u);
  mpfr_clear (s);

  return inex;
}

#if 0
/* evaluates erfc(x) using the expansion at x=infinity:

   sqrt(Pi)*x*exp(x^2)*erfc(x) = 1 + sum((-1)^k*(1*3*...*(2k-1))/(2x^2)^k,k>=1)

   Assumes x is neither NaN nor infinite nor zero.
   Assumes also that e*x^2 > n (target precision).
   Since n >= 2, we have x >= sqrt(2/e), and since
   f(x) := sqrt(Pi)*x*exp(x^2)*erfc(x) is increasing, we have
   f(x) >= f(sqrt(2/e)) ~ 0.7142767512, thus the final partial sum
   should be > 0.5, and MPFR_EXP(s) should always be >= 0.
 */
int
mpfr_erfc_inf (mpfr_ptr res, mpfr_srcptr x, mp_rnd_t rnd)
{
  mp_prec_t n, m;
  mpfr_t y, s, t;
  unsigned long k;
  double tauk;
  long log2tauk;
  mp_exp_t sigmak, nuk;
  double xf = mpfr_get_d1 (x);

  n = MPFR_PREC(res); /* target precision */

  mpfr_init2 (y, 2);
  mpfr_init2 (s, 2);
  mpfr_init2 (t, 2);

  m = n; /* working precision */

  xf = xf * xf; /* approximation of x^2 */

  do
    {
      m += __gmpfr_ceil_log2 ((double) n);

      /* check that 2 * (EXP(x) - 1) * x^2 > m, which ensures the smallest
         term is less than 2^(-m) */
      if (2.0 * (double) (MPFR_EXP(x) - 1) * xf <= (double) m)
        {
          mpfr_clear (y);
          mpfr_clear (s);
          mpfr_clear (t);

          return mpfr_erf_0 (res, x, rnd);
        }

      mpfr_set_prec (y, m);
      mpfr_set_prec (s, m);
      mpfr_set_prec (t, m);

      mpfr_mul (y, x, x, GMP_RNDD); /* err <= 1 ulp */
      MPFR_EXP(y) ++;               /* exact */
      mpfr_set_ui (s, 1, GMP_RNDN); /* exact */
      mpfr_set_ui (t, 1, GMP_RNDN); /* exact */

      tauk = 0.0;
      for (k = 1; k <= (unsigned long) xf; k++)
        {
          mpfr_mul_ui (t, t, 2 * k - 1, GMP_RNDU);
          mpfr_div (t, t, y, GMP_RNDU);
          sigmak = MPFR_EXP(s);
          if (k % 2)
            mpfr_sub (s, s, t, GMP_RNDN);
          else
            mpfr_add (s, s, t, GMP_RNDN);
          sigmak -= MPFR_EXP(s);
          nuk = MPFR_EXP(t) - MPFR_EXP(s);

          if (nuk < - (mp_exp_t) m)
            break;

          /* tauk <- 1/2 + tauk * 2^sigmak + 2^(2k+2+nuk) */
          tauk = 0.5 + mul_2exp (tauk, sigmak)
            + mul_2exp (1.0, 2 * k + 2 + nuk);
        }

      if (nuk >= - (mp_exp_t) m)
        abort();

      mpfr_add_one_ulp (y, GMP_RNDU); /* x^2 rounded up */
      nuk = MPFR_EXP(y);
      mpfr_exp (t, y, GMP_RNDU);
      mpfr_mul (t, t, x, GMP_RNDU);
      mpfr_const_pi (y, GMP_RNDD);
      mpfr_sqrt (y, y, GMP_RNDD);
      mpfr_mul (t, t, y, GMP_RNDN);
      mpfr_div (s, s, t, GMP_RNDN);

      /* final error bound on s */
      tauk = mul_2exp (3.0, nuk + 5) + 2.0 * tauk + 115.0;
      log2tauk = __gmpfr_ceil_log2 (tauk);
    }
  while (mpfr_can_round (s, m - log2tauk, GMP_RNDN, rnd, n) == 0);

  mpfr_set (res, s, rnd);

  mpfr_clear (y);
  mpfr_clear (s);
  mpfr_clear (t);

  return 1;
}

/* evaluates erf(x) using the expansion at x=infinity for erfc(x) = 1 - erf(x).
   Assumes x is neither NaN nor infinite nor zero.
   Assumes also that e*x^2 > n (target precision).
 */
int
mpfr_erf_inf (mpfr_ptr res, mpfr_srcptr x, mp_rnd_t rnd)
{
  mp_prec_t n, m;
  mpfr_t tmp;
  mp_exp_t sh;

  n = MPFR_PREC(res); /* target precision */
  m = n;

  mpfr_init2 (tmp, 2);

  do
    {
      m += __gmpfr_ceil_log2 ((double) n);
      mpfr_set_prec (tmp, m);
      mpfr_erfc_inf (tmp, x, GMP_RNDN);   /* err <= 1/2 ulp */
      sh = MPFR_EXP(tmp);
      mpfr_ui_sub (tmp, 1, tmp, GMP_RNDN); /* err <= 1/2 + 1/2*2^sh */
      sh -= MPFR_EXP(tmp);
      /* the final error is bounded by 2^max(sh, 0) */
      if (sh < 0)
        sh = 0;
    }
  while (mpfr_can_round (tmp, m - sh, GMP_RNDN, rnd, n) == 0);

  mpfr_set (res, tmp, rnd);

  mpfr_clear (tmp);

  return 1;
}
#endif
