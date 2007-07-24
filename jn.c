/* mpfr_j0, mpfr_j1, mpfr_jn -- Bessel functions of 1st kind, integer order.
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

static int mpfr_jn_asympt (mpfr_ptr, long, mpfr_srcptr, mp_rnd_t);

int
mpfr_j0 (mpfr_ptr res, mpfr_srcptr z, mp_rnd_t r)
{
  return mpfr_jn (res, 0, z, r);
}

int
mpfr_j1 (mpfr_ptr res, mpfr_srcptr z, mp_rnd_t r)
{
  return mpfr_jn (res, 1, z, r);
}

int
mpfr_jn (mpfr_ptr res, long n, mpfr_srcptr z, mp_rnd_t r)
{
  int inex;
  unsigned long absn;
  mp_prec_t prec, err;
  mp_exp_t exps, expT;
  mpfr_t y, s, t;
  unsigned long k, zz;
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

  /* check for tiny input for j0: j0(z) = 1 - z^2/4 + ..., more precisely
     |j0(z) - 1| <= z^2/4 for -1 <= z <= 1. */
  if (n == 0)
    MPFR_FAST_COMPUTE_IF_SMALL_INPUT (res, __gmpfr_one, -2 * MPFR_GET_EXP (z),
                                      2, 0, r, return _inexact);

  /* idem for j1: j1(z) = z/2 - z^3/16 + ..., more precisely
     |j1(z) - z/2| <= |z^3|/16 for -1 <= z <= 1, with the sign of j1(z) - z/2
     being the opposite of that of z. */
  if (n == 1)
    /* we first compute 2j1(z) = z - z^3/8 + ..., then divide by 2 using
       the "extra" argument of MPFR_FAST_COMPUTE_IF_SMALL_INPUT. */
    MPFR_FAST_COMPUTE_IF_SMALL_INPUT (res, z, -2 * MPFR_GET_EXP (z), 3,
                                      0, r, mpfr_div_2ui (res, res, 1, r));

  /* we can use the asymptotic expansion as soon as |z| > p log(2)/2,
     but to get some margin we use it for |z| > p/2 */
  if (mpfr_cmp_ui (z, MPFR_PREC(res) / 2 + 3) > 0 ||
      mpfr_cmp_si (z, - ((long) MPFR_PREC(res) / 2 + 3)) < 0)
    {
      inex = mpfr_jn_asympt (res, n, z, r);
      if (inex != 0)
        return inex;
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

  inex = ((n >= 0) || ((n & 1) == 0)) ? mpfr_set (res, s, r)
                                      : mpfr_neg (res, s, r);

  mpfr_clear (y);
  mpfr_clear (s);
  mpfr_clear (t);

  return inex;
}

/* Implements asymptotic expansion (formula 9.2.5 from Abramowitz & Stegun).
   Assumes z > p log(2)/2, where p is the target precision.
   Return 0 if the expansion does not converge enough (the value 0 as inexact
   flag should not happen for normal input).
*/
static int
mpfr_jn_asympt (mpfr_ptr res, long n, mpfr_srcptr z, mp_rnd_t r)
{
  mpfr_t s, c, P, Q, t, iz, err_t, err_s, err_u;
  mp_prec_t w;
  long k;
  int inex, stop, diverge = 0;
  mp_exp_t err2, err;
  MPFR_ZIV_DECL (loop);

  mpfr_init (c);

  w = MPFR_PREC(res) + MPFR_INT_CEIL_LOG2(MPFR_PREC(res)) + 4;

  MPFR_ZIV_INIT (loop, w);
  for (;;)
    {
      mpfr_set_prec (c, w);
      mpfr_init2 (s, w);
      mpfr_init2 (P, w);
      mpfr_init2 (Q, w);
      mpfr_init2 (t, w);
      mpfr_init2 (iz, w);
      mpfr_init2 (err_t, 31);
      mpfr_init2 (err_s, 31);
      mpfr_init2 (err_u, 31);

      /* Approximate sin(z) and cos(z). In the following, err <= k means that
         the approximate value y and the true value x are related by 
         y = x * (1 + u)^k with |u| <= 2^(-w), following Higham's method. */
      mpfr_sin_cos (s, c, z, GMP_RNDN);
      if (MPFR_IS_NEG(z))
        mpfr_neg (s, s, GMP_RNDN); /* compute jn(|z|), fix sign later */
      /* The absolute error on s/c is bounded by 1/2 ulp(1/2) <= 2^(-w-1). */
      mpfr_add (t, s, c, GMP_RNDN);
      mpfr_sub (c, s, c, GMP_RNDN);
      mpfr_swap (s, t);
      /* now s approximates sin(z)+cos(z), and c approximates sin(z)-cos(z),
         with total absolute error bounded by 2^(1-w). */

      /* precompute 1/(8|z|) */
      mpfr_si_div (iz, MPFR_IS_POS(z) ? 1 : -1, z, GMP_RNDN);   /* err <= 1 */
      mpfr_div_2ui (iz, iz, 3, GMP_RNDN);

      /* compute P and Q */
      mpfr_set_ui (P, 1, GMP_RNDN);
      mpfr_set_ui (Q, 0, GMP_RNDN);
      mpfr_set_ui (t, 1, GMP_RNDN); /* current term */
      mpfr_set_ui (err_t, 0, GMP_RNDN); /* error on t */
      mpfr_set_ui (err_s, 0, GMP_RNDN); /* error on P and Q (sum of errors) */
      for (k = 1, stop = 0; stop < 4; k++)
        {
          /* compute next term: t(k)/t(k-1) = (2n+2k-1)(2n-2k+1)/(8kz) */
          mpfr_mul_si (t, t, 2 * (n + k) - 1, GMP_RNDN); /* err <= err_k + 1 */
          mpfr_mul_si (t, t, 2 * (n - k) + 1, GMP_RNDN); /* err <= err_k + 2 */
          mpfr_div_ui (t, t, k, GMP_RNDN);               /* err <= err_k + 3 */
          mpfr_mul (t, t, iz, GMP_RNDN);                 /* err <= err_k + 5 */
          /* the relative error on t is bounded by (1+u)^(5k)-1, which is
             bounded by 6ku for 6ku <= 0.02: first |5 log(1+u)| <= |5.5u|
             for |u| <= 0.15, then |exp(5.5u)-1| <= 6u for |u| <= 0.02. */
          mpfr_mul_ui (err_t, t, 6 * k, MPFR_IS_POS(t) ? GMP_RNDU : GMP_RNDD);
          mpfr_abs (err_t, err_t, GMP_RNDN); /* exact */
          /* the absolute error on t is bounded by err_t * 2^(-w) */
          mpfr_abs (err_u, t, GMP_RNDU);
          mpfr_mul_2ui (err_u, err_u, w, GMP_RNDU); /* t * 2^w */
          mpfr_add (err_u, err_u, err_t, GMP_RNDU); /* max|t| * 2^w */
          if (stop >= 2)
            {
              /* take into account the neglected terms: t * 2^(-w) */
              mpfr_mul_2ui (err_s, err_s, w, GMP_RNDU);
              if (MPFR_IS_POS(t))
                mpfr_add (err_s, err_s, t, GMP_RNDU);
              else
                mpfr_sub (err_s, err_s, t, GMP_RNDU);
              mpfr_div_2ui (err_s, err_s, w, GMP_RNDU);
              stop ++;
            }
          /* if k is odd, add to Q, otherwise to P */
          else if (k & 1)
            {
              /* if k = 1 mod 4, add, otherwise subtract */
              if ((k & 2) == 0)
                mpfr_add (Q, Q, t, GMP_RNDN);
              else
                mpfr_sub (Q, Q, t, GMP_RNDN);
              /* check if the next term is smaller than ulp(Q): if EXP(err_u)
                 <= EXP(Q), since the current term is bounded by
                 err_u * 2^(-w), it is bounded by ulp(Q) */
              if (MPFR_EXP(err_u) <= MPFR_EXP(Q))
                stop ++;
              else
                stop = 0;
            }
          else
            {
              /* if k = 0 mod 4, add, otherwise subtract */
              if ((k & 2) == 0)
                mpfr_add (P, P, t, GMP_RNDN);
              else
                mpfr_sub (P, P, t, GMP_RNDN);
              /* check if the next term is smaller than ulp(P) */
              if (MPFR_EXP(err_u) <= MPFR_EXP(P))
                stop ++;
              else
                stop = 0;
            }
          mpfr_add (err_s, err_s, err_t, GMP_RNDU);
          /* the sum of the rounding errors on P and Q is bounded by
             err_s * 2^(-w) */

          /* stop when start to diverge */
          if (stop < 2 &&
              ((MPFR_IS_POS(z) && mpfr_cmp_ui (z, (k + 1) / 2) < 0) ||
               (MPFR_IS_NEG(z) && mpfr_cmp_si (z, - ((k + 1) / 2)) > 0)))
            {
              /* if we have to stop the series because it diverges, then
                 increasing the precision will most probably fail, since
                 we will stop to the same point, and thus compute a very
                 similar approximation */
              diverge = 1;
              stop = 2; /* force stop */
            }
        }
      /* the sum of the total errors on P and Q is bounded by err_s * 2^(-w) */

      /* Now combine: the sum of the rounding errors on P and Q is bounded by
         err_s * 2^(-w), and the absolute error on s/c is bounded by 2^(1-w) */
      if ((n & 1) == 0) /* n even: P * (sin + cos) + Q (cos - sin) */
        {
          mpfr_mul (c, c, Q, GMP_RNDN); /* Q * (sin - cos) */
          err = MPFR_EXP(c);
          mpfr_mul (s, s, P, GMP_RNDN); /* P * (sin + cos) */
          if (MPFR_EXP(s) > err)
            err = MPFR_EXP(s);
          mpfr_sub (s, s, c, GMP_RNDN);
        }
      else /* n odd: P * (sin - cos) + Q (cos + sin) */
        {
          mpfr_mul (c, c, P, GMP_RNDN); /* P * (sin - cos) */
          err = MPFR_EXP(c);
          mpfr_mul (s, s, Q, GMP_RNDN); /* Q * (sin + cos) */
          if (MPFR_EXP(s) > err)
            err = MPFR_EXP(s);
          mpfr_add (s, s, c, GMP_RNDN);
        }
      if ((n & 2) != 0)
        mpfr_neg (s, s, GMP_RNDN);
      if (MPFR_EXP(s) > err)
        err = MPFR_EXP(s);
      /* the absolute error on s is bounded by P*err(s/c) + Q*err(s/c)
         + err(P)*(s/c) + err(Q)*(s/c) + 3 * 2^(err - w - 1)
         <= (|P|+|Q|) * 2^(1-w) + err_s * 2^(1-w) + 2^err * 2^(1-w),
         since |c|, |old_s| <= 2. */
      err2 = (MPFR_EXP(P) >= MPFR_EXP(Q)) ? MPFR_EXP(P) + 2 : MPFR_EXP(Q) + 2;
      /* (|P| + |Q|) * 2^(1 - w) <= 2^(err2 - w) */
      err = MPFR_EXP(err_s) >= err ? MPFR_EXP(err_s) + 2 : err + 2;
      /* err_s * 2^(1-w) + 2^old_err * 2^(1-w) <= 2^err * 2^(-w) */
      err2 = (err >= err2) ? err + 1 : err2 + 1;
      /* now the absolute error on s is bounded by 2^(err2 - w) */
      
      /* multiply by sqrt(1/(Pi*z)) */
      mpfr_const_pi (c, GMP_RNDN);     /* Pi, err <= 1 */
      mpfr_mul (c, c, z, GMP_RNDN);    /* err <= 2 */
      mpfr_si_div (c, MPFR_IS_POS(z) ? 1 : -1, c, GMP_RNDN); /* err <= 3 */
      mpfr_sqrt (c, c, GMP_RNDN);      /* err<=5/2, thus the absolute error is
                                          bounded by 3*u*|c| for |u| <= 0.25 */
      mpfr_mul (err_t, c, s, MPFR_SIGN(c)==MPFR_SIGN(s) ? GMP_RNDU : GMP_RNDD);
      mpfr_abs (err_t, err_t, GMP_RNDU);
      mpfr_mul_ui (err_t, err_t, 3, GMP_RNDU);
      /* 3*2^(-w)*|old_c|*|s| [see below] is bounded by err_t * 2^(-w) */
      err2 += MPFR_EXP(c);
      /* |old_c| * 2^(err2 - w) [see below] is bounded by 2^(err2-w) */
      mpfr_mul (c, c, s, GMP_RNDN);    /* the absolute error on c is bounded by
                                          1/2 ulp(c) + 3*2^(-w)*|old_c|*|s|
                                          + |old_c| * 2^(err2 - w) */
      /* compute err_t * 2^(-w) + 1/2 ulp(c) = (err_t + 2^EXP(c)) * 2^(-w) */
      err = (MPFR_EXP(err_t) > MPFR_EXP(c)) ? MPFR_EXP(err_t) + 1 : MPFR_EXP(c) + 1;
      /* err_t * 2^(-w) + 1/2 ulp(c) <= 2^(err - w) */
      /* now err_t * 2^(-w) bounds 1/2 ulp(c) + 3*2^(-w)*|old_c|*|s| */
      err = (err >= err2) ? err + 1 : err2 + 1;
      /* the absolute error on c is bounded by 2^(err - w) */

      mpfr_clear (s);
      mpfr_clear (P);
      mpfr_clear (Q);
      mpfr_clear (t);
      mpfr_clear (iz);
      mpfr_clear (err_t);
      mpfr_clear (err_s);
      mpfr_clear (err_u);

      err -= MPFR_EXP(c);
      if (MPFR_LIKELY (MPFR_CAN_ROUND (c, w - err, MPFR_PREC(res), r)))
        break;
      if (diverge != 0)
        {
          mpfr_set (c, z, r); /* will force inex=0 below, which means the
                               asymptotic expansion failed */
          break;
        }
      MPFR_ZIV_NEXT (loop, w);
    }
  MPFR_ZIV_FREE (loop);

  inex = MPFR_IS_POS(z) ? mpfr_set (res, c, r) : mpfr_neg (res, c, r);
  mpfr_clear (c);

  return inex;
}
