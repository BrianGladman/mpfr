/* mpfr_gamma -- gamma function

Copyright 2001, 2002, 2003, 2004, 2005 Free Software Foundation.

This file is part of the MPFR Library, and was contributed by Mathieu Dutour.

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

/* The error analysis of gamma has been lost.
   As a consequence, we can't change the algorithm...
   But we may compute the exp(A-k) in the inner loop a lot faster
   but less accurate, so it changes the precsion.
   All MPFR tests still pass anyway */
/* #define USE_PRECOMPUTED_EXP */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* We use the reflection formula
  Gamma(1+t) Gamma(1-t) = - Pi t / sin(Pi (1 + t))
  in order to treat the case x <= 1,
  i.e. with x = 1-t, then Gamma(x) = -Pi*(1-x)/sin(Pi*(2-x))/GAMMA(2-x)
*/

int
mpfr_gamma (mpfr_ptr gamma, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  mpfr_t xp, GammaTrial, tmp, tmp2;
  mpz_t fact;
#ifdef USE_PRECOMPUTED_EXP
  mpfr_t *tab;
#endif
  mp_prec_t Prec, realprec, A, k;
  int compared, inex, is_integer;
  MPFR_GROUP_DECL (group);
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC (("x[%#R]=%R rnd=%d", x, x, rnd_mode),
                 ("gamma[%#R]=%R inexact=%d", gamma, gamma, inex));

  /* Trivial cases */
  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    {
      if (MPFR_IS_NAN (x))
        {
          MPFR_SET_NAN (gamma);
          MPFR_RET_NAN;
        }
      else if (MPFR_IS_INF (x))
        {
          if (MPFR_IS_NEG (x))
            {
              MPFR_SET_NAN (gamma);
              MPFR_RET_NAN;
            }
          else
            {
              MPFR_SET_INF (gamma);
              MPFR_SET_POS (gamma);
              MPFR_RET (0);  /* exact */
            }
        }
      else /* x is zero */
        {
          MPFR_ASSERTD(MPFR_IS_ZERO(x));
          MPFR_SET_INF(gamma);
          MPFR_SET_SAME_SIGN(gamma, x);
          MPFR_RET (0);  /* exact */
        }
    }

  is_integer = mpfr_integer_p (x);
  /* gamma(x) for x a negative integer gives NaN */
  if (is_integer && MPFR_IS_NEG(x))
    {
      MPFR_SET_NAN (gamma);
      MPFR_RET_NAN;
    }

  compared = mpfr_cmp_ui (x, 1);
  if (compared == 0)
    return mpfr_set_ui (gamma, 1, rnd_mode);

  /* if x is an integer that fits into an unsigned long, use mpfr_fac_ui */
  if (is_integer && mpfr_fits_ulong_p (x, GMP_RNDN))
    {
      unsigned long int u;
      u = mpfr_get_ui (x, GMP_RNDN);
      return mpfr_fac_ui (gamma, u - 1, rnd_mode);
    }

  /* check for overflow: according to (6.1.37) in Abramowitz & Stegun,
     gamma(x) >= exp(-x) * x^(x-1/2) * sqrt(2*Pi)
              >= 2 * (x/e)^x / x for x >= 1 */
  if (compared > 0)
    {
      int overflow;

      mpfr_init2 (xp, 53);
      mpfr_set_d (xp, EXPM1, GMP_RNDZ); /* 1/e rounded down */
      mpfr_mul (xp, x, xp, GMP_RNDZ);
      mpfr_pow (xp, xp, x, GMP_RNDZ);
      mpfr_mul_2exp (xp, xp, 1, GMP_RNDZ);
      overflow = MPFR_EXP(xp) > __gmpfr_emax;
      mpfr_clear (xp);
      if (overflow)
        return mpfr_overflow (gamma, rnd_mode, 1);
    }

  MPFR_SAVE_EXPO_MARK (expo);

  /* check for underflow: for x < 1,
     gamma(x) = -Pi*(1-x)/sin(Pi*(2-x))/gamma(2-x).
     Since gamma(2-x) >= 2 * ((2-x)/e)^(2-x) / (2-x), we have
     |gamma(x)| <= Pi*(1-x)*(2-x)/2/((2-x)/e)^(2-x) / |sin(Pi*(2-x))|
                <= 12 * ((2-x)/e)^x / |sin(Pi*(2-x))| */
  if (MPFR_IS_NEG(x))
    {
      int underflow = 0, sgn;

      mpfr_init2 (xp, 53);
      mpfr_init2 (tmp, 53);
      mpfr_init2 (tmp2, 53);
      /* we want an upper bound for 12 * ((2-x)/e)^x;
         since x < 0, y -> y^x is decreasing, thus we need
         a lower bound on (2-x)/e */
      mpfr_ui_sub (xp, 2, x, GMP_RNDZ);
      mpfr_set_d (tmp, EXPM1, GMP_RNDZ); /* 1/e rounded down */
      mpfr_mul (xp, xp, tmp, GMP_RNDZ);
      mpfr_pow (xp, xp, x, GMP_RNDU);
      mpfr_mul_ui (xp, xp, 12, GMP_RNDU);

      /* we need an upper bound on 1/|sin(Pi*(2-x))|,
         thus a lower bound on |sin(Pi*(2-x))|.
         If 2-x is exact, then the error of Pi*(2-x) is (1+u)^2 with u = 2^(-p)
         thus the error on sin(Pi*(2-x)) is less than 1/2ulp + 3Pi(2-x)u,
         assuming u <= 1, thus <= u + 3Pi(2-x)u */
      while (mpfr_ui_sub (tmp, 2, x, GMP_RNDN) != 0)
        {
          mpfr_set_prec (tmp, mpfr_get_prec (tmp) * 3 / 2);
          mpfr_set_prec (tmp2, mpfr_get_prec (tmp));
        }
      mpfr_const_pi (tmp2, GMP_RNDN);
      mpfr_mul (tmp2, tmp2, tmp, GMP_RNDN); /* Pi*(2-x) */
      mpfr_sin (tmp, tmp2, GMP_RNDN); /* sin(Pi*(2-x)) */
      mpfr_mul_ui (tmp2, tmp2, 3, GMP_RNDU); /* 3Pi(2-x) */
      mpfr_add_ui (tmp2, tmp2, 1, GMP_RNDU); /* 3Pi(2-x)+1 */
      mpfr_div_2exp (tmp2, tmp2, mpfr_get_prec (tmp), GMP_RNDU);
      /* if tmp2<|tmp|, we get a lower bound */
      sgn = mpfr_sgn (tmp);
      mpfr_abs (tmp, tmp, GMP_RNDN);
      if (mpfr_cmp (tmp2, tmp) < 0)
        {
          mpfr_sub (tmp, tmp, tmp2, GMP_RNDZ);
          mpfr_div (xp, xp, tmp, GMP_RNDU);
          underflow = MPFR_EXP(xp) <= expo.saved_emin - 2;
        }

      mpfr_clear (xp);
      mpfr_clear (tmp);
      mpfr_clear (tmp2);
      if (underflow) /* the sign is the opposite of that of sin(Pi*(2-x)) */
        {
          MPFR_SAVE_EXPO_FREE (expo);
          return mpfr_underflow (gamma, (rnd_mode == GMP_RNDN) ? GMP_RNDZ : rnd_mode, -sgn);
        }
    }

  realprec = MPFR_PREC (gamma);
  realprec = realprec + MPFR_INT_CEIL_LOG2 (realprec) + 10;

  MPFR_GROUP_INIT_4 (group, realprec + BITS_PER_MP_LIMB,
                     xp, tmp, tmp2, GammaTrial);
  mpz_init (fact);
  MPFR_ZIV_INIT (loop, realprec);
  for (;;)
    {
      /* If compared < 0, we use the reflection formula */
      /* Precision stuff */
      Prec = realprec + 2 * (compared < 0);
      /* A   = (prec_nec-0.5)*CST
         CST = ln(2)/(ln(2*pi))) = 0.38
         This strange formula is just to avoid any overflow */
      A = (Prec/100)*38 + ((Prec%100)*38+100-38/2)/100 - 1;
      /* Estimated_cancel is the amount of bit that will be flushed:
         estimated_cancel = A + ecCST * A;
         ecCST = {1+sup_{x\in [0,1]} x*ln((1-x)/x)}/ln(2) = 1.84
         This strange formula is just to avoid any overflow */
      Prec += 16 + (A + (A + (A/100)*84 + ((A%100)*84)/100));

      /* FIXME: for x near 0, we want 1-x to be exact since the error
         term does not seem to take into account the possible cancellation
         here. Warning: if x < 0, we need one more bit. */
      if (MPFR_EXP(x) < 0 && Prec <= MPFR_PREC(x) - MPFR_EXP(x))
        Prec = MPFR_PREC(x) - MPFR_EXP(x) + 1;

      MPFR_GROUP_REPREC_4 (group, Prec, xp, tmp, tmp2, GammaTrial);

      if (compared < 0)
        mpfr_sub (xp, __gmpfr_one, x, GMP_RNDN);
      else
        mpfr_sub (xp, x, __gmpfr_one, GMP_RNDN);
      mpfr_set_ui (GammaTrial, 0, GMP_RNDN);
      mpz_set_ui (fact, 1);

      /* It is faster to compute exp from k=1 to A, but
         it changes the order of the sum, which change the error
         analysis... Don't change it too much until we recompute
         the error analysis.
         Another trick is to compute factorial k in a MPFR rather than a MPZ
         (Once again it changes the error analysis) */
#ifdef USE_PRECOMPUTED_EXP
      tab = (*__gmp_allocate_func) (sizeof (mpfr_t)*(A-1));
      mpfr_init2 (tab[0], Prec);
      mpfr_exp (tab[0], __gmpfr_one, GMP_RNDN);
      for (k = 1; k < A-1; k++)
        {
          mpfr_init2 (tab[k], Prec);
          mpfr_mul (tab[k], tab[k-1], tab[0], GMP_RNDN);
        }
#endif

      for (k = 1; k < A; k++)
        {
#ifdef USE_PRECOMPUTED_EXP
          mpfr_set (tmp2, tab[A-k-1], GMP_RNDN);
#else
          mpfr_set_ui (tmp, A - k, GMP_RNDN);
          mpfr_exp (tmp2, tmp, GMP_RNDN);
#endif
          /* tmp2 = exp(A-k) */
          mpfr_ui_pow_ui (tmp, A - k, k - 1, GMP_RNDN);
          /* tmp = (A-k)^(k-1) */
          mpfr_mul (tmp2, tmp2, tmp, GMP_RNDN);
          /* tmp2 = (A-k)^(k-1) * exp(A-k) */
          mpfr_sqrt_ui (tmp, A - k, GMP_RNDN);
          mpfr_mul (tmp2, tmp2, tmp, GMP_RNDN);
          /* tmp2 = sqrt(A-k) * (A-k)^(k-1) * exp(A-k) */
          if (k >= 3)
            {
              /* mpfr_fac_ui (tmp, k-1, GMP_RNDN); */
              mpz_mul_ui (fact, fact, k - 1);
              /* fact = (k-1)! */
              mpfr_set_z (tmp, fact, GMP_RNDN);
              mpfr_div (tmp2, tmp2, tmp, GMP_RNDN);
            }
          /* tmp2 = sqrt(A-k) * (A-k)^(k-1) * exp(A-k) / (k-1)! */
          mpfr_add_ui (tmp, xp, k, GMP_RNDN);
          mpfr_div (tmp2, tmp2, tmp, GMP_RNDN);
          /* tmp2 = sqrt(A-k) * (A-k)^(k-1) * exp(A-k) / (k-1)! / (x+k) */
          if ((k & 1) == 0)
            MPFR_CHANGE_SIGN (tmp2);
          /* (-1)^(k-1)*sqrt(A-k) * (A-k)^(k-1) * exp(A-k) / (k-1)! / (x+k) */
          mpfr_add (GammaTrial, GammaTrial, tmp2, GMP_RNDN);
        }
#ifdef DEBUG
      printf("GammaTrial =");
      mpfr_out_str (stdout, 10, 0, GammaTrial, GMP_RNDD);
      printf ("\n");
#endif
      mpfr_const_pi (tmp, GMP_RNDN);
      mpfr_mul_2ui (tmp, tmp, 1, GMP_RNDN); /* 2*Pi */
      mpfr_sqrt (tmp, tmp, GMP_RNDN); /* sqrt(2*Pi) */
      mpfr_add (GammaTrial, GammaTrial, tmp, GMP_RNDN);
      /* sqrt(2*Pi) + sum((-1)^(k-1)*sqrt(A-k)*(A-k)^(k-1)*exp(A-k)
                          /(k-1)!/(xp+k),k=1..A-1) */

      mpfr_add_ui (tmp2, xp, A, GMP_RNDN); /* xp+A */
      mpfr_set_ui_2exp (tmp, 1, -1, GMP_RNDN); /* tmp= 1/2 */
      mpfr_add (tmp, tmp, xp, GMP_RNDN); /* xp+1/2 */
      mpfr_pow (tmp, tmp2, tmp, GMP_RNDN); /* (xp+A)^(xp+1/2) */
      mpfr_mul (GammaTrial, GammaTrial, tmp, GMP_RNDN);
      /* (xp+A)^(xp+1/2) * [sqrt(2*Pi) +
           sum((-1)^(k-1)*(A-k)^(k-1/2)*exp(A-k)/(k-1)!/(xp+k),k=1..A-1)] */

      mpfr_neg (tmp, tmp2, GMP_RNDN); /* -(xp+A) */
      mpfr_exp (tmp, tmp, GMP_RNDN); /* exp(-xp-A) */
      mpfr_mul (GammaTrial, GammaTrial, tmp, GMP_RNDN);

      if (compared < 0)
        {
          mpfr_const_pi (tmp2, GMP_RNDN);
          mpfr_sub (tmp, x, __gmpfr_one, GMP_RNDN);
          mpfr_mul (tmp, tmp2, tmp, GMP_RNDN);
          mpfr_div (GammaTrial, tmp, GammaTrial, GMP_RNDN);
          mpfr_sin (tmp, tmp, GMP_RNDN);
          mpfr_div (GammaTrial, GammaTrial, tmp, GMP_RNDN);
        }
#ifdef DEBUG
      printf("GammaTrial =");
      mpfr_out_str (stdout, 10, 0, GammaTrial, GMP_RNDD);
      printf ("\n");
#endif
#ifdef USE_PRECOMPUTED_EXP
      for (k = 0 ; k < A-1 ; k++)
        mpfr_clear (tab[k]);
      (*__gmp_free_func) (tab, sizeof (mpfr_t) * (A-1) );
#endif
      if (mpfr_can_round (GammaTrial, realprec, GMP_RNDD, GMP_RNDZ,
                          MPFR_PREC(gamma) + (rnd_mode == GMP_RNDN)))
        break;
      MPFR_ZIV_NEXT (loop, realprec);
    }
  MPFR_ZIV_FREE (loop);

  inex = mpfr_set (gamma, GammaTrial, rnd_mode);
  MPFR_GROUP_CLEAR (group);
  mpz_clear (fact);

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (gamma, inex, rnd_mode);
}
