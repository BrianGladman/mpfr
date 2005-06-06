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

#ifdef DEBUG
#include <stdio.h>
#include <stdlib.h>
#endif

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
  i.e. if x = 1-t, then Gamma(x) = -Pi*(1-x)/sin(Pi*(2-x))/GAMMA(2-x)
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

  MPFR_SAVE_EXPO_MARK (expo);

  realprec = MPFR_PREC (gamma);
  realprec = realprec + MPFR_INT_CEIL_LOG2 (realprec) + 10;

  MPFR_GROUP_INIT_4 (group, realprec + BITS_PER_MP_LIMB,
                     xp, tmp, tmp2, GammaTrial);
  mpz_init (fact);
  MPFR_ZIV_INIT (loop, realprec);
  for (;;)
    {
      /* If compared < 0, we use the reflexion formula */
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
         Another trick is to compute factoriel k in a MPFR rather than a MPZ
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
#ifndef USE_PRECOMPUTED_EXP
          mpfr_set_ui (tmp, A - k, GMP_RNDN);
          mpfr_exp (tmp2, tmp, GMP_RNDN);
#else
          mpfr_set (tmp2, tab[A-k-1], GMP_RNDN);
#endif
          mpfr_ui_pow_ui (tmp, A - k, k - 1, GMP_RNDN);
          mpfr_mul (tmp2, tmp2, tmp, GMP_RNDN);
          mpfr_sqrt_ui (tmp, A - k, GMP_RNDN);
          mpfr_mul (tmp2, tmp2, tmp, GMP_RNDN);
          if (k >= 3)
            {
              /* mpfr_fac_ui (tmp, k-1, GMP_RNDN); */
              mpz_mul_ui (fact, fact, k-1);
              mpfr_set_z (tmp, fact, GMP_RNDN);
              mpfr_div (tmp2, tmp2, tmp, GMP_RNDN);
            }
          mpfr_add_ui (tmp, xp, k, GMP_RNDN);
          mpfr_div (tmp2, tmp2, tmp, GMP_RNDN);
          if ((k & 1) == 0)
            MPFR_CHANGE_SIGN (tmp2);
          mpfr_add (GammaTrial, GammaTrial, tmp2, GMP_RNDN);
        }
#ifdef DEBUG
      printf("GammaTrial =");
      mpfr_out_str (stdout, 10, 0, GammaTrial, GMP_RNDD);
      printf ("\n");
#endif
      mpfr_const_pi (tmp, GMP_RNDN);
      mpfr_mul_2ui (tmp, tmp, 1, GMP_RNDN);
      mpfr_sqrt (tmp, tmp, GMP_RNDN);
      mpfr_add (GammaTrial, GammaTrial, tmp, GMP_RNDN);

      mpfr_add_ui (tmp2, xp, A, GMP_RNDN);
      mpfr_set_ui_2exp (tmp, 1, -1, GMP_RNDN); /* tmp= 1/2 */
      mpfr_add (tmp, tmp, xp, GMP_RNDN);
      mpfr_pow (tmp, tmp2, tmp, GMP_RNDN);
      mpfr_mul (GammaTrial, GammaTrial, tmp, GMP_RNDN);

      mpfr_neg (tmp, tmp2, GMP_RNDN);
      mpfr_exp (tmp, tmp, GMP_RNDN);
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
  return mpfr_check_range(gamma, inex, rnd_mode);
}
