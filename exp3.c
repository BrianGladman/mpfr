/* mpfr_exp -- exponential of a floating-point number

Copyright 1999, 2001, 2002, 2003, 2004, 2005 Free Software Foundation, Inc.

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

#include <limits.h>

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

static void
mpfr_exp_rational (mpfr_ptr y, mpz_srcptr p, long r, int m)
{
  long n;
  int i, j, k, l;
  mpz_t *P, *S, *ptoj;
  mp_exp_t diff, expo;
  mp_prec_t precy = MPFR_PREC(y), prec_i_have, accu;
  mp_prec_t *mult, *nb_terms;

  MPFR_ASSERTN ((size_t) m < sizeof (long) * CHAR_BIT - 1);
  n = 1L << m;

  /* Allocate TAB */
  P    = (mpz_t*) (*__gmp_allocate_func) (3*(m+1)*sizeof(mpz_t));
  S    = P + (m+1);
  ptoj = P + 2*(m+1);                     /* ptoj[i] = mantissa^(2^i) */
 
  mult     = (mp_prec_t*) (*__gmp_allocate_func) (2*(m+1)*sizeof(mp_prec_t));
  nb_terms = mult + (m+1);
 
  /* Init var */
  for (i = 0; i <= m; i++)
    {
      mpz_init (P[i]);
      mpz_init (S[i]);
      mpz_init (ptoj[i]);
    }

  /* Set initial var */
  mult[0] = 0;
  mpz_set (ptoj[0], p);
  for (i = 1; i < m; i++)
    mpz_mul (ptoj[i], ptoj[i-1], ptoj[i-1]);
  mpz_set_ui (P[0], 1);
  mpz_set_ui (S[0], 1);
  k = 0;
  nb_terms[0] = 1;
  prec_i_have = 0; 

  /* Main Loop */
  for (i = 1; (prec_i_have < precy) && (i < n); i++)
    {
      /* invariant: P[0]*P[1]*...*P[k] equals i! */
      k++;
      nb_terms[k] = 1;
      mpz_set_ui (P[k], i + 1);
      mpz_set (S[k], P[k]);
      j = i + 1;
      l = 0;
      while ((j & 1) == 0)
        {
          mpz_mul (S[k], S[k], ptoj[l]);
          mpz_mul (S[k-1], S[k-1], P[k]);
          mpz_mul_2exp (S[k-1], S[k-1], r * (1 << l));
          mpz_add (S[k-1], S[k-1], S[k]);
          mpz_mul (P[k-1], P[k-1], P[k]);
          nb_terms[k-1] += nb_terms[k];
	  MPFR_MPZ_SIZEINBASE2 (prec_i_have, P[k]);
	  mult[k] = mult[k-1] + (1 << l) * (r >> 2) + prec_i_have - 1;
          prec_i_have = mult[k];
          /* since mult[k] >= mult[k-1] + nbits(P[k]),
             we have P[0]*...*P[k] <= 2^mult[k] = 2^prec_i_have */
          l++;
          j >>= 1;
          k--;
        }
    }

  /* accumulate all products in P[0] */
  l = 0;
  accu = 0;
  while (k > 0)
    {
      mpz_mul (S[k], S[k], ptoj[MPFR_INT_CEIL_LOG2 (nb_terms[k])]);
      mpz_mul (S[k-1], S[k-1], P[k]);
      accu += nb_terms[k];
      mpz_mul_2exp (S[k-1], S[k-1], r * accu);
      mpz_add (S[k-1], S[k-1], S[k]);
      mpz_mul (P[k-1], P[k-1], P[k]);     
      l++;
      k--;
    }

  /* P[0] now equals i! */
  MPFR_MPZ_SIZEINBASE2 (prec_i_have, S[0]);
  diff = (mp_exp_t) prec_i_have - 2 * (mp_exp_t) precy;
  expo = diff;
  if (diff >= 0)
    mpz_div_2exp (S[0], S[0], diff);
  else 
    mpz_mul_2exp (S[0], S[0], -diff);

  MPFR_MPZ_SIZEINBASE2 (prec_i_have, P[0]);
  diff = (mp_exp_t) prec_i_have - (mp_prec_t) precy;
  expo -= diff;
  if (diff > 0)
    mpz_div_2exp (P[0], P[0], diff);
  else
    mpz_mul_2exp (P[0], P[0], -diff);

  mpz_tdiv_q (S[0], S[0], P[0]);
  mpfr_set_z (y, S[0], GMP_RNDD);
  MPFR_SET_EXP (y, MPFR_GET_EXP (y) + expo - r * (i - 1) );

  for (i = 0; i <= m; i++)
    {
      mpz_clear (P[i]);
      mpz_clear (S[i]);
      mpz_clear (ptoj[i]);
    }
  (*__gmp_free_func) (P, 3*(m+1)*sizeof(mpz_t));
  (*__gmp_free_func) (mult, 2*(m+1)*sizeof(mp_prec_t));
}

#define shift (BITS_PER_MP_LIMB/2)

int
mpfr_exp_3 (mpfr_ptr y, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  mpfr_t t, x_copy, tmp;
  int i, k, loop;
  mpz_t uk;
  mp_exp_t ttt, shift_x;
  unsigned long twopoweri;
  int prec_x;
  mp_prec_t realprec, Prec;
  int iter;
  int inexact = 0;
  MPFR_ZIV_DECL (ziv_loop);

  /* decompose x */
  /* we first write x = 1.xxxxxxxxxxxxx
     ----- k bits -- */
  prec_x = MPFR_INT_CEIL_LOG2 (MPFR_PREC (x)) - MPFR_LOG2_BITS_PER_MP_LIMB;
  if (prec_x < 0)
    prec_x = 0;

  ttt = MPFR_GET_EXP (x);
  mpfr_init2 (x_copy, MPFR_PREC(x));
  mpfr_set (x_copy, x, GMP_RNDD);

  /* we shift to get a number less than 1 */
  if (ttt > 0) 
    {
      shift_x = ttt;
      mpfr_div_2ui (x_copy, x, ttt, GMP_RNDN);
      ttt = MPFR_GET_EXP (x_copy);
    }
  else
    shift_x = 0;
  MPFR_ASSERTD (ttt <= 0);

  /* Init prec and vars */
  realprec = MPFR_PREC (y) + MPFR_INT_CEIL_LOG2 (prec_x + MPFR_PREC (y));
  Prec = realprec + shift + 2 + shift_x;
  mpfr_init2 (t, Prec);
  mpfr_init2 (tmp, Prec);
  mpz_init (uk);

  /* Main loop */
  MPFR_ZIV_INIT (ziv_loop, realprec);
  for (;;)
    {
      k = MPFR_INT_CEIL_LOG2 (Prec) - MPFR_LOG2_BITS_PER_MP_LIMB;

      /* now we have to extract */
      twopoweri = BITS_PER_MP_LIMB;

      /* Particular case for i==0 */
      mpfr_extract (uk, x_copy, 0);
      mpfr_exp_rational (tmp, uk, shift + twopoweri - ttt, k + 1);
      for (loop = 0 ; loop < shift; loop++)
	mpfr_mul (tmp, tmp, tmp, GMP_RNDD);
      twopoweri *=2;

      /* General case */
      iter = (k <= prec_x) ? k : prec_x;
      for (i = 1 ; i <= iter; i++)
        {
          mpfr_extract (uk, x_copy, i);
	  mpfr_exp_rational (t, uk, twopoweri - ttt, k  - i + 1);
	  mpfr_mul (tmp, tmp, t, GMP_RNDD); 
          MPFR_ASSERTN (twopoweri <= LONG_MAX/2);
          twopoweri *=2;
        }

      mpfr_clear_overflow ();
      mpfr_clear_underflow ();
      for (loop = 0 ; loop < shift_x; loop++)
	mpfr_mul (tmp, tmp, tmp, GMP_RNDD);

      if (MPFR_UNLIKELY (mpfr_overflow_p ()))
	{
	  /* We hack to set a FP number outside the valid range so that
	     mpfr_check_range properly generates an overflow */
	  mpfr_setmax (y, __gmpfr_emax);
	  MPFR_EXP (y) ++; 
	  inexact = 1;
	  break;
	}
      else if (MPFR_UNLIKELY (mpfr_underflow_p ()))
	{
	  /* We hack to set a FP number outside the valid range so that
             mpfr_check_range properly generates an underflow.
             Note that the range has been increased to allow a safe
             detection of underflow (MPFR_EMIN_MIN-3 in exp.c) even for
             RNDN */
          mpfr_setmax (y, MPFR_EMIN_MIN-2);
          inexact = -1;
          break;
	}
      else if (mpfr_can_round (tmp, realprec, GMP_RNDD, GMP_RNDZ,
                          MPFR_PREC(y) + (rnd_mode == GMP_RNDN)))
	{
	  inexact = mpfr_set (y, tmp, rnd_mode);	  
	  break;
	}
      MPFR_ZIV_NEXT (ziv_loop, realprec);
      Prec = realprec + shift + 2 + shift_x;
      mpfr_set_prec (t, Prec);
      mpfr_set_prec (tmp, Prec);
    }
  MPFR_ZIV_FREE (ziv_loop);
  
  mpz_clear (uk);
  mpfr_clear (tmp);
  mpfr_clear (t);
  mpfr_clear (x_copy);
  return inexact;
} 
