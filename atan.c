/* mpfr_atan -- arc-tangent of a floating-point number

Copyright 2001, 2002, 2003, 2004 Free Software Foundation.

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
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

#define CST   2.27  /* CST=1+ln(2.4)/ln(2) */

/*
#define A
#define A1 1
#define A2 2
#define C
#define C1  3
#define C2  2
#define NO_FACTORIAL
#define GENERIC mpfr_atan_aux
#include "generic.c"
#undef C
#undef C1
#undef C2
#undef A
#undef A1
#undef A2
#undef NO_FACTORIAL
#undef GENERIC
*/
static void
mpfr_atan_aux (mpfr_ptr y, mpz_srcptr p, long r, int m)
{
  unsigned long n,i,k,j,l;
  mpz_t *S, *T, *ptoj;
  mp_exp_t diff, expo;
  TMP_DECL (maker);

  TMP_MARK (maker);
  S    = (mpz_t*) TMP_ALLOC (3*(m+1)*sizeof (mpz_t));
  ptoj = S + 1*(m+1);
  T    = S + 2*(m+1);

  for (i = 0 ; i <= m ; i++) {
    mpz_init (S[i]);
    mpz_init (ptoj[i]);
    mpz_init (T[i]);
  }

  mpz_mul_2exp (ptoj[0], p, 1);
  for (i=1;i<m;i++) 
    mpz_mul (ptoj[i], ptoj[i-1], ptoj[i-1]);

  mpz_set_ui (S[0], 1);
  mpz_set_ui (T[0], 1);

  k = 0;
  n = 1UL << m;
  for (i = 1 ; i < n ; i++) {
    k++;
    mpz_set_ui (T[k], 1 + 2*i);
    mpz_set_ui (S[k], 1);

    for (j = i+1, l = 0 ; (j & 1) == 0 ; l++, j>>=1, k--) {
      mpz_mul (S[k], S[k], ptoj[l]);
      mpz_mul (S[k], S[k], T[k-1]);
      mpz_mul (S[k-1], S[k-1], T[k]);
      mpz_mul_2exp (S[k-1], S[k-1], (r+1)*(1<<l));
      mpz_add (S[k-1], S[k-1], S[k]);
      mpz_mul (T[k-1], T[k-1], T[k]);
    }
  }

  diff = mpz_sizeinbase (S[0], 2) - 2*MPFR_PREC (y);
  expo = diff;
  if (diff >=0)
    mpz_fdiv_q_2exp (S[0],S[0],diff);
  else
    mpz_mul_2exp (S[0],S[0],-diff);

  mpz_mul_2exp (T[0], T[0], (1<<m) -1);
  diff = mpz_sizeinbase (T[0], 2) - MPFR_PREC (y);
  expo -= diff;
  if (diff >=0)
    mpz_fdiv_q_2exp (T[0], T[0],diff);
  else
    mpz_mul_2exp (T[0], T[0],-diff);

  mpz_tdiv_q (S[0], S[0], T[0]);
  mpfr_set_z (y, S[0], GMP_RNDD);
  MPFR_SET_EXP (y, MPFR_EXP (y) + expo - r*(i-1) );

  for (i = 0 ; i <= m ; i++) {
    mpz_clear (S[i]);
    mpz_clear (ptoj[i]);
    mpz_clear (T[i]);
  }
  TMP_FREE (maker);
}

int
mpfr_atan (mpfr_ptr atan, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  mpfr_t xp;
  mpfr_t arctgt, sk, tmp, tmp2;
  mpz_t  ukz;
  int comparaison, sign, inexact, inexact2;
  mp_exp_t supplement, estimated_delta;
  mp_prec_t prec, realprec;
  mp_exp_t exptol;
  int i, n0;
  unsigned long twopoweri;
  MPFR_SAVE_EXPO_DECL (expo);

  /* Singular cases */
  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    {
      if (MPFR_IS_NAN (x))
	{
	  MPFR_SET_NAN (atan);
	  MPFR_RET_NAN;
	}
      else if (MPFR_IS_INF (x))
	{
	  if (MPFR_IS_POS (x))  /* arctan(+inf) = Pi/2 */
	    inexact = mpfr_const_pi (atan, rnd_mode);
	  else /* arctan(-inf) = -Pi/2 */
	    {
	      inexact = -mpfr_const_pi (atan, 
					MPFR_INVERT_RND (rnd_mode));
	      MPFR_CHANGE_SIGN (atan);
	    }
	  inexact2 = mpfr_div_2ui (atan, atan, 1, rnd_mode);
	  if (MPFR_UNLIKELY (inexact2)) 
	    inexact = inexact2; /* An underflow occurs */
	  MPFR_RET (inexact);
	}
      else /* x is necessarily 0 */
	{
          MPFR_ASSERTD (MPFR_IS_ZERO (x));
	  MPFR_SET_ZERO (atan);
	  MPFR_RET (0);
	}
    }

  /* Set x_p=|x| */
  sign = MPFR_SIGN (x);
  xp[0] = x[0];                 /* Hack to avoid copying X */
  MPFR_SET_POS (xp);

  /* Other simple case arctang(-+1)=-+pi/4 */
  comparaison = mpfr_cmp_ui (xp, 1);
  if (MPFR_UNLIKELY (comparaison == 0))
    {
      inexact = mpfr_const_pi (atan, MPFR_IS_POS_SIGN (sign) ? rnd_mode
                               : MPFR_INVERT_RND (rnd_mode));
      if (MPFR_IS_NEG_SIGN (sign))
        {
          inexact = -inexact;
          MPFR_CHANGE_SIGN (atan);
        }
      inexact2 = mpfr_div_2ui (atan, atan, 2, rnd_mode);
      if (MPFR_UNLIKELY (inexact2))
	inexact = inexact2; /* an underflow occurs */
      return inexact;
    }

  supplement = 2 - (comparaison > 0) ? 0 : MPFR_GET_EXP (xp);
  realprec = MPFR_PREC (atan) + MPFR_INT_CEIL_LOG2 (MPFR_PREC (atan)) + 4;
  prec = realprec + supplement + BITS_PER_MP_LIMB;

  MPFR_SAVE_EXPO_MARK (expo);

  mpz_init (ukz);

  /* Initialisation    */
  mpfr_init2 (sk, prec);
  mpfr_init2 (tmp2, prec);
  mpfr_init2 (tmp, prec);
  mpfr_init2 (arctgt, prec);

  for (;;)
    {
      n0 = __gmpfr_ceil_log2 ((double) realprec + supplement + CST);
      MPFR_ASSERTD (3*n0 > 2);
      estimated_delta = 1 + supplement + MPFR_INT_CEIL_LOG2 (3*n0-2);
      prec = realprec + estimated_delta;

      /* Initialisation */
      mpfr_set_prec (sk, prec);
      mpfr_set_prec (tmp2, prec);
      mpfr_set_prec (tmp, prec);
      mpfr_set_prec (arctgt, prec);

      if (comparaison > 0)
	mpfr_ui_div (sk, 1, xp, GMP_RNDN);
      else
	mpfr_set (sk, xp, GMP_RNDN);

      /* sk is 1/|x| if |x| > 1, and |x| otherwise, i.e. min(|x|, 1/|x|) */

      /* Assignation  */
      MPFR_SET_ZERO (arctgt);
      twopoweri = 1;
      for (i = 0; i <= n0; i++)
        {
          /* Calculation of  trunc(tmp) --> mpz */
          mpfr_mul_2ui (tmp, sk, twopoweri, GMP_RNDN);
          mpfr_trunc (tmp, tmp);
          exptol = mpfr_get_z_exp (ukz, tmp);
          /* since the s_k are decreasing (see algorithms.tex),
             and s_0 = min(|x|, 1/|x|) < 1, we have sk < 1,
             thus exptol < 0 */
          MPFR_ASSERTD (exptol < 0);
          mpz_tdiv_q_2exp (ukz, ukz, (unsigned long int) (-exptol));

          /* Calculation of arctan(Ak) */
	  mpfr_set_z (tmp, ukz, GMP_RNDN);
	  mpfr_div_2ui (tmp, tmp, twopoweri, GMP_RNDN);
          mpz_mul (ukz, ukz, ukz);
	  mpz_neg (ukz, ukz);
          mpfr_atan_aux (tmp2, ukz, 2*twopoweri, n0 - i);
          mpfr_mul (tmp2, tmp2, tmp, GMP_RNDN);

          /* Addition and iteration */
          mpfr_add (arctgt, arctgt, tmp2, GMP_RNDN);
          if (MPFR_LIKELY (i < n0))
            {
              mpfr_sub (tmp2, sk, tmp, GMP_RNDN);
              mpfr_mul (sk, sk, tmp, GMP_RNDN);
              mpfr_add_ui (sk, sk, 1, GMP_RNDN);
              mpfr_div (sk, tmp2, sk, GMP_RNDN);
              twopoweri <<= 1;
            }
        }

      if (comparaison > 0)
	{
	  mpfr_const_pi (tmp, GMP_RNDN);
	  mpfr_div_2ui (tmp, tmp, 1, GMP_RNDN);
	  mpfr_sub (arctgt, tmp, arctgt, GMP_RNDN);
	}
      MPFR_SET_POS (arctgt);

      if (mpfr_can_round (arctgt, realprec, GMP_RNDN, GMP_RNDZ,
                          MPFR_PREC (atan) + (rnd_mode == GMP_RNDN)))
	break;
      realprec += BITS_PER_MP_LIMB;
    }

  inexact = mpfr_set4 (atan, arctgt, rnd_mode, sign);

  mpfr_clear (arctgt);
  mpfr_clear (tmp);
  mpfr_clear (tmp2);
  mpfr_clear (sk);

  mpz_clear (ukz);

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (arctgt, inexact, rnd_mode);
}
