/* mpfr_atan -- arc-tangent of a floating-point number

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
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include <stdio.h>

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

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
/* This is the code of 'generic.c' slighty optimized for mpfr_atan */   
static void
mpfr_atan_aux (mpfr_ptr y, mpz_ptr p, long r, int m, mpz_t *tab, 
	       unsigned int *itab)
{
  unsigned long n,i,k,j,l;
  mpz_t *S, *T, *ptoj;
  mp_exp_t diff, expo;
  unsigned int *P2tab;
  unsigned int *S2tab;
  int neg;
  mp_limb_t *d;

  /* Set Tables */
  S    = tab;           /* S */
  ptoj = S + 1*(m+1);   /* p^2^j Precomputed table */
  T    = S + 2*(m+1);   /* Product of Odd integer  table  */
  P2tab = itab;         /* Real p[l] = ptoj[l]*2^P2tab[l] */
  S2tab = itab + m + 1; /* Real s[k] = s[k] * 2^S2tab[k] */

  /* Init S[0] and T[0] */
  mpz_set_ui (S[0], 1);
  S2tab[0] = 0;
  mpz_set_ui (T[0], 1);

  /* Normalize p */
  d = PTR (p);
  for (n = 0 ; MPFR_UNLIKELY (*d == 0) ; d++, n+= BITS_PER_MP_LIMB);
  MPFR_ASSERTD (*d != 0);
  count_trailing_zeros (neg, *d);
  P2tab[0] = n + neg + 1;
  if (n+neg > 0)
    mpz_tdiv_q_2exp (p, p, n+neg);

  /* Extract sign */
  neg = mpz_sgn (p) < 0;
  mpz_abs (p, p);

  /* Check if P==1 (Special case) */
  if (mpz_cmp_ui (p, 1) != 0) {
    /* P!= 1: Precomputed ptoj table */
    mpz_set (ptoj[0], p);
    for (i = 1 ; i < m ; i++) {
      mpz_mul (ptoj[i], ptoj[i-1], ptoj[i-1]);
      P2tab[i] = 2*P2tab[i-1];
    }
    /* Main loop */
    k = 0;
    n = 1UL << m;
    for (i = 1 ; i < n ; i++) {
      k++;
      mpz_set_ui (T[k], 1 + 2*i);
      mpz_set_ui (S[k], 1);
      S2tab[k] = 0;
 
      for (j = i+1, l = 0 ; (j & 1) == 0 ; l++, j>>=1, k--) {
	if (l == 0) {
	  mpz_mul (S[k], ptoj[l], T[k-1]);
	  if (neg)
	    mpz_neg (S[k], S[k]);
	} else {
	  mpz_mul (S[k], S[k], ptoj[l]);
	  mpz_mul (S[k], S[k], T[k-1]);
	}
	mpz_mul (S[k-1], S[k-1], T[k]);
	S2tab[k] += P2tab[l];
	S2tab[k-1] += (r+1)*(1<<l);
	MPFR_ASSERTD (S2tab[k-1] > S2tab[k]);
	mpz_mul_2exp (S[k-1], S[k-1], S2tab[k-1]-S2tab[k]);
	S2tab[k-1] = S2tab[k];
	mpz_add (S[k-1], S[k-1], S[k]);
	mpz_mul (T[k-1], T[k-1], T[k]);
      }
    }
  } else {
    for (i = 1 ; i < m ; i++)
      P2tab[i] = 2*P2tab[i-1];
    k = 0;
    n = 1UL << m;
    for (i = 1 ; i < n ; i++) {
      k++;
      mpz_set_ui (T[k], 1 + 2*i);
      mpz_set_ui (S[k], 1);
      S2tab[k] = 0;
      
      for (j = i+1, l = 0 ; (j & 1) == 0 ; l++, j>>=1, k--) {
	if (l == 0)
	  (neg ? mpz_neg : mpz_set) (S[k], T[k-1]);
	else
	  mpz_mul (S[k], S[k], T[k-1]);
	mpz_mul (S[k-1], S[k-1], T[k]);
	S2tab[k] += P2tab[l];
	S2tab[k-1] += (r+1)*(1<<l);
	MPFR_ASSERTD (S2tab[k-1] > S2tab[k]);
	mpz_mul_2exp (S[k-1], S[k-1], S2tab[k-1]-S2tab[k]);
	S2tab[k-1] = S2tab[k];
	mpz_add (S[k-1], S[k-1], S[k]);
	mpz_mul (T[k-1], T[k-1], T[k]);
      }
    }
  }

  /* printf ("S2tab[0]=%lu S[0]=", S2tab[0]); mpz_out_str (stdout, 2, S[0]); putchar ('\n');
     printf ("T[0]="); mpz_out_str (stdout, 2, T[0]); putchar ('\n'); */

  MPFR_MPZ_SIZEINBASE2 (diff, S[0]);
  diff -= 2*MPFR_PREC (y);
  expo = diff + S2tab[0];
  if (diff >=0)
    mpz_tdiv_q_2exp (S[0], S[0], diff);
  else
    mpz_mul_2exp (S[0], S[0], -diff);

  MPFR_MPZ_SIZEINBASE2 (diff, T[0]);
  diff -= MPFR_PREC (y);
  expo -= (diff + n -1);
  if (diff >= 0)
    mpz_tdiv_q_2exp (T[0], T[0],diff);
  else
    mpz_mul_2exp (T[0], T[0],-diff);

  mpz_tdiv_q (S[0], S[0], T[0]);
  mpfr_set_z (y, S[0], GMP_RNDD);
  MPFR_SET_EXP (y, MPFR_EXP (y) + expo - r*(i-1) );
}

int
mpfr_atan (mpfr_ptr atan, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  mpfr_t xp;
  mpfr_t arctgt, sk, tmp, tmp2;
  mpz_t  ukz;
  int comparaison, inexact, inexact2;
  mp_exp_t estimated_delta;
  mp_prec_t prec, realprec;
  mp_exp_t exptol;
  int i, n0, oldn0;
  unsigned long twopoweri;
  mpz_t *tabz;
  unsigned int *tabi;
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC (("x[%#R]=%R rnd=%d", x, x, rnd_mode),
		 ("atan[%#R]=%R inexact=%d", atan, atan, inexact));

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
  MPFR_TMP_INIT_ABS (xp, x);

  /* Other simple case arctang(-+1)=-+pi/4 */
  comparaison = mpfr_cmp_ui (xp, 1);
  if (MPFR_UNLIKELY (comparaison == 0))
    {
      int neg = MPFR_IS_NEG (x);
      inexact = mpfr_const_pi (atan, MPFR_IS_POS (x) ? rnd_mode
                               : MPFR_INVERT_RND (rnd_mode));
      if (neg)
        {
          inexact = -inexact;
          MPFR_CHANGE_SIGN (atan);
        }
      inexact2 = mpfr_div_2ui (atan, atan, 2, rnd_mode);
      if (MPFR_UNLIKELY (inexact2))
	inexact = inexact2; /* an underflow occurs */
      return inexact;
    }

  realprec = MPFR_PREC (atan) + MPFR_INT_CEIL_LOG2 (MPFR_PREC (atan)) + 4;
  if (MPFR_PREC (atan) + 5 > MPFR_PREC (x) && MPFR_GET_EXP (x) < 0)
    {
      mpfr_uexp_t ue = (mpfr_uexp_t) (-MPFR_GET_EXP (x));
      if (realprec < 2*ue + 5)
	realprec = 2*ue + 5;
    }
  prec = realprec + BITS_PER_MP_LIMB;

  MPFR_SAVE_EXPO_MARK (expo);

  /* Initialisation    */
  mpz_init (ukz);
  mpfr_init2 (sk, prec);
  mpfr_init2 (tmp2, prec);
  mpfr_init2 (tmp, prec);
  mpfr_init2 (arctgt, prec);
  oldn0 = 0;
  tabz = NULL;
  tabi = NULL;

  MPFR_ZIV_INIT (loop, prec);
  for (;;)
    {
      /* n0 = ceil(log2(realprec + 2 + 1+ln(2.4)/ln(2))) */
      n0 = MPFR_INT_CEIL_LOG2 (realprec + 2 + 3);
      MPFR_ASSERTD (3*n0 > 2);
      estimated_delta = 1 + 2 + MPFR_INT_CEIL_LOG2 (3*n0-2);
      prec = realprec + estimated_delta;

      /* Initialisation */
      mpfr_set_prec (sk, prec);
      mpfr_set_prec (tmp2, prec);
      mpfr_set_prec (tmp, prec);
      mpfr_set_prec (arctgt, prec);
      if (oldn0 < 3*n0+1) {
	if (MPFR_LIKELY (oldn0 == 0)) {
	  tabz = (mpz_t *) (*__gmp_allocate_func) (3*(n0+1)*sizeof (mpz_t));
	  tabi = (unsigned int *) (*__gmp_allocate_func) 
	    (2*(n0+1)*sizeof(unsigned int));
	} else {
	  tabz = (mpz_t *) (*__gmp_reallocate_func) 
	    (tabz, oldn0*sizeof (mpz_t), 3*(n0+1)*sizeof (mpz_t));
	  tabi = (unsigned int *) (*__gmp_reallocate_func)
	    (tabi, oldn0/3*2*sizeof (unsigned int),
	     2*(n0+1)*sizeof(unsigned int));
	}
	for (i = oldn0 ; i < 3*(n0+1) ; i++)
	  mpz_init (tabz[i]);
	oldn0 = 3*(n0+1);
      }
      
      if (comparaison > 0)
	mpfr_ui_div (sk, 1, xp, GMP_RNDN);
      else
	mpfr_set (sk, xp, GMP_RNDN);

      /* sk is 1/|x| if |x| > 1, and |x| otherwise, i.e. min(|x|, 1/|x|) */

      /* Assignation  */
      MPFR_SET_ZERO (arctgt);
      twopoweri = 1<<0;
      MPFR_ASSERTD (n0 >= 4);
      for (i = 0 ; i <= n0; i++)
        {
          /* Calculation of  trunc(tmp) --> mpz */
          mpfr_mul_2ui (tmp, sk, twopoweri, GMP_RNDN);
          mpfr_trunc (tmp, tmp);
	  if (!MPFR_IS_ZERO (tmp))
	    {
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
	      mpfr_atan_aux (tmp2, ukz, 2*twopoweri, n0 - i, tabz, tabi);
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
	  else
	    twopoweri <<= 1;
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
      MPFR_ZIV_NEXT (loop, realprec);
    }
  MPFR_ZIV_FREE (loop);

  inexact = mpfr_set4 (atan, arctgt, rnd_mode, MPFR_SIGN (x));

  for (i = 0 ; i < oldn0 ; i++)
    mpz_clear (tabz[i]);
  (*__gmp_free_func) (tabz, oldn0*sizeof (mpz_t));
  (*__gmp_free_func) (tabi, oldn0/3*2*sizeof (unsigned int));

  mpfr_clear (arctgt);
  mpfr_clear (tmp);
  mpfr_clear (tmp2);
  mpfr_clear (sk);

  mpz_clear (ukz);

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (arctgt, inexact, rnd_mode);
}
