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
#define CST2  1.45  /* CST2=1/ln(2) */

static int mpfr_atan_aux _MPFR_PROTO((mpfr_ptr, mpz_srcptr, long, int));

#undef B
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

int
mpfr_atan (mpfr_ptr arctangent, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  mpfr_t Pisur2;
  mpfr_t xp;
  mpfr_t arctgt;

  int comparaison, sign, supplement, inexact;

  mpfr_t t_arctan;
  int i;
  mpz_t ukz;
  mpfr_t ukf;
  mpfr_t sk,Ak;
  mpz_t square;
  mpfr_t tmp_arctan;
  mpfr_t tmp, tmp2;
#ifdef DEBUG
  mpfr_t tst;
#endif
  int twopoweri;
  int Prec;
  int prec_x;
  int prec_arctan;
  int realprec;
  int estimated_delta;
  /* calculation of the floor */
  mp_exp_t exptol;

  int N0;
  int logn;

  /* Trivial cases */
  if (MPFR_UNLIKELY( MPFR_IS_SINGULAR(x) ))
    {
      if (MPFR_IS_NAN(x))
	{
	  MPFR_SET_NAN(arctangent);
	  MPFR_RET_NAN;
	}
      else if (MPFR_IS_INF(x))
	{
	  MPFR_CLEAR_FLAGS(arctangent);
	  if (MPFR_IS_POS(x))
	    /* arctan(+inf) = Pi/2 */
	    inexact = mpfr_const_pi (arctangent, rnd_mode);
	  else 
	    /* arctan(-inf) = -Pi/2 */
	    {
	      if (rnd_mode == GMP_RNDU)
		rnd_mode = GMP_RNDD;
	      else if (rnd_mode == GMP_RNDD)
		rnd_mode = GMP_RNDU;
	      inexact = -mpfr_const_pi (arctangent, rnd_mode);
	      MPFR_CHANGE_SIGN (arctangent);
	    }
	  MPFR_SET_EXP (arctangent, MPFR_GET_EXP (arctangent) - 1);
	  return inexact;
	}
      else /* x is necessarily 0 */
	{
          MPFR_ASSERTD(MPFR_IS_ZERO(x));
	  mpfr_set_ui (arctangent, 0, GMP_RNDN);
	  return 0; /* exact result */
	}
    }
  MPFR_CLEAR_FLAGS(arctangent);

  sign = MPFR_SIGN(x);
  prec_arctan = MPFR_PREC(arctangent);

  /* Set x_p=|x| */
  mpfr_init2 (xp, MPFR_PREC(x));
  mpfr_abs (xp, x, rnd_mode);

  /* Other simple case arctang(-+1)=-+pi/4 */
  comparaison = mpfr_cmp_ui (xp, 1);
  if (MPFR_UNLIKELY(comparaison == 0))
    {
      inexact = mpfr_const_pi (arctangent, MPFR_IS_POS_SIGN(sign) ? rnd_mode
                               : MPFR_INVERT_RND(rnd_mode));
      MPFR_SET_EXP (arctangent, MPFR_GET_EXP (arctangent) - 2);
      if (MPFR_IS_NEG_SIGN( sign ))
        {
          inexact = -inexact;
          MPFR_CHANGE_SIGN(arctangent);
        }
      mpfr_clear (xp);
      return inexact;
    }

  if (comparaison > 0)
    supplement = 2;
  else
    supplement = 2 - MPFR_GET_EXP (xp);

  mpfr_save_emin_emax ();

  /* FIXME: Could use MPFR_INT_CEIL_LOG2? */
  prec_x = __gmpfr_ceil_log2 ((double) MPFR_PREC(x) / BITS_PER_MP_LIMB);
  logn   = MPFR_INT_CEIL_LOG2 ( prec_x);
  if (logn < 2) 
    logn = 2;
  realprec = prec_arctan + MPFR_INT_CEIL_LOG2 (prec_arctan) + 4;

  mpz_init (ukz);
  mpz_init (square);

  /* Initialisation    */
  mpfr_init (sk);
  mpfr_init (ukf);
  mpfr_init (t_arctan);
  mpfr_init (tmp_arctan);
  mpfr_init (tmp);
  mpfr_init (tmp2);
  mpfr_init (Ak);
  mpfr_init (arctgt);
  mpfr_init (Pisur2);

  while (1)
    {
      N0 = MPFR_INT_CEIL_LOG2 (realprec + supplement + CST);
      estimated_delta = 1 + supplement + MPFR_INT_CEIL_LOG2 (3*N0-2);
      Prec = realprec + estimated_delta;

      /* Initialisation    */
      mpfr_set_prec (sk,Prec);
      mpfr_set_prec (ukf, Prec);
      mpfr_set_prec (t_arctan, Prec);
      mpfr_set_prec (tmp_arctan, Prec);
      mpfr_set_prec (tmp, Prec);
      mpfr_set_prec (tmp2, Prec);
      mpfr_set_prec (Ak, Prec);
      mpfr_set_prec (arctgt, Prec);

      if (comparaison > 0)
        {
          mpfr_set_prec (Pisur2, Prec);
          mpfr_const_pi (Pisur2, GMP_RNDN);
          mpfr_div_2ui (Pisur2, Pisur2, 1, GMP_RNDN);
          mpfr_ui_div (sk, 1, xp, GMP_RNDN);
        }
      else
	mpfr_set (sk, xp, GMP_RNDN);

      /* sk is 1/|x| if |x| > 1, and |x| otherwise, i.e. min(|x|, 1/|x|) */

      /* Assignation  */
      mpfr_set_ui (tmp_arctan, 0, GMP_RNDN);
      twopoweri = 1;
      for (i = 0; i <= N0; i++)
        {
          mpfr_mul_2ui (tmp, sk, twopoweri, GMP_RNDN);
          /* Calculation of  trunc(tmp) --> mpz */
          mpfr_trunc (ukf, tmp);
          exptol = mpfr_get_z_exp (ukz, ukf);
          /* since the s_k are decreasing (see algorithms.tex),
             and s_0 = min(|x|, 1/|x|) < 1, we have sk < 1,
             thus exptol < 0 */
          MPFR_ASSERTD(exptol < 0);
          mpz_tdiv_q_2exp (ukz, ukz, (unsigned long int) (-exptol));

          /* Calculation of arctan(Ak) */
          mpz_mul (square, ukz, ukz);
	  mpz_neg (square, square);
          mpfr_atan_aux (t_arctan, square, 2*twopoweri, N0 - i);
          mpfr_set_z (Ak, ukz, GMP_RNDN);
          mpfr_div_2ui (Ak, Ak, twopoweri, GMP_RNDN);
          mpfr_mul (t_arctan, t_arctan, Ak, GMP_RNDN);

          /* Addition and iteration */
          mpfr_add (tmp_arctan, tmp_arctan, t_arctan, GMP_RNDN);
          if (i < N0)
            {
              mpfr_sub (tmp, sk, Ak, GMP_RNDN);
              mpfr_mul (tmp2, sk, Ak, GMP_RNDN);
              mpfr_add_ui (tmp2, tmp2, 1, GMP_RNDN);
              mpfr_div (sk, tmp, tmp2, GMP_RNDN);
              twopoweri <<= 1;
            }
        }

      if (comparaison > 0)
	mpfr_sub(arctgt, Pisur2, tmp_arctan, GMP_RNDN);
      else
	mpfr_set(arctgt, tmp_arctan, GMP_RNDN);
      MPFR_SET_POS(arctgt);

      if (!mpfr_can_round (arctgt, realprec, GMP_RNDN, GMP_RNDZ,
                          MPFR_PREC (arctangent) + (rnd_mode == GMP_RNDN)))
	realprec += MPFR_INT_CEIL_LOG2 (realprec);
      else
	break;
    }

  inexact = MPFR_IS_POS_SIGN(sign) ? mpfr_set (arctangent, arctgt, rnd_mode)
    : mpfr_neg (arctangent, arctgt, rnd_mode);

  mpfr_clear (sk);
  mpfr_clear (ukf);
  mpfr_clear (t_arctan);
  mpfr_clear (tmp_arctan);
  mpfr_clear (tmp);
  mpfr_clear (tmp2);
  mpfr_clear (Ak);
  mpfr_clear (arctgt);

  mpfr_clear (Pisur2);

  mpfr_clear (xp);
  mpz_clear (ukz);
  mpz_clear (square);

  mpfr_restore_emin_emax ();
  return mpfr_check_range (arctgt, inexact, rnd_mode);
}
