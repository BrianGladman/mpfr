/* mpfr_ai2 -- implementation of Airy Ai function by Smith algorithm.

Copyright 2010 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

#define SMALL_PRECISION 32

static unsigned long isqrt (unsigned long n)
{
  mpz_t tmp;
  unsigned int res;
  mpz_init (tmp);
  mpz_set_ui (tmp, n);
  mpz_sqrt (tmp, tmp);
  res = mpz_get_ui (tmp);
  mpz_clear (tmp);
  return res;
}

/* Reminder and notations:
   -----------------------

   Ai is the solution of:
   / y'' - x*y = 0
   {    Ai (0)   = 1/ ( 9^(1/3) * Gamma (2/3) )
   \  Ai' (0)   = -1/ ( 3^(1/3) * Gamma (1/3) )

   Series development:
   Ai (x) = sum (a_i*x^i)
   = sum (t_i)

   Recurrences:
   a_(i+3) = a_i / ((i+2)*(i+3))
   t_(i+3) = t_i * x^3 / ((i+2)*(i+3))

   Values:
   a_0 = Ai (0)  ~  0.355
   a_1 = Ai' (0) ~ -0.259
*/


/* Airy function Ai evaluated by Smith algorithm */
int
mpfr_ai2 (mpfr_ptr y, mpfr_srcptr x, mpfr_rnd_t rnd)
{
  MPFR_ZIV_DECL (loop);
  MPFR_SAVE_EXPO_DECL (expo);
  mpfr_prec_t wprec;             /* working precision */
  mpfr_prec_t prec;              /* target precision */
  mpfr_prec_t err;               /* used to estimate the evaluation error */
  mpfr_prec_t correctBits;       /* estimates the number of correct bits*/
  unsigned long int i, j, L, t;
  unsigned long int cond;        /* condition number of the series */
  unsigned long int assumed_exponent;  /* used as a lowerbound of |EXP (Ai (x))| */
  int r;                         /* returned ternary value */
  mpfr_t s;                      /* used to store the partial sum */
  mpfr_t u0, u1;
  mpfr_t *z;                     /* used to store the (x^3j) */
  mpfr_t result;
  mpfr_t tmp_sp, tmp2_sp;        /* small precision variables */
  unsigned long int x3u;         /* used to store ceil (x^3) */
  mpfr_t temp1, temp2;
  int test0, test1;

  /* Logging */
  MPFR_LOG_FUNC ( ("x[%#R]=%R rnd=%d", x, x, rnd), ("y[%#R]=%R", y, y) );

  /* Special cases */
  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    {
      if (MPFR_IS_NAN (x))
        {
          MPFR_SET_NAN (y);
          MPFR_RET_NAN;
        }
      else if (MPFR_IS_INF (x))
        return mpfr_set_ui (y, 0, rnd);
    }

  /* Save current exponents range */
  MPFR_SAVE_EXPO_MARK (expo);

  /* FIXME: underflow for large values of |x| */


  /* Set initial precision */
  /* See the analysis for the naive evaluation */

  /* We begin with 11 guard bits */
  prec = MPFR_PREC (y) + 11;
  MPFR_ZIV_INIT (loop, prec);

  mpfr_init2 (tmp_sp, SMALL_PRECISION);
  mpfr_init2 (tmp2_sp, SMALL_PRECISION);
  mpfr_abs (tmp_sp, x, MPFR_RNDU);
  mpfr_pow_ui (tmp_sp, tmp_sp, 3, MPFR_RNDU);
  mpfr_sqrt (tmp_sp, tmp_sp, MPFR_RNDU); /* tmp_sp ~ x^3/2 */

  /* 0.96179669392597567 >~ 2/3 * log2(e). See algorithms.tex */
  mpfr_set_str(tmp2_sp, "0.96179669392597567", 10, MPFR_RNDU);
  mpfr_mul (tmp2_sp, tmp_sp, tmp2_sp, MPFR_RNDU);

  /* cond represents the number of lost bits in the evaluation of the sum */
  if ( (MPFR_IS_ZERO(x)) || (MPFR_GET_EXP (x) <= 0) )
    cond = 0;
  else
    cond = mpfr_get_ui (tmp2_sp, MPFR_RNDU) - (MPFR_GET_EXP (x) - 1)/4 - 1;

  /* This variable is used to store the maximal assumed exponent of       */
  /* Ai (x). More precisely, we assume that |Ai (x)| will be greater than */
  /* 2^{-assumedExp}.                                                     */
  if (MPFR_IS_ZERO(x)) assumed_exponent = 2;
  else 
    {
      if (MPFR_IS_POS (x))
	{
	  if (MPFR_GET_EXP (x) <= 0)
	    assumed_exponent = 3;
	  else
	    assumed_exponent = 2 + (MPFR_GET_EXP (x)/4 + 1) + mpfr_get_ui (tmp2_sp, MPFR_RNDU);
	}
      /* We do not know Ai (x) yet */
      /* We cover the case when EXP (Ai (x))>=-10 */
      else
	assumed_exponent = 10;
    }

  wprec = prec + MPFR_INT_CEIL_LOG2 (prec) + 6 + cond + assumed_exponent;

  /* We assume that the truncation rank will be ~ prec */
  L = isqrt (prec);
  MPFR_LOG_MSG (("size of blocks L = %lu\n", L));

  z = (mpfr_t *) (*__gmp_allocate_func) ( (L + 1) * sizeof (mpfr_t) );
  MPFR_ASSERTN (z != NULL);
  for (j=0; j<=L; j++)
    mpfr_init (z[j]);

  mpfr_init (s);
  mpfr_init (u0); mpfr_init (u1);
  mpfr_init (result);
  mpfr_init (temp1);
  mpfr_init (temp2);

  /* ZIV loop */
  for (;;)
    {
      MPFR_LOG_MSG (("working precision: %Pu\n", wprec));

      for (j=0; j<=L; j++)
        mpfr_set_prec (z[j], wprec);
      mpfr_set_prec (s, wprec);
      mpfr_set_prec (u0, wprec); mpfr_set_prec (u1, wprec);
      mpfr_set_prec (result, wprec);

      mpfr_set_ui (u0, 1, MPFR_RNDN);
      mpfr_set (u1, x, MPFR_RNDN);

      mpfr_set_ui (z[0], 1, MPFR_RNDU);
      mpfr_sqr (z[1], u1, MPFR_RNDU);
      mpfr_mul (z[1], z[1], x, (MPFR_IS_POS (x) ? MPFR_RNDU : MPFR_RNDD) );

      if (MPFR_IS_NEG (x))
        MPFR_CHANGE_SIGN (z[1]);
      x3u = mpfr_get_ui (z[1], MPFR_RNDU);   /* x3u >= ceil (x^3) */
      if (MPFR_IS_NEG (x))
        MPFR_CHANGE_SIGN (z[1]);

      for (j=2; j<=L ;j++)
        {
          if (j%2 == 0)
            mpfr_sqr (z[j], z[j/2], MPFR_RNDN);
          else
            mpfr_mul (z[j], z[j-1], z[1], MPFR_RNDN);
        }

      mpfr_gamma_one_and_two_third (temp1, temp2, wprec);
      mpfr_set_ui (u0, 9, MPFR_RNDN);
      mpfr_cbrt (u0, u0, MPFR_RNDN);
      mpfr_mul (u0, u0, temp2, MPFR_RNDN);
      mpfr_ui_div (u0, 1, u0 , MPFR_RNDN); /* u0 = 1 / ( Gamma (2/3)*9^(1/3) ) */

      mpfr_set_ui (u1, 3, MPFR_RNDN);
      mpfr_cbrt (u1, u1, MPFR_RNDN);
      mpfr_mul (u1, u1, temp1, MPFR_RNDN);
      mpfr_neg (u1, u1, MPFR_RNDN);
      mpfr_div (u1, x, u1, MPFR_RNDN); /* u1 = -x/(Gamma (1/3)*3^(1/3)) */

      mpfr_set_ui (result, 0, MPFR_RNDN);
      t = 0;

      /* Evaluation of the series by Smith' method    */
      for (i=0; ; i++)
        {
          t += 3 * L;

          /* k = 0 */
          t -= 3;
          mpfr_set (s, z[L-1], MPFR_RNDN);
          for (j=L-2; ; j--)
            {
              t -= 3;
              mpfr_div_ui2 (s, s, (t+2), (t+3), MPFR_RNDN);
              mpfr_add (s, s, z[j], MPFR_RNDN);
              if (j==0)
                break;
            }
          mpfr_mul (s, s, u0, MPFR_RNDN);
          mpfr_add (result, result, s, MPFR_RNDN);

          mpfr_mul (u0, u0, z[L], MPFR_RNDN);
          for (j=0; j<=L-1; j++)
            {
              mpfr_div_ui2 (u0, u0, (t + 2), (t + 3), MPFR_RNDN);
              t += 3;
            }

          t++;

          /* k = 1 */
          t -= 3;
          mpfr_set (s, z[L-1], MPFR_RNDN);
          for (j=L-2; ; j--)
            {
              t -= 3;
              mpfr_div_ui2 (s, s, (t + 2), (t + 3), MPFR_RNDN);
              mpfr_add (s, s, z[j], MPFR_RNDN);
              if (j==0)
                break;
            }
          mpfr_mul (s, s, u1, MPFR_RNDN);
          mpfr_add (result, result, s, MPFR_RNDN);

          mpfr_mul (u1, u1, z[L], MPFR_RNDN);
          for (j=0; j<=L-1; j++)
            {
              mpfr_div_ui2 (u1, u1, (t + 2), (t + 3), MPFR_RNDN);
              t += 3;
            }

          t++;

          /* k = 2 */
          t++;

          /* End of the loop over k */
          t -= 3;

          test0 = MPFR_IS_ZERO (u0) ||
            (MPFR_GET_EXP (u0) + (mp_exp_t)prec + 4 <= MPFR_GET_EXP (result));
          test1 = MPFR_IS_ZERO (u1) ||
            (MPFR_GET_EXP (u1) + (mp_exp_t)prec + 4 <= MPFR_GET_EXP (result));

          if ( test0 && test1 && (x3u <= (t + 2) * (t + 3) / 2) )
            break;
        }

      MPFR_LOG_MSG (("Truncation rank: %lu\n", t));

      err = 5 + MPFR_INT_CEIL_LOG2 (L+1) + MPFR_INT_CEIL_LOG2 (i+1) + cond - MPFR_GET_EXP (result);

      /* err is the number of bits lost due to the evaluation error */
      /* wprec-(prec+1): the number of bits lost due to the approximation error */
      MPFR_LOG_MSG (("Roundoff error: %Pu\n", err));
      MPFR_LOG_MSG (("Approxim error: %Pu\n", wprec - prec - 1));

      if (wprec < err+1)
        correctBits = 0;
      else
        {
          if (wprec < err+prec+1)
            correctBits = wprec - err - 1;
          else
            correctBits = prec;
        }

      if (MPFR_LIKELY (MPFR_CAN_ROUND (result, correctBits, MPFR_PREC (y), rnd)))
        break;

      for (j=0; j<=L; j++)
        mpfr_clear (z[j]);
      (*__gmp_free_func) (z, (L + 1) * sizeof (mpfr_t));
      L = isqrt (t);
      MPFR_LOG_MSG (("size of blocks L = %lu\n", L));
      z = (mpfr_t *) (*__gmp_allocate_func) ( (L + 1) * sizeof (mpfr_t));
      MPFR_ASSERTN (z != NULL);
      for (j=0; j<=L; j++)
        mpfr_init (z[j]);

      if (correctBits == 0)
        {
          assumed_exponent *= 2;
          MPFR_LOG_MSG (("Not a single bit correct (assumed_exponent=%lu)\n", assumed_exponent));
          wprec = prec + 6 + MPFR_INT_CEIL_LOG2 (t) + cond + assumed_exponent;
        }
    else
      {
        if (correctBits < prec)
          { /* The precision was badly chosen */
            MPFR_LOG_MSG (("Bad assumption on the exponent of Ai (x) (E=%ld)\n", (long)(MPFR_GET_EXP (result))));
            wprec = prec + err + 1;
          }
        else
          { /* We are really in a bad case of the TMD */
            MPFR_ZIV_NEXT (loop, prec);

            /* We update wprec */
            /* We assume that t will not be multiplied by more than 4 */
            wprec = prec + (MPFR_INT_CEIL_LOG2 (t) + 2) + 6 + cond - MPFR_GET_EXP (result);
          }
      }
    } /* End of ZIV loop */

  MPFR_ZIV_FREE (loop);
  MPFR_SAVE_EXPO_FREE (expo);

  r = mpfr_set (y, result, rnd);

  mpfr_clear (tmp_sp);
  mpfr_clear (tmp2_sp);
  for (j=0; j<=L; j++)
    mpfr_clear (z[j]);
  (*__gmp_free_func) (z, (L + 1) * sizeof (mpfr_t));

  mpfr_clear (s);
  mpfr_clear (u0); mpfr_clear (u1);
  mpfr_clear (result);
  mpfr_clear (temp1);
  mpfr_clear (temp2);

  return r;
}
