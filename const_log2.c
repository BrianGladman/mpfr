/* mpfr_const_log2 -- compute natural logarithm of 2

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
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

static int mpfr_aux_log2 (mpfr_ptr, mpz_srcptr, long, int);
static int mpfr_const_aux_log2 (mpfr_ptr, mp_rnd_t);

#define A
#define A1 1
#define A2 1
#undef B
#define C
#define C1 2
#define C2 1
#define NO_FACTORIAL
#undef R_IS_RATIONAL
#define GENERIC mpfr_aux_log2
#include "generic.c"
#undef A
#undef A1
#undef A2
#undef NO_FACTORIAL
#undef GENERIC
#undef C
#undef C1
#undef C2

static int
mpfr_const_aux_log2 (mpfr_ptr mylog, mp_rnd_t rnd_mode)
{
  mp_prec_t prec;
  mpfr_t tmp1, tmp2, result,tmp3;
  mpz_t cst;
  int logn;
  mp_prec_t prec_i_want;
  mp_prec_t prec_x;
  int inexact = 0;  /* here, 0 means not set */

  mpz_init (cst);
  prec_i_want = MPFR_PREC(mylog);
  logn = MPFR_INT_CEIL_LOG2 (prec_i_want);
  prec_x = prec_i_want + logn;
  while (!inexact)
    {
      prec = MPFR_INT_CEIL_LOG2 (prec_x);
      mpfr_init2 (tmp1, prec_x);
      mpfr_init2 (result, prec_x);
      mpfr_init2 (tmp2, prec_x);
      mpfr_init2 (tmp3, prec_x);
      mpz_set_ui (cst, 1);
      mpfr_aux_log2 (tmp1, cst, 4, prec-2);
      mpfr_div_2ui (tmp1, tmp1, 4, GMP_RNDD);
      mpfr_mul_ui (tmp1, tmp1, 15, GMP_RNDD);

      mpz_set_ui (cst, 3);
      mpfr_aux_log2 (tmp2, cst, 7, prec-2);
      mpfr_div_2ui (tmp2, tmp2, 7, GMP_RNDD);
      mpfr_mul_ui (tmp2, tmp2, 5*3, GMP_RNDD);
      mpfr_sub (result, tmp1, tmp2, GMP_RNDD);

      mpz_set_ui (cst, 13);
      mpfr_aux_log2 (tmp3, cst, 8, prec-2);
      mpfr_div_2ui (tmp3, tmp3, 8, GMP_RNDD);
      mpfr_mul_ui (tmp3, tmp3, 3*13, GMP_RNDD);
      mpfr_sub (result, result, tmp3, GMP_RNDD);

      mpfr_clear (tmp1);
      mpfr_clear (tmp2);
      mpfr_clear (tmp3);
      if (mpfr_can_round (result, prec_x - 2, GMP_RNDD, GMP_RNDZ,
                          prec_i_want + (rnd_mode == GMP_RNDN)))
        {
          inexact = mpfr_set (mylog, result, rnd_mode);
          MPFR_ASSERTN (inexact != 0);
        }
      else
        {
          prec_x += logn;
        }

      mpfr_clear (result);
    }
  mpz_clear (cst);
  return inexact;
}

/* Cross-over point from nai"ve Taylor series to binary splitting,
   obtained experimentally on a Pentium II. Optimal value for
   target machine should be determined by tuneup. */
#define LOG2_THRESHOLD 25000

/* set x to log(2) rounded to precision MPFR_PREC(x) with direction rnd_mode

   use formula log(2) = sum(1/k/2^k, k=1..infinity)

   whence 2^N*log(2) = S(N) + R(N)

   where S(N) = sum(2^(N-k)/k, k=1..N-1)
   and   R(N) = sum(1/k/2^(k-N), k=N..infinity) < 2/N

   Let S'(N) = sum(floor(2^(N-k)/k), k=1..N-1)

   Then 2^N*log(2)-S'(N) <= N-1+2/N <= N for N>=2.
*/
int
(mpfr_const_log2) (mpfr_ptr x, mp_rnd_t rnd_mode)
{
  mp_prec_t N, k, precx;
  mpz_t s, t, u;
  int inexact;
  MPFR_SAVE_EXPO_DECL (expo);

  MPFR_SAVE_EXPO_MARK (expo);
  precx = MPFR_PREC(x);

  /* need to recompute */
  if (precx < LOG2_THRESHOLD) /* use nai"ve Taylor series evaluation */
    {
      /* the following was checked by exhaustive search to give a correct
         result for all 4 rounding modes up to precx = 13500 */
      N = precx + 2 * MPFR_INT_CEIL_LOG2 (precx) + 1;

      mpz_init (s); /* set to zero */
      mpz_init (u);
      mpz_init_set_ui (t, 1);

      /* use log(2) = sum((6*k-1)/(2*k^2-k)/2^(2*k+1), k=1..infinity) */
      mpz_mul_2exp (t, t, N-1);
      for (k=1; k<=N/2; k++)
        {
          mpz_div_2exp (t, t, 2);
          mpz_mul_ui (u, t, 6*k-1);
          mpz_fdiv_q_ui (u, u, k*(2*k-1));
          mpz_add (s, s, u);
        }

      inexact = mpfr_set_z (x, s, rnd_mode); /* Can overflow => save_emin */
      MPFR_SET_EXP (x, MPFR_GET_EXP (x) - N);
      mpz_clear (s);
      mpz_clear (t);
      mpz_clear (u);
    }
  else /* use binary splitting method */
    inexact = mpfr_const_aux_log2 (x, rnd_mode);

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (x, inexact, rnd_mode);
}
