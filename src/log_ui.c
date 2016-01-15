/* mpfr_log_ui -- compute natural logarithm of an unsigned long

Copyright 2014-2016 Free Software Foundation, Inc.
Contributed by the AriC and Caramel projects, INRIA.

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

/* Auxiliary function: Compute using binary splitting the sum
   -sum((-x)^(n+1-n1)/n, n = n1..n2-1) where x = p/2^k,
   -1/3 <= x <= 1/3. For n1=1 we get -sum((-x)^n/n, n = 1..n2-1).
   Numerator is P[0], denominator is Q[0],
   Need 1+ceil(log(n2-n1)/log(2)) cells in P[],Q[].
*/
static void
S (mpz_t *P, mpz_t *Q, unsigned long n1, unsigned long n2,
   long p, unsigned long k)
{
  if (n2 == n1 + 1)
    {
      mpz_set_si (P[0], p);
      mpz_set_ui (Q[0], n1);
      mpz_mul_2exp (Q[0], Q[0], k);
    }
  else
    {
      unsigned long m = (n1 / 2) + (n2 / 2) + (n1 & 1UL & n2);
      /* m = floor((n1+n2)/2) */

      S (P, Q, n1, m, p, k);
      S (P + 1, Q + 1, m, n2, p, k);
      /* P1/Q1 + (-p)^(m-n1)/2^(k(m-n1)) P2/Q2 =
         (P1*2^(k(m-n1))*Q2 + Q1*(-p)^(m-n1)*P2)/(Q1*Q2*2^(k(m-n1)))
         We know 2^(k(m-n1)) divides Q1, thus:
         (P1*Q2 + (Q1/2^(k(m-n1)))*(-p)^(m-n1)*P2)/(Q1*Q2) */
      mpz_mul (P[0], P[0], Q[1]);       /* P1*Q2 */
      mpz_mul (Q[1], Q[0], Q[1]);       /* Q1*Q2 */
      mpz_tdiv_q_2exp (Q[0], Q[0], k * (m - n1)); /* Q1/2^(k(m-n1)) */
      mpz_mul (P[1], P[1], Q[0]);       /* Q1*P2/2^(k(m-n1)) */
      mpz_swap (Q[0], Q[1]);            /* Q1*Q2 */
      if (p > 0)
        {
          mpz_ui_pow_ui (Q[1], p, m - n1);       /* p^m */
          if ((m - n1) & 1)
            mpz_neg (Q[1], Q[1]);
        }
      else
        mpz_ui_pow_ui (Q[1], -p, m - n1);       /* (-p)^m */
      mpz_mul (P[1], P[1], Q[1]);       /* Q1*Q2*2^(km) */
      mpz_add (P[0], P[0], P[1]);
    }
}

int
mpfr_log_ui (mpfr_ptr x, unsigned long n, mpfr_rnd_t rnd_mode)
{
  unsigned long k;
  mpfr_prec_t w; /* working precision */
  mpz_t *P, *Q;
  mpfr_t t, q;
  int inexact;
  unsigned long N, lgN, i;
  long p;
  MPFR_GROUP_DECL(group);
  MPFR_TMP_DECL(marker);
  MPFR_ZIV_DECL(loop);
  MPFR_SAVE_EXPO_DECL (expo);

  if (n <= 2)
    {
      if (n == 0)
        {
          MPFR_SET_INF (x);
          MPFR_SET_NEG (x);
          MPFR_SET_DIVBY0 ();
          MPFR_RET (0); /* log(0) is an exact -infinity */
        }
      else if (n == 1)
        {
          MPFR_SET_ZERO (x);
          MPFR_SET_POS (x);
          MPFR_RET (0); /* only "normal" case where the result is exact */
        }
      /* now n=2 */
      return mpfr_const_log2 (x, rnd_mode);
    }

  /* here n >= 3 */

  /* argument reduction: compute k such that 2/3 <= n/2^k < 4/3,
     i.e., 2^(k+1) <= 3n < 2^(k+2) */

  k = __gmpfr_ceil_log2 (3.0 * (double) n) - 2;
  MPFR_ASSERTD (k >= 2);

  /* the reduced argument is n/2^k - 1 = (n-2^k)/2^k */
  p = (long) n - (1L << k);  /* FIXME: integer overflow for large n */

  MPFR_TMP_MARK(marker);
  w = MPFR_PREC(x) + 10;
  MPFR_GROUP_INIT_2(group, w, t, q);
  MPFR_SAVE_EXPO_MARK (expo);

  MPFR_ZIV_INIT (loop, w);
  for (;;)
    {
      mpfr_t tmp;
      unsigned int err;
      /* we need at most w*log(2)/log(3) terms for an accuracy of w bits */
      mpfr_init2 (tmp, 32);
      /* 1354911329/2^31 is a 32-bit upper bound for log(2)/log(3) */
      mpfr_set_ui_2exp (tmp, 1354911329, -31, MPFR_RNDU);
      mpfr_mul_ui (tmp, tmp, w, MPFR_RNDU);
      N = mpfr_get_ui (tmp, MPFR_RNDU);
      lgN = MPFR_INT_CEIL_LOG2 (N) + 1;
      mpfr_clear (tmp);
      P = (mpz_t *) MPFR_TMP_ALLOC (2 * lgN * sizeof (mpz_t));
      Q = P + lgN;
      for (i = 0; i < lgN; i++)
        {
          mpz_init (P[i]);
          mpz_init (Q[i]);
        }
      S (P, Q, 1, N, p, k);
      mpfr_set_z (t, P[0], MPFR_RNDN); /* t = P[0] * (1 + theta_1) */
      mpfr_set_z (q, Q[0], MPFR_RNDN); /* q = Q[0] * (1 + theta_2) */
      mpfr_div (t, t, q, MPFR_RNDN);   /* t = P[0]/Q[0] * (1 + theta_3)^3
                                            = log(n/2^k) * (1 + theta_4)^4
                                            for |theta_i| < 2^(-w) */
      /* argument reconstruction: add k*log(2) */
      mpfr_const_log2 (q, MPFR_RNDN);
      mpfr_mul_ui (q, q, k, MPFR_RNDN);
      mpfr_add (t, t, q, MPFR_RNDN);
      for (i = 0; i < lgN; i++)
        {
          mpz_clear (P[i]);
          mpz_clear (Q[i]);
        }
      /* the maximal error is 5 ulps for P/Q, since |(1+/-u)^4 - 1| < 5*u
         for u < 2^(-12), k ulps for k*log(2), and 1 ulp for the addition,
         thus at most k+6 ulps */
      err = MPFR_INT_CEIL_LOG2 (k + 6);
      if (MPFR_LIKELY (MPFR_CAN_ROUND (t, w - err, MPFR_PREC(x), rnd_mode)))
        break;

      MPFR_ZIV_NEXT (loop, w);
      MPFR_GROUP_REPREC_2(group, w, t, q);
    }
  MPFR_ZIV_FREE (loop);

  inexact = mpfr_set (x, t, rnd_mode);

  MPFR_GROUP_CLEAR(group);
  MPFR_TMP_FREE(marker);

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (x, inexact, rnd_mode);
}
