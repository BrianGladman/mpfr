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

/* Declare the cache */
MPFR_DECL_INIT_CACHE(__gmpfr_cache_const_log2, mpfr_const_log2_internal);

/* Set User interface */
#undef mpfr_const_log2
int
mpfr_const_log2 (mpfr_ptr x, mp_rnd_t rnd_mode) {
  return mpfr_cache (x, __gmpfr_cache_const_log2, rnd_mode);
}

/* Auxiliary function: Compute the terms from n1 to n2 (excluded) 
   3/4*sum((-1)^n*n!^2/2^n/(2*n+1)!, n = n1..n2-1).
   Numerator is T, denominator is Q.
   Compute P only when need_P is non-zero.
*/
static void
S (mpz_t T, mpz_t P, mpz_t Q, unsigned long n1, unsigned long n2, int need_P)
{
  if (n2 == n1 + 1)
    {
      if (n1 == 0)
	mpz_set_ui (P, 3);
      else
	{
	  mpz_set_ui (P, n1);
	  mpz_neg (P, P);
	}
      if (n1 <= (ULONG_MAX / 4 - 1) / 2)
	mpz_set_ui (Q, 4 * (2 * n1 + 1));
      else /* to avoid overflow in 4 * (2 * n1 + 1) */
	{
	  mpz_set_ui (Q, n1);
	  mpz_mul_2exp (Q, Q, 1);
	  mpz_add_ui (Q, Q, 1);
	  mpz_mul_2exp (Q, Q, 2);
	}
      mpz_set (T, P);
    }
  else
    {
      unsigned long m = (n1 / 2) + (n2 / 2) + (n1 & 1UL & n2);
      mpz_t T2, P2, Q2;
      unsigned long v, w;
      S (T, P, Q, n1, m, 1);
      mpz_init (T2);
      mpz_init (P2);
      mpz_init (Q2);
      S (T2, P2, Q2, m, n2, need_P);
      mpz_mul (T, T, Q2);
      mpz_mul (T2, T2, P);
      mpz_add (T, T, T2);
      if (need_P)
	mpz_mul (P, P, P2);
      mpz_mul (Q, Q, Q2);

      /* remove common trailing zeroes if any */
      v = mpz_scan1 (T, 0);
      if (v > 0)
	{
	  w = mpz_scan1 (Q, 0);
	  if (w < v)
	    v = w;
	  if (need_P)
	    {
	      w = mpz_scan1 (P, 0);
	      if (w < v)
		v = w;
	    }
	  /* now v = min(val(T), val(Q), val(P)) */
	  if (v > 0)
	    {
	      mpz_div_2exp (T, T, v);
	      mpz_div_2exp (Q, Q, v);
	      if (need_P)
		mpz_div_2exp (P, P, v);
	    }
	}

      mpz_clear (T2);
      mpz_clear (P2);
      mpz_clear (Q2);
    }
}

int
mpfr_const_log2_internal (mpfr_ptr x, mp_rnd_t rnd_mode)
{
  unsigned long n = mpfr_get_prec (x);
  mp_prec_t w; /* working precision */
  unsigned long N;
  mpz_t T, P, Q;
  mpfr_t t, q;
  int inexact;
  int ok = 1; /* ensures that the 1st try will give correct rounding */
  MPFR_SAVE_EXPO_DECL (expo);

  MPFR_SAVE_EXPO_MARK (expo);
  mpz_init (T);
  mpz_init (P);
  mpz_init (Q);
  mpfr_init2 (t, MPFR_PREC_MIN);
  mpfr_init2 (q, MPFR_PREC_MIN);

  if (n < 1253)
    w = n + 10; /* ensures correct rounding for the four rounding modes,
		   together with N = w / 3 + 1 (see below). */
  else if (n < 2571)
    w = n + 11; /* idem */
  else if (n < 3983)
    w = n + 12;
  else if (n < 4854)
    w = n + 13;
  else if (n < 26248)
    w = n + 14;
  else
    {
      w = n + 15;
      ok = 0;
    }

  do
    {
      N = w / 3 + 1; /* Warning: do not change that (even increasing N!)
			without checking correct rounding in the above
			ranges for n. */

      /* the following are needed for error analysis (see algorithms.tex) */
      MPFR_ASSERTD(w >= 3 && N >= 2);

      S (T, P, Q, 0, N, 0);

      mpfr_set_prec (t, w);
      mpfr_set_prec (q, w);

      mpfr_set_z (t, T, GMP_RNDN);
      mpfr_set_z (q, Q, GMP_RNDN);
      mpfr_div (t, t, q, GMP_RNDN);

      if (ok == 0)
	{
	  ok = mpfr_can_round (t, w - 2, GMP_RNDN, rnd_mode, n);
	  if (ok == 0)
	    w += __gmpfr_ceil_log2 ((double) w);
	}
    }
  while (ok == 0);

  inexact = mpfr_set (x, t, rnd_mode);

  mpz_clear (T);
  mpz_clear (P);
  mpz_clear (Q);

  mpfr_clear (t);
  mpfr_clear (q);

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (x, inexact, rnd_mode);
}
