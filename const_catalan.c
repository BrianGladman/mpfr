/* mpfr_const_catalan -- compute Catalan's constant.

Copyright 2005 Free Software Foundation.

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

#include "mpfr-impl.h"

/* Declare the cache */
MPFR_DECL_INIT_CACHE(__gmpfr_cache_const_catalan, mpfr_const_catalan_internal);

/* Set User Interface */
#undef mpfr_const_catalan
int
mpfr_const_catalan (mpfr_ptr x, mp_rnd_t rnd_mode) {
  return mpfr_cache (x, __gmpfr_cache_const_catalan, rnd_mode);
}

/* return T, Q such that T/Q = sum(k!^2/(2k)!/(2k+1)^2, k=n1..n2-1) */
static void
S (mpz_t T, mpz_t P, mpz_t Q, unsigned long n1, unsigned long n2)
{
  if (n2 == n1 + 1)
    {
      if (n1 == 0)
        {
          mpz_set_ui (P, 1);
          mpz_set_ui (Q, 1);
        }
      else
        {
          mpz_set_ui (P, 2 * n1 - 1);
          mpz_mul_ui (P, P, n1);
          mpz_ui_pow_ui (Q, 2 * n1 + 1, 2);
          mpz_mul_2exp (Q, Q, 1);
        }
      mpz_set (T, P);
    }
  else
    {
      unsigned long m = (n1 + n2) / 2;
      mpz_t T2, P2, Q2;
      S (T, P, Q, n1, m);
      mpz_init (T2);
      mpz_init (P2);
      mpz_init (Q2);
      S (T2, P2, Q2, m, n2);
      mpz_mul (T, T, Q2);
      mpz_mul (T2, T2, P);
      mpz_add (T, T, T2);
      mpz_mul (P, P, P2);
      mpz_mul (Q, Q, Q2);
      mpz_clear (T2);
      mpz_clear (P2);
      mpz_clear (Q2);
    }
}

/* Don't need to save/restore exponent range: the cache does it.
   Catalan's constant is G = sum((-1)^k/(2*k+1)^2, k=0..infinity).
   We compute it using formula (31) of Victor Adamchik's page
   "33 representations for Catalan's constant"
   http://www-2.cs.cmu.edu/~adamchik/articles/catalan/catalan.htm

   G = Pi/8*log(2+sqrt(3)) + 3/8*sum(k!^2/(2k)!/(2k+1)^2,k=0..infinity)
*/
int
mpfr_const_catalan_internal (mpfr_ptr g, mp_rnd_t rnd_mode)
{
  mpfr_t x, y, z;
  mpz_t T, P, Q;
  mp_prec_t pg, p, t;
  MPFR_ZIV_DECL (loop);
  int inex;

  MPFR_LOG_FUNC (("rnd_mode=%d", rnd_mode), ("g[%#R]=%R inex=%d", g, g, inex));

  pg = MPFR_PREC (g);

  p = pg + 8; /* pg + 7 avoids failure up for pg < 912
		 pg + 8 gives no failure up to pg = 10000 */

  /* add about log2(p) bits */
  for (t = p; t; t >>= 1, p++);
  
  mpfr_init2 (x, p);
  mpfr_init2 (y, p);
  mpfr_init2 (z, p);
  mpz_init (T);
  mpz_init (P);
  mpz_init (Q);
  
  MPFR_ZIV_INIT (loop, p);
  for (;;) {    
    mpfr_sqrt_ui (x, 3, GMP_RNDU);
    mpfr_add_ui (x, x, 2, GMP_RNDU);
    mpfr_log (x, x, GMP_RNDU);
    mpfr_const_pi (y, GMP_RNDU);
    mpfr_mul (x, x, y, GMP_RNDN);
    S (T, P, Q, 0, (p - 1) / 2);
    mpz_mul_ui (T, T, 3);
    mpfr_set_z (y, T, GMP_RNDU);
    mpfr_set_z (z, Q, GMP_RNDD);
    mpfr_div (y, y, z, GMP_RNDN);
    mpfr_add (x, x, y, GMP_RNDN);
    mpfr_div_2exp (x, x, 3, GMP_RNDN);

    if (MPFR_LIKELY (MPFR_CAN_ROUND (x, p - 4, pg, rnd_mode)))
      break;

    MPFR_ZIV_NEXT (loop, p);
    mpfr_set_prec (x, p);
    mpfr_set_prec (y, p);
    mpfr_set_prec (z, p);
  }
  MPFR_ZIV_FREE (loop);
  inex = mpfr_set (g, x, rnd_mode);

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
  mpz_clear (T);
  mpz_clear (P);
  mpz_clear (Q);

  return inex;
}
