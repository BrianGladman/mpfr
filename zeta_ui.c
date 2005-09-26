/* mpfr_zeta_ui -- compute the Riemann Zeta function for integer argument.

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
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA. */

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

int
mpfr_zeta_ui (mpfr_ptr z, unsigned long m, mp_rnd_t r)
{
  if (m == 0)
    {
      mpfr_set_ui (z, 1, r);
      mpfr_div_2ui (z, z, 1, r);
      MPFR_CHANGE_SIGN (z);
      MPFR_RET (0);
    }
  else if (m == 1)
    {
      MPFR_SET_INF (z);
      MPFR_SET_POS (z);
      return 0;
    }
  else /* m >= 2 */
    {
      mp_prec_t p = MPFR_PREC(z);
      unsigned long n, k, err, kbits;
      mpz_t d, t, s, q;
      mpfr_t y;
      int inex;
      int st = 0;

      if (m >= p) /* 2^(-m) < ulp(1) = 2^(1-p). This means that
		     2^(-m) <= 1/2*ulp(1). We have 3^(-m)+4^(-m)+... < 2^(-m)
		     i.e. zeta(m) < 1+2*2^(-m) for m >= 3 */
		     
	{
	  if (m == 2) /* necessarily p=2 */
	    return mpfr_set_ui_2exp (z, 13, -3, r);
	  else if (r == GMP_RNDZ || r == GMP_RNDD || (r == GMP_RNDN && m > p))
	    {
	      mpfr_set_ui (z, 1, r);
	      return -1;
	    }
	  else
	    {
	      mpfr_set_ui (z, 1, r);
	      mpfr_nextabove (z);
	      return 1;
	    }
	}

      /* FIXME: treat also the case where zeta(m) - (1+1/2^m) < ulp(1).
	 For m >= 6, we have zeta(m) - (1+1/2^m) < 2^(-1.5*m),
	 thus if p <= 1.5*m, then zeta(m) - (1+1/2^m) < 2^(-p)
	 and the result is either 1 or 1+. */

      mpz_init (s);
      mpz_init (d);
      mpz_init (t);
      mpz_init (q);
      mpfr_init2 (y, MPFR_PREC_MIN);

      p += MPFR_INT_CEIL_LOG2(p); /* account of the n term in the error */

      do
	{
	  p += MPFR_INT_CEIL_LOG2(p) + 15;
	  /* 0.39321985067869744 = log(2)/log(3+sqrt(8)) */
	  n = 1 + (unsigned long) (0.39321985067869744 * (double) p);
	  err = n + 4;

	  mpfr_set_prec (y, p);

	  /* computation of the d[k] */
	  mpz_set_ui (s, 0);
	  mpz_set_ui (t, 1);
	  mpz_mul_2exp (t, t, 2 * n - 1); /* t[n] */
	  mpz_set (d, t);
	  st = 0;
	  for (k = n; k > 0; k--)
	    {
	      count_leading_zeros (kbits, k);
	      kbits = m * (BITS_PER_MP_LIMB - kbits);
	      /* if k^m is too large, use mpz_tdiv_q */
	      if (kbits > 2 * BITS_PER_MP_LIMB)
		{
		  mpz_ui_pow_ui (q, k, m);
		  mpz_tdiv_q (q, d, q);
		}
	      else /* use several mpz_tdiv_q_ui calls */
		{
		  unsigned long km = k, mm = m - 1;
		  while (mm > 0 && km < ULONG_MAX / k)
		    {
		      km *= k;
		      mm --;
		    }
		  mpz_tdiv_q_ui (q, d, km);
		  while (mm > 0)
		    {
		      km = k;
		      mm --;
		      while (mm > 0 && km < ULONG_MAX / k)
			{
			  km *= k;
			  mm --;
			}
		      mpz_tdiv_q_ui (q, q, km);
		    }
		}
	      if (k % 2)
		mpz_add (s, s, q);
	      else
		mpz_sub (s, s, q);

	      /* we have d[k] = sum(t[i], i=k+1..n)
		 with t[i] = n*(n+i-1)!*4^i/(n-i)!/(2i)!
		 t[k-1]/t[k] = (2k)*(2k-1)/(n-k+1)/(n+k-1)/4 */
#if (BITS_PER_MP_LIMB == 32)
#define KMAX 46341 /* max k such that k*(2k-1) < 2^32 */
#elif (BITS_PER_MP_LIMB == 64)
#define KMAX 3037000500
#endif
#ifdef KMAX
	      if (k <= KMAX)
		mpz_mul_ui (t, t, k * (2 * k - 1));
	      else
#endif
		{
		  mpz_mul_ui (t, t, k);
		  mpz_mul_ui (t, t, 2 * k - 1);
		}
	      mpz_div_2exp (t, t, 1);
	      if (n < 1UL << (BITS_PER_MP_LIMB / 2))
		/* (n - k + 1) * (n + k - 1) < n^2 */
		mpz_divexact_ui (t, t, (n - k + 1) * (n + k - 1));
	      else
		{
		  mpz_divexact_ui (t, t, n - k + 1);
		  mpz_divexact_ui (t, t, n + k - 1);
		}
	      mpz_add (d, d, t);
      	    }

	  /* multiply by 1/(1-2^(1-m)) = 1 + 2^(1-m) + 2^(2-m) + ... */
	  mpz_div_2exp (t, s, m - 1);
	  do
	    {
	      err ++;
	      mpz_add (s, s, t);
	      mpz_div_2exp (t, t, m - 1);
	    }
	  while (mpz_cmp_ui (t, 0) > 0);

	  /* divide by d[n] */
	  mpz_mul_2exp (s, s, p);
	  mpz_tdiv_q (s, s, d);
	  mpfr_set_z (y, s, GMP_RNDN);
	  mpfr_div_2exp (y, y, p, GMP_RNDN);

	  err = MPFR_INT_CEIL_LOG2 (err);

	}
      while (MPFR_CAN_ROUND (y, p - err, MPFR_PREC(z), r) == 0);
      mpz_clear (d);
      mpz_clear (t);
      mpz_clear (q);
      mpz_clear (s);
      inex = mpfr_set (z, y, r);
      mpfr_clear (y);
      return inex;
    }
}
