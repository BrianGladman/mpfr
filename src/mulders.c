/* Mulders' MulHigh function (short product)

Copyright 2005, 2006, 2007, 2008, 2009, 2010 Free Software Foundation, Inc.
Contributed by the Arenaire and Caramel projects, INRIA.

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

#ifndef MUL_FFT_THRESHOLD
#define MUL_FFT_THRESHOLD 8448
#endif

/* Don't use MPFR_MULHIGH_SIZE since it is handled by tuneup */
#ifdef MPFR_MULHIGH_TAB_SIZE
static short mulhigh_ktab[MPFR_MULHIGH_TAB_SIZE];
#else
static short mulhigh_ktab[] = {MPFR_MULHIGH_TAB};
#define MPFR_MULHIGH_TAB_SIZE \
  ((mp_size_t) (sizeof(mulhigh_ktab) / sizeof(mulhigh_ktab[0])))
#endif

/* Put in  rp[n..2n-1] an approximation of the n high limbs
   of {mp, n} * {np, n}.
   The error is at worst of ln(n) for rp[n] and rp[n-1] is totally wrong */
static void
mpfr_mulhigh_n_basecase (mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n)
{
  mp_size_t i;

  rp += n-1;
  umul_ppmm (rp[1], rp[0], up[n-1], vp[0]);
  for (i = 1 ; i < n ; i++)
    rp[i+1] = mpn_addmul_1 (rp, up + (n - i - 1), i+1, vp[i]);
}

void
mpfr_mulhigh_n (mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
  mp_size_t k;

  MPFR_ASSERTN (MPFR_MULHIGH_TAB_SIZE >= 8); /* so that 3*(n/4) > n/2 */
  k = MPFR_LIKELY (n < MPFR_MULHIGH_TAB_SIZE) ? mulhigh_ktab[n] : 3*(n/4);
  MPFR_ASSERTD (k == -1 || k == 0 || (k > n/2 && k < n));
  if (k < 0)
    mpn_mul_basecase (rp, np, n, mp, n);
  else if (k == 0)
    mpfr_mulhigh_n_basecase (rp, np, mp, n);
  else if (n > MUL_FFT_THRESHOLD)
    mpn_mul_n (rp, np, mp, n);
  else
    {
      mp_size_t l = n - k;
      mp_limb_t cy;

      mpn_mul_n (rp + 2 * l, np + l, mp + l, k); /* fills rp[2l..2n-1] */
      mpfr_mulhigh_n (rp, np + k, mp, l);          /* fills rp[l-1..2l-1] */
      cy = mpn_add_n (rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
      mpfr_mulhigh_n (rp, np, mp + k, l);          /* fills rp[l-1..2l-1] */
      cy += mpn_add_n (rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
      mpn_add_1 (rp + n + l, rp + n + l, k, cy); /* propagate carry */
    }
}

#ifdef MPFR_SQRHIGH_TAB_SIZE
static short sqrhigh_ktab[MPFR_SQRHIGH_TAB_SIZE];
#else
static short sqrhigh_ktab[] = {MPFR_SQRHIGH_TAB};
#define MPFR_SQRHIGH_TAB_SIZE (sizeof(sqrhigh_ktab) / sizeof(sqrhigh_ktab[0]))
#endif

void
mpfr_sqrhigh_n (mp_ptr rp, mp_srcptr np, mp_size_t n)
{
  mp_size_t k;

  MPFR_ASSERTN (MPFR_SQRHIGH_TAB_SIZE > 2); /* ensures k < n */
  k = MPFR_LIKELY (n < MPFR_SQRHIGH_TAB_SIZE) ? sqrhigh_ktab[n]
    : n/2 + n/50 + 1; /* +1 ensures that k > n/2 */
  MPFR_ASSERTD (k == -1 || k == 0 || (k > n/2 && k < n));
  if (k < 0)
    /* we can't use mpn_sqr_basecase here, since it requires
       n <= SQR_KARATSUBA_THRESHOLD, where SQR_KARATSUBA_THRESHOLD
       is not exported by GMP */
    mpn_sqr_n (rp, np, n);
  else if (k == 0)
    mpfr_mulhigh_n_basecase (rp, np, np, n);
  else
    {
      mp_size_t l = n - k;
      mp_limb_t cy;

      mpn_sqr_n (rp + 2 * l, np + l, k);          /* fills rp[2l..2n-1] */
      mpfr_mulhigh_n (rp, np, np + k, l);         /* fills rp[l-1..2l-1] */
      /* {rp+n-1,l+1} += 2 * {rp+l-1,l+1} */
      /* mpn_lshift is 30% to 50% faster than mpn_add_n on Core 2 */
#if 0
      cy = mpn_add_n (rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
#else
      cy = mpn_lshift (rp + l - 1, rp + l - 1, l + 1, 1);
#endif
      cy += mpn_add_n (rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
      mpn_add_1 (rp + n + l, rp + n + l, k, cy); /* propagate carry */
    }
}
