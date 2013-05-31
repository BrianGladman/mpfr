/* mpfr_nrandom (rop, state, rnd_mode) -- Generate a normal deviate with mean 0
   and variance 1 and round it to the precision of rop according to the given
   rounding mode.

Copyright 2013 Free Software Foundation, Inc.
Contributed by Charles Karney <charles@karney.com>, SRI International.

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

/*
 * Sampling from the normal distribution with zero mean and unit variance.
 * This uses Algorithm N given in Charles F. F. Karney, Sampling exactly from
 * the normal distribution (March 2013), http://arxiv.org/abs/1303.6257v1 .
 *
 * The implementation here closely follows the C++ one given in the paper
 * above.  However, here, C is simplified by using gmp_urandomm_ui; the initial
 * rejection step in H just tests the leading bit of p; and the assignment of
 * the sign to the deviate using gmp_urandomb_ui.  Finally, the C++
 * implementation benefits from caching temporary random deviates across calls.
 * This isn't possible in C without additional machinery which would complicate
 * the interface.
 */

#include "random_deviate.h"

/* Algorithm H: true with probability exp(-1/2). */
static int
H (gmp_randstate_t r, mpfr_random_deviate_t p, mpfr_random_deviate_t q)
{
  /* p and q are temporaries */
  mpfr_random_deviate_reset (p);
  if (mpfr_random_deviate_tstbit (p, 1, r))
    return 1;
  for (;;)
    {
      mpfr_random_deviate_reset (q);
      if (!mpfr_random_deviate_less (q, p, r))
        return 0;
      mpfr_random_deviate_reset (p);
      if (!mpfr_random_deviate_less (p, q, r))
        return 1;
    }
}

/* Step N1: return n >= 0 with prob. exp(-n/2) * (1 - exp(-1/2)). */
static unsigned long
G (gmp_randstate_t r, mpfr_random_deviate_t p, mpfr_random_deviate_t q)
{
  /* p and q are temporaries */
  unsigned long n = 0;

  while (H (r, p, q))
    ++n;
  return n;
}

/* Step N2: true with probability exp(-n/2). */
static int
P (unsigned long n, gmp_randstate_t r,
   mpfr_random_deviate_t p, mpfr_random_deviate_t q)
{
  /* p and q are temporaries */
  while (n--)
    {
      if (!H (r, p, q))
        return 0;
    }
  return 1;
}

/* Algorithm C: return (-1, 0, 1) with prob (1/m, 1/m, 1-2/m). */
static int
C (unsigned long m, gmp_randstate_t r)
{
  unsigned long n =  gmp_urandomm_ui (r, m);
  return n == 0 ? -1 : (n == 1 ? 0 : 1);
}

/* Algorithm B: true with prob exp(-x * (2*k + x) / (2*k + 2)). */
static int
B (unsigned long k, mpfr_random_deviate_t x, gmp_randstate_t r,
   mpfr_random_deviate_t p, mpfr_random_deviate_t q)
{
  /* p and q are temporaries */
  unsigned long n = 0, m = 2 * k + 2;
  int  f;

  for (;; ++n)
    {
      if ( ((f = k ? 0 : C (m, r)) < 0) ||
           (mpfr_random_deviate_reset (q),
            !mpfr_random_deviate_less (q, n ? p : x, r)) ||
           ((f = k ? C (m, r) : f) < 0) ||
           (f == 0 &&
            (mpfr_random_deviate_reset (p),
             !mpfr_random_deviate_less (p, x, r))) )
        break;
    mpfr_random_deviate_swap (p, q); /* an efficient way of doing p = q */
    }
  return (n & 1UL) == 0;
}

/* return a normal random deviate with mean 0 and variance 1 as a MPFR  */
int
mpfr_nrandom (mpfr_t z, gmp_randstate_t r, mpfr_rnd_t rnd)
{
  mpfr_random_deviate_t x, p, q;
  int inex;
  unsigned long k, j;

  mpfr_random_deviate_init (x);
  mpfr_random_deviate_init (p);
  mpfr_random_deviate_init (q);
  for (;;)
    {
      k = G (r, p, q);                               /* step 1 */
      if (!P (k * (k - 1), r, p, q))
        continue;                                    /* step 2 */
      mpfr_random_deviate_reset (x);                 /* step 3 */
      for (j = 0; j <= k && B (k, x, r, p, q); ++j); /* step 4 */
      if (j > k)
        break;
    }
  mpfr_random_deviate_clear (q);
  mpfr_random_deviate_clear (p);
  /* steps 5, 6, 7 */
  inex = mpfr_random_deviate_value (gmp_urandomb_ui (r, 1), k, x, z, r, rnd);
  mpfr_random_deviate_clear (x);
  return inex;
}
