/* mpfr_erandom (rop, state, rnd_mode) -- Generate an exponential deviate with
   mean 1 and round it to the precision of rop according to the given rounding
   mode.

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
 * Sampling from the exponential distribution with unit mean using the method
 * given in John von Neumann, Various techniques used in connection with random
 * digits, in A. S. Householder, G. E. Forsythe, and H. H. Germond, editors,
 * "Monte Carlo Method", number 12 in Applied Mathematics Series, pp. 36-38
 * (NBS, Washington, DC, 1951), Proceedings of a symposium held June 29-July 1,
 * 1949, in Los Angeles.
 *
 * A modification to this algorithm is given in Charles F. F. Karney, Sampling
 * exactly from the normal distribution (March 2013),
 * http://arxiv.org/abs/1303.6257v1 .  Although this improves the bit
 * efficiency, in practice, it results in a slightly slower algorithm for MPFR.
 * So here the original von Neumann algorithm is used.
 */

#include "random_deviate.h"

/* true with prob exp(-x) */
static int
E (mpfr_random_deviate_t x, gmp_randstate_t r,
   mpfr_random_deviate_t p, mpfr_random_deviate_t q)
{
  /* p and q are temporaries */
  mpfr_random_deviate_reset (p);
  if (!mpfr_random_deviate_less (p, x, r))
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

/* return an exponential random deviate with mean 1 as a MPFR  */
int
mpfr_erandom (mpfr_t z, gmp_randstate_t r, mpfr_rnd_t rnd)
{
  mpfr_random_deviate_t x, p, q;
  int inex;
  unsigned long k = 0;

  mpfr_random_deviate_init (x);
  mpfr_random_deviate_init (p);
  mpfr_random_deviate_init (q);
  while (!E(x, r, p, q))
    {
      ++k;
      mpfr_random_deviate_reset (x);
    }
  mpfr_random_deviate_clear (q);
  mpfr_random_deviate_clear (p);
  inex = mpfr_random_deviate_value (0, k, x, z, r, rnd);
  mpfr_random_deviate_clear (x);
  return inex;
}
