/* Declarations of random_deviate routines for mpfr_erandom and mpfr_nrandom.

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

#if !defined(MPFR_RANDOM_DEVIATE_H)
#define MPFR_RANDOM_DEVIATE_H 1

#include "mpfr-impl.h"

#if defined(__cplusplus)
extern "C" {
#endif

  /* This should be an unsigned type with a width of at least 32 and capable of
   * representing at least 2*MPFR_PREC_MAX.  This is used to count the bits in
   * the fraction of a mpfr_random_deviate_t.  See the checks made on this type
   * in random_deviate_generate. */
  typedef unsigned long mpfr_random_size_t;

  typedef struct {
    mpfr_random_size_t e;       /* total bits in the fraction */
    unsigned long h;            /* the high W bits of the fraction */
    mpz_t f;                    /* the rest of the fraction */
  } mpfr_random_deviate_t[1];

  /* allocate and set to (0,1) */
  void mpfr_random_deviate_init(mpfr_random_deviate_t x);

  /* reset to (0,1) */
  void mpfr_random_deviate_reset(mpfr_random_deviate_t x);

  /* deallocate */
  void mpfr_random_deviate_clear(mpfr_random_deviate_t x);

  /* swap two random deviates */
  void mpfr_random_deviate_swap(mpfr_random_deviate_t x,
                                mpfr_random_deviate_t y);

  /* return kth bit of fraction, representing 2^-k */
  int mpfr_random_deviate_tstbit(mpfr_random_deviate_t x, mpfr_random_size_t k,
                                 gmp_randstate_t r);

  /* compare two random deviates, x < y */
  int mpfr_random_deviate_less(mpfr_random_deviate_t x,
                               mpfr_random_deviate_t y,
                               gmp_randstate_t r);

  /* set mpfr_t z = (neg ? -1 : 1) * (n + x) */
  int mpfr_random_deviate_value(int neg, unsigned long n,
                                mpfr_random_deviate_t x, mpfr_t z,
                                gmp_randstate_t r, mpfr_rnd_t rnd);

#if defined(__cplusplus)
}
#endif

#endif
