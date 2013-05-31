/* Declarations for mpfr_erandom and mpfr_nrandom.  Testing routine
   mpfr_urandom_alt is also declared.

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

#if !defined(MPFR_RANDOM_H)
#define MPFR_RANDOM_H 1

#include <mpfr.h>

#if defined(__cplusplus)
extern "C" {
#endif

  /* return an exponential random deviate with mean 1 as a MPFR */
  int mpfr_erandom(mpfr_t z, gmp_randstate_t r, mpfr_rnd_t rnd);

  /* return a normal random deviate with mean 0 and variance 1 as a MPFR */
  int mpfr_nrandom(mpfr_t z, gmp_randstate_t r, mpfr_rnd_t rnd);

#if MPFR_ALT_RANDOM
  /* mimic the behavior of mpfr_urandom */
  int mpfr_urandom_alt(mpfr_t z, gmp_randstate_t r, mpfr_rnd_t rnd);

  /* mimic the behavior of mpfr_grandom */
  int mpfr_grandom_alt(mpfr_t z1, mpfr_t z2,
                       gmp_randstate_t r, mpfr_rnd_t rnd);
#endif

#if defined(__cplusplus)
}
#endif

#endif
