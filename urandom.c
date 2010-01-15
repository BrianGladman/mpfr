/* mpfr_urandom (rop, state, rnd_mode) -- Generate a uniform pseudorandom
   real number between 0 and 1 (exclusive) and round it to the precision of rop
   according to the given rounding mode.

Copyright 2000, 2001, 2002, 2003, 2004, 2006, 2007, 2008, 2009 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

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

/* FIXME.
 * - The algorithm itself looks wrong (the bits should be determined MSB
 *   first in fixed point, so that there's no reason to determine the
 *   rounding bit first). Please include a proof.
 */

int
mpfr_urandom (mpfr_ptr rop, gmp_randstate_t rstate, mpfr_rnd_t rnd_mode)
{
  mp_ptr rp;
  mp_prec_t nbits;
  mp_size_t nlimbs;
  mp_size_t n;
  mp_size_t k; /* number of high zero limbs */
  mp_exp_t exp;
  int rndbit = 0;
  int cnt;
  int inex;

  rp = MPFR_MANT (rop);
  nbits = MPFR_PREC (rop);
  nlimbs = MPFR_LIMB_SIZE (rop);
  MPFR_SET_POS (rop);
  exp = 0;


  /* Rounding bit */
  if (rnd_mode == MPFR_RNDN)
    {
      mpfr_rand_raw (rp, rstate, GMP_NUMB_BITS);
      rndbit = rp[0] & 1;
    }


  /* Exponent */
  k = nlimbs;
  while (k == nlimbs)
    {
      mpfr_rand_raw (rp, rstate, nlimbs * GMP_NUMB_BITS);

      /* Count the null significant limbs k and remaining limbs n */
      k = 0;
      n = nlimbs;
      while (k != nlimbs && rp[n - 1] == 0)
        {
          k ++;
          n --;
        }

      if (exp < MPFR_EXP_MIN + k * GMP_NUMB_BITS)
        {
          /* To get here, we have been drawing more than 2^31 zeros in
             a raw (very unlucky). */
          MPFR_SET_ZERO (rop);
          if (rnd_mode == MPFR_RNDU || rnd_mode == MPFR_RNDA
              || (rnd_mode == MPFR_RNDN && rndbit))
            {
              MPFR_SAVE_EXPO_DECL (expo);
              MPFR_SAVE_EXPO_MARK (expo);
              mpfr_nextabove (rop);
              MPFR_SAVE_EXPO_FREE (expo);              
              return mpfr_check_range (rop, +1, rnd_mode);
            }
          return -1;
        }

      exp -= k * GMP_NUMB_BITS;
    }
  count_leading_zeros (cnt, rp[n - 1]);
  MPFR_EXP (rop) = exp - cnt;


  /* Significand */
  mpfr_rand_raw (rp, rstate, nlimbs * GMP_NUMB_BITS);

  /* Set the msb to 1 since it was fixed by the exponent choice */
  rp[nlimbs - 1] |= MPFR_LIMB_HIGHBIT;

  /* If nbits isn't a multiple of GMP_NUMB_BITS, mask the low bits */
  n = nlimbs * GMP_NUMB_BITS - nbits;
  if (MPFR_LIKELY (n != 0))
    rp[0] &= ~MPFR_LIMB_MASK (n);


  /* Rounding */
  if (rnd_mode == MPFR_RNDU || rnd_mode == MPFR_RNDA
      || (rnd_mode == MPFR_RNDN && rndbit))
    {
      MPFR_SAVE_EXPO_DECL (expo);
      MPFR_SAVE_EXPO_MARK (expo);
      mpfr_nextabove (rop);
      inex = +1;
      MPFR_SAVE_EXPO_FREE (expo);              
    }
  else
    inex = -1;

  return mpfr_check_range (rop, inex, rnd_mode);
}
