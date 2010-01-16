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


static int
random_rounding_bit (gmp_randstate_t rstate)
{
  mp_limb_t r;

  mpfr_rand_raw (&r, rstate, 1);
  return r & MPFR_LIMB_ONE;
}


int
mpfr_urandom (mpfr_ptr rop, gmp_randstate_t rstate, mpfr_rnd_t rnd_mode)
{
  mp_ptr rp;
  mp_prec_t nbits;
  mp_size_t nlimbs;
  mp_size_t n;
  mp_exp_t exp;
  mp_exp_t emin;
  int cnt;
  int inex;

  rp = MPFR_MANT (rop);
  nbits = MPFR_PREC (rop);
  nlimbs = MPFR_LIMB_SIZE (rop);
  MPFR_SET_POS (rop);
  exp = 0;
  emin = mpfr_get_emin ();


  /* Exponent */
  /* We first get a random limb rp[0] that cannot be zero, then the
     first non-zero bit determine the exponent. If only the very last
     bit is set, loop again. */
  /* FIXME: count_leading_zeros has an undefined behavior on 0.
     Also (probably related), GMP_NUMB_BITS - 1 seems incorrect. */
  cnt = GMP_NUMB_BITS - 1;
  while (cnt == GMP_NUMB_BITS - 1)
    {
      mpfr_rand_raw (rp, rstate, GMP_NUMB_BITS);
      count_leading_zeros (cnt, rp[0]);

      if (MPFR_UNLIKELY (exp < emin + cnt))
        {
          /* To get here, we have been drawing more than emin zeros in
             a raw, then return 0 or the smallest representable
             number.

             The rounding to nearest mode is subtle:
             If exp - cnt == emin - 1, the rounding bit is set except
             if cnt-th bit in the limb is the less significant bit. */
          MPFR_SET_ZERO (rop);
          if (rnd_mode == MPFR_RNDU || rnd_mode == MPFR_RNDA
              || (rnd_mode == MPFR_RNDN && cnt == exp - emin - 1
                  && (cnt!=GMP_NUMB_BITS-1 || random_rounding_bit (rstate))))
            {
              mpfr_set_ui_2exp (rop, 1, emin - 1, rnd_mode);
              return +1;
            }
          return -1;
        }
      exp -= cnt;
    }
  MPFR_EXP (rop) = exp; /* Warning: may be outside the current
                           exponent range */


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
      || (rnd_mode == MPFR_RNDN && random_rounding_bit (rstate)))
    {
      /* Take care of the exponent range: it may have been reduced */
      if (exp < emin)
        mpfr_set_ui_2exp (rop, 1, emin - 1, rnd_mode);
      else if (exp > mpfr_get_emax ())
        mpfr_set_inf (rop, +1); /* overflow, flag set by mpfr_check_range */
      else
        mpfr_nextabove (rop);
      inex = +1;
    }
  else
    inex = -1;

  return mpfr_check_range (rop, inex, rnd_mode);
}
