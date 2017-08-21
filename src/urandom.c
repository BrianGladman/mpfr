/* mpfr_urandom (rop, state, rnd_mode) -- Generate a uniform pseudorandom
   real number between 0 and 1 (exclusive) and round it to the precision of rop
   according to the given rounding mode.

Copyright 2000-2004, 2006-2017 Free Software Foundation, Inc.
Contributed by the AriC and Caramba projects, INRIA.

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


/* generate one random bit */
static int
random_rounding_bit (gmp_randstate_t rstate)
{
  mp_limb_t r;

  mpfr_rand_raw (&r, rstate, 1);
  return r & MPFR_LIMB_ONE;
}

/* FIXME: Generate exceptions (flags) when need be, as usual? The current
   behavior is inconsistent as when reached, mpfr_check_range will detect
   some exceptions and set the corresponding flag. */

/* The mpfr_urandom() function is implemented in the following way for
   the generic case.
   1. One determines the exponent exp: 0 with probability 1/2, -1 with
      probability 1/4, -2 with probability 1/8, etc.
   2. One draws a 1-ulp interval ]a,b[ containing the exact result (the
      interval can be regarded as open since it has the same measure as
      the closed interval).
   3. Rounding is done. For the directed rounding modes, the rounded value
      is uniquely determined. For rounding to nearest, ]a,m[ and ]m,b[,
      where m = (a+b)/2, have the same measure, so that one gets a or b
      with equal probabilities. */

int
mpfr_urandom (mpfr_ptr rop, gmp_randstate_t rstate, mpfr_rnd_t rnd_mode)
{
  mpfr_limb_ptr rp;
  mpfr_prec_t nbits;
  mp_size_t nlimbs;
  mp_size_t n;
  mpfr_exp_t exp;
  int cnt;
  int inex;

  rp = MPFR_MANT (rop);
  nbits = MPFR_PREC (rop);
  MPFR_SET_POS (rop);

  if (MPFR_UNLIKELY (__gmpfr_emin > 0))
    {
      /* The minimum positive representable number 2^(emin-1) is >= 1,
         so that we need to round to +0 or 2^(emin-1). For the directed
         rounding modes, the rounded value is uniquely determined. For
         rounding to nearest: if emin = 1, one has probability 1/2 for
         each; otherwise (i.e. if emin > 1), the rounded value is 0. */
      if (rnd_mode == MPFR_RNDU || rnd_mode == MPFR_RNDA
          || (__gmpfr_emin == 1 && rnd_mode == MPFR_RNDN
              && random_rounding_bit (rstate)))
        {
          mpfr_set_ui_2exp (rop, 1, __gmpfr_emin - 1, rnd_mode);
          return +1;
        }
      else
        {
          MPFR_SET_ZERO (rop);
          return -1;
        }
    }

  exp = 0;
  MPFR_ASSERTD (exp >= __gmpfr_emin);

  /* Step 1 (exponent). */
#define DRAW_BITS 8 /* we draw DRAW_BITS at a time */
  cnt = DRAW_BITS;
  MPFR_STAT_STATIC_ASSERT (DRAW_BITS <= GMP_NUMB_BITS);
  while (cnt == DRAW_BITS)
    {
      /* generate DRAW_BITS in rp[0] */
      mpfr_rand_raw (rp, rstate, DRAW_BITS);
      if (MPFR_UNLIKELY (rp[0] == 0))
        cnt = DRAW_BITS;
      else
        {
          count_leading_zeros (cnt, rp[0]);
          cnt -= GMP_NUMB_BITS - DRAW_BITS;
        }
      exp -= cnt;  /* no integer overflow */

      if (MPFR_UNLIKELY (exp < __gmpfr_emin))
        {
          /* To get here, we have been drawing more than -emin zeros
             in a row, then return 0 or the smallest representable
             positive number.

             The rounding to nearest mode is subtle:
             If exp == emin - 1, the rounding bit is set, except
             if cnt == DRAW_BITS in which case the rounding bit is
             outside rp[0] and must be generated. */
          if (rnd_mode == MPFR_RNDU || rnd_mode == MPFR_RNDA
              || (rnd_mode == MPFR_RNDN && exp == __gmpfr_emin - 1
                  && (cnt != DRAW_BITS || random_rounding_bit (rstate))))
            {
              mpfr_set_ui_2exp (rop, 1, __gmpfr_emin - 1, rnd_mode);
              return +1;
            }
          else
            {
              MPFR_SET_ZERO (rop);
              return -1;
            }
        }
      MPFR_ASSERTD (exp >= __gmpfr_emin);
    }

  MPFR_ASSERTD (exp >= __gmpfr_emin);
  MPFR_EXP (rop) = exp; /* Warning: may be larger than emax */

  /* special code for nbits = 1 */
  if (nbits == 1)
    {
      rp[0] = MPFR_LIMB_HIGHBIT;
      goto rounding;
    }

  /* Step 2 (significand): we need generate only nbits-1 bits, since the
     most significant bit is 1. */
  mpfr_rand_raw (rp, rstate, nbits - 1);
  nlimbs = MPFR_LIMB_SIZE (rop);
  n = nlimbs * GMP_NUMB_BITS - nbits;
  if (MPFR_LIKELY (n != 0)) /* this will put the low bits to zero */
    mpn_lshift (rp, rp, nlimbs, n);
  rp[nlimbs - 1] |= MPFR_LIMB_HIGHBIT;

 rounding:
  /* Rounding */
  if (rnd_mode == MPFR_RNDU || rnd_mode == MPFR_RNDA
      || (rnd_mode == MPFR_RNDN && random_rounding_bit (rstate)))
    {
      if (MPFR_UNLIKELY (exp > __gmpfr_emax))
        mpfr_set_inf (rop, +1);  /* overflow */
      else
        mpfr_nextabove (rop);
      inex = +1;
      /* There is an overflow in the first case and possibly in the second
         case. If this occurs, the flag will be set by mpfr_check_range. */
    }
  else
    inex = -1;

  return mpfr_check_range (rop, inex, rnd_mode);
}
