/* mpfr_set_z -- set a floating-point number from a multiple-precision integer

Copyright 1999, 2000, 2001, 2002, 2003, 2004 Free Software Foundation, Inc.

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

#include <stdio.h>

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* set f to the integer z */
int 
mpfr_set_z (mpfr_ptr f, mpz_srcptr z, mp_rnd_t rnd_mode)
{
  mp_size_t fn, zn, dif;
  int k, sign_z, inex;
  mp_limb_t *fp, *zp;
  mp_exp_t exp;
  
  sign_z = mpz_sgn (z);
  if (MPFR_UNLIKELY (sign_z == 0))
    {
      MPFR_SET_ZERO(f);
      MPFR_SET_POS(f);
      MPFR_RET(0);
    }
  MPFR_ASSERTD (sign_z == MPFR_SIGN_POS || sign_z == MPFR_SIGN_NEG);  

  fp = MPFR_MANT(f);
  fn = MPFR_LIMB_SIZE(f);
  zn = ABS(SIZ(z));
  MPFR_ASSERTD (zn >= 1);
  dif = zn - fn;
  zp = PTR(z);
  count_leading_zeros (k, zp[zn-1]);

  if (MPFR_UNLIKELY (zn > MPFR_EMAX_MAX / BITS_PER_MP_LIMB + 1))
    return mpfr_set_overflow(f, rnd_mode, sign_z);
  /* because zn >= __gmpfr_emax / BITS_PER_MP_LIMB + 2
     and zn * BITS_PER_MP_LIMB >= __gmpfr_emax + BITS_PER_MP_LIMB + 1
     and exp = zn * BITS_PER_MP_LIMB - k > __gmpfr_emax */

  /* now zn <= MPFR_EMAX_MAX / BITS_PER_MP_LIMB + 1
     thus zn * BITS_PER_MP_LIMB <= MPFR_EMAX_MAX + BITS_PER_MP_LIMB
     and exp = zn * BITS_PER_MP_LIMB - k
             <= MPFR_EMAX_MAX + BITS_PER_MP_LIMB */
  exp = (mp_prec_t) zn * BITS_PER_MP_LIMB - k;
  /* The exponent will be exp or exp + 1 (due to rounding) */
  if (MPFR_UNLIKELY (exp > __gmpfr_emax))
    return mpfr_set_overflow (f, rnd_mode, sign_z);
  if (MPFR_UNLIKELY (exp + 1 < __gmpfr_emin))
    return mpfr_set_underflow(f, rnd_mode == GMP_RNDN ? GMP_RNDZ : rnd_mode,
                              sign_z);

  if (MPFR_LIKELY (dif >= 0))
    {
      mp_limb_t rb, sb, ulp;
      int sh;

      /* number has to be truncated */
      if (MPFR_LIKELY (k != 0))
        {
          mpn_lshift (fp, &zp[dif], fn, k);
          if (MPFR_LIKELY (dif > 0))
            fp[0] |= zp[dif - 1] >> (BITS_PER_MP_LIMB - k);
        }
      else
        MPN_COPY (fp, zp + dif, fn);
      
      /* Compute Rounding Bit and Sticky Bit */
      MPFR_UNSIGNED_MINUS_MODULO (sh, MPFR_PREC (f) );
      if (MPFR_LIKELY (sh != 0)) 
	{
	  mp_limb_t mask = MPFR_LIMB_ONE << (sh-1);
	  mp_limb_t limb = fp[0];
	  rb = limb & mask;
	  sb = limb & (mask-1);
	  ulp = 2*mask;
	  fp[0] = limb & ~(ulp-1);
	}
      else /* sh == 0 */
	{
	  mp_limb_t mask = MPFR_LIMB_ONE << (BITS_PER_MP_LIMB - 1 - k);
	  if (MPFR_LIKELY (dif > 0))
	    {
	      rb = zp[--dif] & mask;
	      sb = zp[dif] & (mask-1);
	    }
	  else
	    rb = sb = 0;
	  k = 0;
	  ulp = MPFR_LIMB_ONE;
	}
      if (MPFR_UNLIKELY (sb == 0) && MPFR_LIKELY (dif > 0))
	{
	  sb = zp[--dif];
	  if (MPFR_LIKELY (k != 0))
	    sb &= MPFR_LIMB_MASK (BITS_PER_MP_LIMB - k);
	  if (MPFR_UNLIKELY (sb == 0) && MPFR_LIKELY (dif > 0)) 
	    do {
	      sb = zp[--dif];
	    } while (dif > 0 && sb == 0);
	}

      /* Rounding */
      if (MPFR_LIKELY (rnd_mode == GMP_RNDN))
        {
	  if (rb == 0 || MPFR_UNLIKELY (sb == 0 && (fp[0] & ulp) == 0))
	    goto trunc;
	  else
	    goto addoneulp;
	}
      else /* Not Nearest */
	{
	  if (MPFR_LIKELY (MPFR_IS_LIKE_RNDZ (rnd_mode, sign_z < 0))
	      || MPFR_UNLIKELY ( (sb|rb) == 0 ))
	    goto trunc;
	  else
	    goto addoneulp;
	}

    trunc:
      inex = MPFR_LIKELY ((sb | rb) != 0) ? -1 : 0;
      goto end;

    addoneulp:
      inex = 1;
      if (MPFR_UNLIKELY (mpn_add_1 (fp, fp, fn, ulp))) 
	{ 
	  /* Pow 2 case */
	  if (MPFR_UNLIKELY (exp == __gmpfr_emax))
	    return mpfr_set_overflow (f, rnd_mode, sign_z);
	  exp ++;
	  fp[fn-1] = MPFR_LIMB_HIGHBIT;
	}
    end:
      (void) 0;
    }
  else   /* dif < 0: Mantissa F is strictly bigger than z's one */
    {
      if (MPFR_LIKELY (k != 0))
        mpn_lshift (fp - dif, zp, zn, k);
      else
        MPN_COPY (fp - dif, zp, zn);
      /* fill with zeroes */
      MPN_ZERO (fp, -dif);
      inex = 0; /* result is exact */
    }

  if (MPFR_UNLIKELY (exp < __gmpfr_emin))
    {
      if (rnd_mode == GMP_RNDN && inex == 0 && mpfr_powerof2_raw (f))
        rnd_mode = GMP_RNDZ;
      return mpfr_set_underflow(f, rnd_mode, sign_z);
    }

  MPFR_SET_EXP (f, exp);
  MPFR_SET_SIGN (f, sign_z);
  MPFR_RET (inex*sign_z);
}

