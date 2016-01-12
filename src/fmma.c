/* mpfr_fmma, mpfr_fmms -- Compute a*b +/- c*d

Copyright (C) 2014-2016 INRIA.
Contributed by the AriC and Caramel projects, INRIA.

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

#include "mpfr-impl.h"

static int
mpfr_fmma_slow (mpfr_ptr z, mpfr_srcptr a, mpfr_srcptr b, mpfr_srcptr c,
                   mpfr_srcptr d, mpfr_rnd_t rnd)
{
  mpfr_t ab, cd;
  int inex;

  mpfr_init2 (ab, MPFR_PREC(a) + MPFR_PREC(b));
  mpfr_init2 (cd, MPFR_PREC(c) + MPFR_PREC(d));
  mpfr_mul (ab, a, b, MPFR_RNDZ); /* exact */
  mpfr_mul (cd, c, d, MPFR_RNDZ); /* exact */
  inex = mpfr_add (z, ab, cd, rnd);
  mpfr_clear (ab);
  mpfr_clear (cd);
  return inex;
}

/* z <- a*b + c*d */
static int
mpfr_fmma_fast (mpfr_ptr z, mpfr_srcptr a, mpfr_srcptr b, mpfr_srcptr c,
                mpfr_srcptr d, mpfr_rnd_t rnd)
{
   /* Assumes that a, b, c, d are finite and non-zero; so any multiplication
      of two of them yielding an infinity is an overflow, and a
      multiplication yielding 0 is an underflow.
      Assumes further that z is distinct from a, b, c, d. */

   int inex;
   mpfr_t u, v;
   mp_size_t an, bn, cn, dn;
   mp_ptr up, vp;
   MPFR_TMP_DECL(marker);
   MPFR_SAVE_EXPO_DECL (expo);

   MPFR_TMP_MARK(marker);
   MPFR_SAVE_EXPO_MARK (expo);

   /* u=a*b, v=c*d exactly */
   an = MPFR_LIMB_SIZE(a);
   bn = MPFR_LIMB_SIZE(b);
   cn = MPFR_LIMB_SIZE(c);
   dn = MPFR_LIMB_SIZE(d);
   MPFR_TMP_INIT (up, u, (an + bn) * GMP_NUMB_BITS, an + bn);
   MPFR_TMP_INIT (vp, v, (cn + dn) * GMP_NUMB_BITS, cn + dn);

   /* u <- a*b */
   if (an >= bn)
     mpn_mul (up, MPFR_MANT(a), an, MPFR_MANT(b), bn);
   else
     mpn_mul (up, MPFR_MANT(b), bn, MPFR_MANT(a), an);
   if ((up[an + bn - 1] & MPFR_LIMB_HIGHBIT) == 0)
     {
       mpn_lshift (up, up, an + bn, 1);
       MPFR_EXP(u) = MPFR_EXP(a) + MPFR_EXP(b) - 1;
     }
   else
     MPFR_EXP(u) = MPFR_EXP(a) + MPFR_EXP(b);

   /* v <- c*d */
   if (cn >= dn)
     mpn_mul (vp, MPFR_MANT(c), cn, MPFR_MANT(d), dn);
   else
     mpn_mul (vp, MPFR_MANT(d), dn, MPFR_MANT(c), cn);
   if ((vp[cn + dn - 1] & MPFR_LIMB_HIGHBIT) == 0)
     {
       mpn_lshift (vp, vp, cn + dn, 1);
       MPFR_EXP(v) = MPFR_EXP(c) + MPFR_EXP(d) - 1;
     }
   else
     MPFR_EXP(v) = MPFR_EXP(c) + MPFR_EXP(d);

   MPFR_PREC(u) = (an + bn) * GMP_NUMB_BITS;
   MPFR_PREC(v) = (cn + dn) * GMP_NUMB_BITS;
   MPFR_SIGN(u) = MPFR_MULT_SIGN(MPFR_SIGN(a), MPFR_SIGN(b));
   MPFR_SIGN(v) = MPFR_MULT_SIGN(MPFR_SIGN(c), MPFR_SIGN(d));

   /* tentatively compute z as u+v; here we need z to be distinct
      from a, b, c, d to avoid losing the input values in case we
      need to call mpfr_fmma_slow */
   inex = mpfr_add (z, u, v, rnd);

   MPFR_TMP_FREE(marker);
   MPFR_SAVE_EXPO_FREE (expo);

   return mpfr_check_range (z, inex, rnd);
}

/* z <- a*b + c*d */
int
mpfr_fmma (mpfr_ptr z, mpfr_srcptr a, mpfr_srcptr b, mpfr_srcptr c,
           mpfr_srcptr d, mpfr_rnd_t rnd)
{
  mp_ptr zp = MPFR_MANT(z);

  return (mpfr_regular_p (a) && mpfr_regular_p (b) && mpfr_regular_p (c) &&
          mpfr_regular_p (d) && zp != MPFR_MANT(a) && zp != MPFR_MANT(b) &&
          zp != MPFR_MANT(c) && zp != MPFR_MANT(d))
    ? mpfr_fmma_fast (z, a, b, c, d, rnd)
    : mpfr_fmma_slow (z, a, b, c, d, rnd);
}

/* z <- a*b - c*d */
int
mpfr_fmms (mpfr_ptr z, mpfr_srcptr a, mpfr_srcptr b, mpfr_srcptr c,
           mpfr_srcptr d, mpfr_rnd_t rnd)
{
  mpfr_t minus_c;

  MPFR_ALIAS (minus_c, c, -MPFR_SIGN(c), MPFR_EXP(c));
  return mpfr_fmma (z, a, b, minus_c, d, rnd);
}
