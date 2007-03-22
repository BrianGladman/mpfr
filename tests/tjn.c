/* tjn -- test file for the Bessel function of first kind

Copyright 2007 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

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
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>

#include "mpfr-test.h"

int
main (int argc, char *argv[])
{
  mpfr_t x, y;

  tests_start_mpfr ();

  mpfr_init (x);
  mpfr_init (y);

  /* special values */
  mpfr_set_nan (x);
  mpfr_jn (y, x, 17, GMP_RNDN);
  MPFR_ASSERTN(mpfr_nan_p (y));

  mpfr_set_inf (x, 1); /* +Inf */
  mpfr_jn (y, x, 17, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (y, 0) == 0 && MPFR_IS_POS (y));

  mpfr_set_inf (x, -1); /* -Inf */
  mpfr_jn (y, x, 17, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (y, 0) == 0 && MPFR_IS_POS (y));

  mpfr_set_ui (x, 0, GMP_RNDN); /* +0 */
  mpfr_jn (y, x, 0, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (y, 1) == 0); /* j0(+0)=1 */
  mpfr_jn (y, x, 17, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (y, 0) == 0 && MPFR_IS_POS (y)); /* j17(+0)=+0 */
  mpfr_jn (y, x, -17, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (y, 0) == 0 && MPFR_IS_NEG (y)); /* j-17(+0)=-0 */
  mpfr_jn (y, x, 42, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (y, 0) == 0 && MPFR_IS_POS (y)); /* j42(+0)=+0 */

  mpfr_set_ui (x, 0, GMP_RNDN);
  mpfr_neg (x, x, GMP_RNDN); /* -0 */
  mpfr_jn (y, x, 0, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (y, 1) == 0); /* j0(-0)=1 */
  mpfr_jn (y, x, 17, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (y, 0) == 0 && MPFR_IS_NEG (y)); /* j17(-0)=-0 */
  mpfr_jn (y, x, -17, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (y, 0) == 0 && MPFR_IS_POS (y)); /* j-17(-0)=+0 */
  mpfr_jn (y, x, 42, GMP_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (y, 0) == 0 && MPFR_IS_POS (y)); /* j42(-0)=+0 */

  mpfr_set_prec (x, 53);
  mpfr_set_prec (y, 53);

  mpfr_set_ui (x, 1, GMP_RNDN);
  mpfr_jn (y, x, 0, GMP_RNDN);
  mpfr_set_str_binary (x, "0.1100001111100011111111101101111010111101110001111");
  if (mpfr_cmp (x, y))
    {
      printf ("Error in mpfr_jn for n=0, x=1, rnd=GMP_RNDN\n");
      printf ("Expected "); mpfr_dump (x);
      printf ("Got      "); mpfr_dump (y);
      exit (1);
    }

  mpfr_set_si (x, -1, GMP_RNDN);
  mpfr_jn (y, x, 0, GMP_RNDN);
  mpfr_set_str_binary (x, "0.1100001111100011111111101101111010111101110001111");
  if (mpfr_cmp (x, y))
    {
      printf ("Error in mpfr_jn for n=0, x=-1, rnd=GMP_RNDN\n");
      printf ("Expected "); mpfr_dump (x);
      printf ("Got      "); mpfr_dump (y);
      exit (1);
    }

  mpfr_set_ui (x, 1, GMP_RNDN);
  mpfr_jn (y, x, 1, GMP_RNDN);
  mpfr_set_str_binary (x, "0.0111000010100111001001111011101001011100001100011011");
  if (mpfr_cmp (x, y))
    {
      printf ("Error in mpfr_jn for n=1, x=1, rnd=GMP_RNDN\n");
      printf ("Expected "); mpfr_dump (x);
      printf ("Got      "); mpfr_dump (y);
      exit (1);
    }

  mpfr_set_ui (x, 1, GMP_RNDN);
  mpfr_jn (y, x, 17, GMP_RNDN);
  mpfr_set_str_binary (x, "0.1100011111001010101001001001000110110000010001011E-65");
  if (mpfr_cmp (x, y))
    {
      printf ("Error in mpfr_jn for n=17, x=1, rnd=GMP_RNDN\n");
      printf ("Expected "); mpfr_dump (x);
      printf ("Got      "); mpfr_dump (y);
      exit (1);
    }

  mpfr_set_ui (x, 1, GMP_RNDN);
  mpfr_jn (y, x, 42, GMP_RNDN);
  mpfr_set_str_binary (x, "0.10000111100011010100111011100111101101000100000001001E-211");
  if (mpfr_cmp (x, y))
    {
      printf ("Error in mpfr_jn for n=42, x=1, rnd=GMP_RNDN\n");
      printf ("Expected "); mpfr_dump (x);
      printf ("Got      "); mpfr_dump (y);
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);

  tests_end_mpfr ();

  return 0;
}
