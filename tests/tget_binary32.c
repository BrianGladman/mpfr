/* Test file for mpfr_get_binary32 and mpfr_set_binary32

Copyright 2009 Free Software Foundation, Inc.
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

#include <stdlib.h>
#include "mpfr-test.h"

int
main (void)
{
  mpfr_t x;
  float f;

  tests_start_mpfr ();

  mpfr_init2 (x, 24);

  mpfr_set_nan (x);
  f = mpfr_get_binary32 (x, MPFR_RNDN);
  if (f == f)
    {
      printf ("Error for mpfr_get_binary32(NaN)\n");
      exit (1);
    }
  mpfr_set_binary32 (x, f, MPFR_RNDN);
  if (mpfr_nan_p (x) == 0)
    {
      printf ("Error for mpfr_set_binary32(NaN)\n");
      exit (1);
    }
  
  mpfr_set_inf (x, 1);
  f = mpfr_get_binary32 (x, MPFR_RNDN);
  mpfr_set_binary32 (x, f, MPFR_RNDN);
  if (mpfr_inf_p (x) == 0 || mpfr_sgn (x) < 0)
    {
      printf ("Error for mpfr_set_binary32(mpfr_get_binary32(+Inf))n");
      exit (1);
    }
  
  mpfr_set_inf (x, -1);
  f = mpfr_get_binary32 (x, MPFR_RNDN);
  mpfr_set_binary32 (x, f, MPFR_RNDN);
  if (mpfr_inf_p (x) == 0 || mpfr_sgn (x) > 0)
    {
      printf ("Error for mpfr_set_binary32(mpfr_get_binary32(-Inf))n");
      exit (1);
    }
  
  mpfr_set_ui (x, 0, MPFR_RNDN);
  f = mpfr_get_binary32 (x, MPFR_RNDN);
  mpfr_set_binary32 (x, f, MPFR_RNDN);
  if (mpfr_zero_p (x) == 0 || MPFR_SIGN (x) < 0)
    {
      printf ("Error for mpfr_set_binary32(mpfr_get_binary32(+0))n");
      exit (1);
    }
  
  mpfr_set_ui (x, 0, MPFR_RNDN);
  mpfr_neg (x, x, MPFR_RNDN);
  f = mpfr_get_binary32 (x, MPFR_RNDN);
  mpfr_set_binary32 (x, f, MPFR_RNDN);
  if (mpfr_zero_p (x) == 0 || MPFR_SIGN (x) > 0)
    {
      printf ("Error for mpfr_set_binary32(mpfr_get_binary32(-0))n");
      exit (1);
    }
  
  mpfr_set_ui (x, 17, MPFR_RNDN);
  f = mpfr_get_binary32 (x, MPFR_RNDN);
  mpfr_set_binary32 (x, f, MPFR_RNDN);
  if (mpfr_cmp_ui (x, 17) != 0)
    {
      printf ("Error for mpfr_set_binary32(mpfr_get_binary32(17))n");
      printf ("expected 17\n");
      printf ("got      ");
      mpfr_dump (x);
      exit (1);
    }
  
  mpfr_set_si (x, -42, MPFR_RNDN);
  f = mpfr_get_binary32 (x, MPFR_RNDN);
  mpfr_set_binary32 (x, f, MPFR_RNDN);
  if (mpfr_cmp_si (x, -42) != 0)
    {
      printf ("Error for mpfr_set_binary32(mpfr_get_binary32(-42))n");
      printf ("expected -42\n");
      printf ("got      ");
      mpfr_dump (x);
      exit (1);
    }
  
  mpfr_clear (x);
  
  tests_end_mpfr ();
  return 0;
}
