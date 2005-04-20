/* Test file for mpfr_const_catalan.

Copyright 2005 Free Software Foundation, Inc.

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

#include <stdlib.h>
#include "mpfr-test.h"

int
main (int argc, char *argv[])
{
  mpfr_t x;

  tests_start_mpfr ();

  mpfr_init2 (x, 32);
  mpfr_const_catalan (x, GMP_RNDN);
  mpfr_mul_2exp (x, x, 32, GMP_RNDN);
  if (mpfr_cmp_ui (x, 3934042271UL))
    {
      printf ("Error in const_catalan for prec=32\n");
      exit (1);
    }
  mpfr_clear (x);
  
  tests_end_mpfr ();
  return 0;
}
