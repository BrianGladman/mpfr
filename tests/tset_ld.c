/* Test file for mpfr_set_ld and mpfr_get_ld.

Copyright 2002 Free Software Foundation, Inc.

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
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include "gmp.h"
#include "mpfr.h"
#include "mpfr-test.h"

int
main (int argc, char *argv[])
{
  /* Test removed as the minimal precision for a long double is something
     like 10 decimal digits. Anyway, there are implementations for which
     long double = double = IEEE double precision. */
#if 0
  long double d, e;
  mpfr_t x;

  tests_start_mpfr ();
  mpfr_test_init ();

  d = 9223372036854775808.0; /* 2^63 */
  e = d + 1.0;
  if (e == d)
    {
      fprintf (stderr, "long double should have >= 64 bits of mantissa\n");
      exit (1);
    }

  mpfr_init2 (x, 64);
  if (mpfr_set_ld (x, e, GMP_RNDN))
    {
      fprintf (stderr, "wrong inexact flag\n");
      exit (1);
    }
  /* check that x is 2^63+1 */
  if (mpfr_sub_ui (x, x, 1, GMP_RNDN))
    {
      fprintf (stderr, "wrong inexact flag\n");
      exit (1);
    }
  if (mpfr_cmp_ui_2exp (x, 1, 63))
    {
      fprintf (stderr, "Error: x should be 2^63\n");
      exit (1);
    }
  mpfr_clear (x);

  tests_end_mpfr ();
#endif

  return 0; 
}
