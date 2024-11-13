/* test file for trigamma function

Copyright 2009-2024 Free Software Foundation, Inc.
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
along with the GNU MPFR Library; see the file COPYING.LESSER.
If not, see <https://www.gnu.org/licenses/>. */

#include "mpfr-test.h"

#define TEST_FUNCTION mpfr_trigamma
#include "tgeneric.c"

static void
special (void)
{
  mpfr_t x, y;

  mpfr_init (x);
  mpfr_init (y);

  mpfr_set_inf (y, -1);
  mpfr_set_inf (x, 1);
  mpfr_trigamma (y, x, MPFR_RNDN);
  mpfr_dump (y);
  if (mpfr_zero_p (y) == 0 || mpfr_signbit (y))
    {
      printf ("error for Trigamma(+Inf)\n");
      printf ("expected +0\n");
      printf ("got      ");
      mpfr_dump (y);
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);
}

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  special ();

  test_generic (MPFR_PREC_MIN, 200, 20);

  // data_check ("data/trigamma", mpfr_digamma, "mpfr_trigamma");

  tests_end_mpfr ();
  return 0;
}
