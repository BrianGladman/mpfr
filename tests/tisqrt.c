/* Test file for __gmpfr_isqrt internal function.

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

static void
tst_isqrt (unsigned long n, unsigned long r)
{
  unsigned long i;

  i = __gmpfr_isqrt (n);
  if (i != r)
    {
      printf ("Error in __gmpfr_isqrt (%lu): got %lu instead of %lu\n",
              n, i, r);
      exit (1);
    }
}

int
main (void)
{
  tst_isqrt (0, 0);
  tst_isqrt (1, 1);
  tst_isqrt (2, 1);
  tst_isqrt (3, 1);
  tst_isqrt (4, 2);
  tst_isqrt (8, 2);
  tst_isqrt (9, 3);
  tst_isqrt (15, 3);
  tst_isqrt (16, 4);
  tst_isqrt (255, 15);
  tst_isqrt (256, 16);
  tst_isqrt (288, 16);
  tst_isqrt (289, 17);
  tst_isqrt (4095, 63);
  tst_isqrt (4096, 64);
  tst_isqrt (65535, 255);
  tst_isqrt (65536, 256);
  tst_isqrt (262143, 511);
  tst_isqrt (262144, 512);
  tst_isqrt (16777215, 4095);
  tst_isqrt (16777216, 4096);
  tst_isqrt (4294967295UL, 65535);

  return 0;
}
