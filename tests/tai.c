/* Test file for mpfr_ai.

Copyright 2010 Free Software Foundation, Inc.
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

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "mpfr-test.h"


#define TEST_FUNCTION mpfr_ai
#define TEST_RANDOM_EMIN -5
#define TEST_RANDOM_EMAX 5
#include "tgeneric.c"

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  test_generic (2, 100, 100);

  tests_end_mpfr ();
  return 0;
}
