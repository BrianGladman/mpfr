/* tinternals -- Test for internals.

Copyright 2005 Free Software Foundation.
Contributed by the Spaces project, INRIA Lorraine.

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
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <stdio.h>

#define MPFR_NEED_LONGLONG_H
#include "mpfr-test.h"

int
main (void)
{
  int i;
  int val[16] = { 0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4 };

  tests_start_mpfr ();

  for (i = 1; i < 17; i++)
    {
      MPFR_ASSERTN (MPFR_INT_CEIL_LOG2 (i) == val[i-1]);
      MPFR_ASSERTN (MPFR_INT_CEIL_LOG2 (i) == __gmpfr_int_ceil_log2 (i));
    }

  tests_end_mpfr ();
  return 0;
}
