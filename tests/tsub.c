/* Test file for mpfr_sub.

Copyright (C) 2001 Free Software Foundation.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Library General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
License for more details.

You should have received a copy of the GNU Library General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"
#include "mpfr-impl.h"
#include "mpfr-test.h"

void check_reuse _PROTO((void));

void
check_reuse ()
{
  mpfr_t x, y;

  mpfr_init2 (x, 5);
  mpfr_init2 (y, 5);
  mpfr_set_str_raw (x, "1e-12");
  mpfr_set_ui (y, 1, GMP_RNDN);
  mpfr_sub (x, y, x, GMP_RNDD);
  mpfr_set_str_raw (y, "0.11111");
  if (mpfr_cmp (x, y))
    {
      fprintf (stderr, "Error in mpfr_sub (x, y, x, GMP_RNDD) for x=2^(-12), y=1\n");
      exit (1);
    }
  mpfr_clear (x);
  mpfr_clear (y);
}

int
main()
{
  check_reuse ();

  return 0;
}
