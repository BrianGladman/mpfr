/* Test file for mpfr_min_prec.

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
along with the GNU MPFR Library; see the file COPYING.LIB.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>

#include "mpfr-test.h"

int
main (int argc, char *argv[])
{
  mpfr_t x;
  mpfr_prec_t ret;

  tests_start_mpfr ();

  mpfr_init (x);

  /* Check special values */
  mpfr_set_nan (x);
  ret = mpfr_min_prec (x);
  MPFR_ASSERTN(ret == 0);

  mpfr_set_inf (x, 1);
  ret = mpfr_min_prec (x);
  MPFR_ASSERTN(ret == 0);

  mpfr_set_inf (x, -1);
  ret = mpfr_min_prec (x);
  MPFR_ASSERTN(ret == 0);

  MPFR_ASSERTN(mpfr_set_ui (x, 0, MPFR_RNDN) == 0);
  ret = mpfr_min_prec (x);
  MPFR_ASSERTN(ret == 0);

  /* Some constants */

  MPFR_ASSERTN(mpfr_set_ui (x, 1, MPFR_RNDN) == 0);
  ret = mpfr_min_prec (x);
  MPFR_ASSERTN(ret == 1);

  MPFR_ASSERTN(mpfr_set_ui (x, 17, MPFR_RNDN) == 0);
  ret = mpfr_min_prec (x);
  MPFR_ASSERTN(ret == 5);

  MPFR_ASSERTN(mpfr_set_ui (x, 42, MPFR_RNDN) == 0);
  ret = mpfr_min_prec (x);
  MPFR_ASSERTN(ret == 5);

  mpfr_clear (x);

  tests_end_mpfr ();
  return 0;
}
