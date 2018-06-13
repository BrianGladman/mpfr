/* Test file for mpfr_set_decimal128 (and mpfr_get_decimal128 when available).

Copyright 2018 Free Software Foundation, Inc.
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
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

/* Needed due to the test on MPFR_WANT_DECIMAL_FLOATS */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#ifdef MPFR_WANT_DECIMAL_FLOATS

#include "mpfr-test.h"

static void
test_set (void)
{
  _Decimal128 d128;
  _Decimal64 d64;
  mpfr_t x;
  int inex;

  mpfr_init2 (x, 53);
  d64 = 1.0D;
  d128 = d64;
  inex = mpfr_set_decimal128 (x, d128, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui (x, 1) == 0);
  MPFR_ASSERTN(inex == 0);
  mpfr_clear (x);
}

int
main (void)
{
  tests_start_mpfr ();
  mpfr_test_init ();

#ifdef MPFR_DEBUG
#ifdef DPD_FORMAT
  printf ("Using DPD format\n");
#else
  printf ("Using BID format\n");
#endif
#endif

  test_set ();

  tests_end_mpfr ();
  return 0;
}

#else /* MPFR_WANT_DECIMAL_FLOATS */

int
main (void)
{
  return 77;
}

#endif /* MPFR_WANT_DECIMAL_FLOATS */

