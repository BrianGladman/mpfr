/* Test file for mpfr_eint.

Copyright 2005, 2006 Free Software Foundation, Inc.

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

#define TEST_FUNCTION mpfr_eint
#include "tgeneric.c"

static void
check_specials (void)
{
  mpfr_t  x, y;

  mpfr_init2 (x, 123L);
  mpfr_init2 (y, 123L);

  mpfr_set_nan (x);
  mpfr_eint (y, x, GMP_RNDN);
  if (! mpfr_nan_p (y))
    {
      printf ("Error: eint(NaN) != NaN\n");
      exit (1);
    }

  mpfr_set_inf (x, 1);
  mpfr_eint (y, x, GMP_RNDN);
  if (! (mpfr_inf_p (y) && mpfr_sgn (y) > 0))
    {
      printf ("Error: eint(+Inf) != +Inf\n");
      exit (1);
    }

  mpfr_set_inf (x, -1);
  mpfr_eint (y, x, GMP_RNDN);
  if (! mpfr_nan_p (y))
    {
      printf ("Error: eint(-Inf) != NaN\n");
      exit (1);
    }

  /* eint(+/-0) = -Inf */
  mpfr_set_ui (x, 0, GMP_RNDN);
  mpfr_eint (y, x, GMP_RNDN);
  if (! (mpfr_inf_p (y) && mpfr_sgn (y) < 0))
    {
      printf ("Error: eint(+0) != -Inf\n");
      exit (1);
    }
  mpfr_neg (x, x, GMP_RNDN);
  mpfr_eint (y, x, GMP_RNDN);
  if (! (mpfr_inf_p (y) && mpfr_sgn (y) < 0))
    {
      printf ("Error: eint(-0) != -Inf\n");
      exit (1);
    }

  /* eint(x) = NaN for x < 0 */
  mpfr_set_si (x, -1, GMP_RNDN);
  mpfr_eint (y, x, GMP_RNDN);
  if (! mpfr_nan_p (y))
    {
      printf ("Error: eint(-1) != NaN\n");
      exit (1);
    }

  mpfr_set_prec (x, 17);
  mpfr_set_prec (y, 17);
  mpfr_set_str_binary (x, "1.0111110100100110e-2");
  mpfr_set_str_binary (y, "-1.0010101001110100e-10");
  mpfr_eint (x, x, GMP_RNDZ);
  if (mpfr_cmp (x, y))
    {
      printf ("Error for x=1.0111110100100110e-2, GMP_RNDZ\n");
      printf ("expected "); mpfr_dump (y);
      printf ("got      "); mpfr_dump (x);
      exit (1);
    }

  mpfr_set_prec (x, 53);
  mpfr_set_prec (y, 53);
  mpfr_set_str_binary (x, "0.10E4");
  mpfr_eint (x, x, GMP_RNDN);
  mpfr_set_str (y, "440.37989953483827", 10, GMP_RNDN);
  if (mpfr_cmp (x, y) != 0)
    {
      printf ("Error for x=0.10E4, GMP_RNDZ\n");
      printf ("expected "); mpfr_dump (y);
      printf ("got      "); mpfr_dump (x);
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);
}

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  if (argc != 1) /* teint x [prec] */
    {
      mpfr_t x;
      mp_prec_t p;
      p = (argc < 3) ? 53 : atoi (argv[2]);
      mpfr_init2 (x, p);
      mpfr_set_str (x, argv[1], 10, GMP_RNDN);
      printf ("eint(");
      mpfr_out_str (stdout, 10, 0, x, GMP_RNDN);
      printf (")=");
      mpfr_eint (x, x, GMP_RNDN);
      mpfr_out_str (stdout, 10, 0, x, GMP_RNDN);
      printf ("\n");
      mpfr_clear (x);
    }
  else
    {
      check_specials ();

      test_generic (2, 100, 100);
    }

  tests_end_mpfr ();
  return 0;
}
