/* Test file for mpfr_get_f.

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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mpfr-test.h"

int
main (void)
{
  mpf_t x;
  mpfr_t y;
  unsigned long i;

  MPFR_TEST_USE_RANDS ();
  tests_start_mpfr ();

  mpfr_init (y);
  mpf_init (x);

  mpfr_set_nan (y);
  if (mpfr_get_f (x, y, GMP_RNDN) == 0)
    {
      printf ("Error: mpfr_get_f(NaN) should fail\n");
      exit (1);
    }

  mpfr_set_inf (y, 1);
  if (mpfr_get_f (x, y, GMP_RNDN) == 0)
    {
      printf ("Error: mpfr_get_f(+Inf) should fail\n");
      exit (1);
    }

  mpfr_set_inf (y, -1);
  if (mpfr_get_f (x, y, GMP_RNDN) == 0)
    {
      printf ("Error: mpfr_get_f(-Inf) should fail\n");
      exit (1);
    }

  mpfr_set_ui (y, 0, GMP_RNDN);
  if (mpfr_get_f (x, y, GMP_RNDN) || mpf_cmp_ui (x, 0))
    {
      printf ("Error: mpfr_get_f(+0) fails\n");
      exit (1);
    }

  mpfr_set_ui (y, 0, GMP_RNDN);
  mpfr_neg (y, y, GMP_RNDN);
  if (mpfr_get_f (x, y, GMP_RNDN) || mpf_cmp_ui (x, 0))
    {
      printf ("Error: mpfr_get_f(-0) fails\n");
      exit (1);
    }

  i = 1;
  while (i)
    {
      mpfr_set_ui (y, i, GMP_RNDN);
      if (mpfr_get_f (x, y, GMP_RNDN) || mpf_cmp_ui (x, i))
        {
          printf ("Error: mpfr_get_f(%lu) fails\n", i);
          exit (1);
        }
      mpfr_set_si (y, (signed long) -i, GMP_RNDN);
      if (mpfr_get_f (x, y, GMP_RNDN) || mpf_cmp_si (x, (signed long) -i))
        {
          printf ("Error: mpfr_get_f(-%lu) fails\n", i);
          exit (1);
        }
      i *= 2;
    }

  /* same tests, but with a larger precision for y, which requires to
     round it */
  mpfr_set_prec (y, 100);
  i = 1;
  while (i)
    {
      mpfr_set_ui (y, i, GMP_RNDN);
      if (mpfr_get_f (x, y, GMP_RNDN) || mpf_cmp_ui (x, i))
        {
          printf ("Error: mpfr_get_f(%lu) fails\n", i);
          exit (1);
        }
      mpfr_set_si (y, (signed long) -i, GMP_RNDN);
      if (mpfr_get_f (x, y, GMP_RNDN) || mpf_cmp_si (x, (signed long) -i))
        {
          printf ("Error: mpfr_get_f(-%lu) fails\n", i);
          exit (1);
        }
      i *= 2;
    }

  mpfr_clear (y);
  mpf_clear (x);

  tests_end_mpfr ();
  return 0;
}
