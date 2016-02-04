/* Test file for mpfr_gamma_inc

Copyright 2016 Free Software Foundation, Inc.
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

#include "mpfr-test.h"

/* do k random tests at precision p */
static void
test_random (mpfr_prec_t p, int k)
{
  mpfr_t a, x, y, z, t;

  mpfr_inits2 (p, a, x, y, z, (mpfr_ptr) 0);
  mpfr_init2 (t, p + 20);
  while (k--)
    {
      do mpfr_urandomb (a, RANDS); while (mpfr_zero_p (a));
      do mpfr_urandomb (x, RANDS); while (mpfr_zero_p (x));
      mpfr_gamma_inc (y, a, x, MPFR_RNDN);
      mpfr_gamma_inc (t, a, x, MPFR_RNDN);
      if (mpfr_can_round (t, mpfr_get_prec (z), MPFR_RNDN, MPFR_RNDN, p))
        {
          mpfr_set (z, t, MPFR_RNDN);
          if (mpfr_cmp (y, z))
            {
              printf ("mpfr_gamma_inc failed for a=");
              mpfr_out_str (stdout, 10, 0, a, MPFR_RNDN);
              printf (" x=");
              mpfr_out_str (stdout, 10, 0, x, MPFR_RNDN);
              printf ("\nexpected ");
              mpfr_out_str (stdout, 10, 0, z, MPFR_RNDN);
              printf ("\ngot      ");
              mpfr_out_str (stdout, 10, 0, y, MPFR_RNDN);
              printf ("\n");
              exit (1);
            }
        }
    }
  mpfr_clears (a, x, y, z, (mpfr_ptr) 0);
  mpfr_clear (t);
}

int
main (int argc, char *argv[])
{
  mpfr_prec_t p;

  tests_start_mpfr ();

  if (argc == 4) /* tgamma_inc a x prec: print gamma_inc(a,x) to prec bits */
    {
      mpfr_prec_t p = atoi (argv[3]);
      mpfr_t a, x;
      mpfr_init2 (a, p);
      mpfr_init2 (x, p);
      mpfr_set_str (a, argv[1], 10, MPFR_RNDN);
      mpfr_set_str (x, argv[2], 10, MPFR_RNDN);
      mpfr_gamma_inc (x, a, x, MPFR_RNDN);
      mpfr_out_str (stdout, 10, 0, x, MPFR_RNDN);
      printf ("\n");
      mpfr_clear (a);
      mpfr_clear (x);
      return 0;
    }

  for (p = MPFR_PREC_MIN; p < 100; p++)
    test_random (p, 10);

  tests_end_mpfr ();
  return 0;
}
