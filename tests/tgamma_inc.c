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

#define TEST_FUNCTION mpfr_gamma_inc
#define TWO_ARGS
#define TEST_RANDOM_POS2 0 /* the 2nd argument is never negative */
#define TGENERIC_NOWARNING 1
#define TEST_RANDOM_EMAX 32
#define TEST_RANDOM_EMIN -32
#include "tgeneric.c"

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
      if (randlimb () & 1)
        mpfr_neg (a, a, MPFR_RNDN);
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

static void
specials (void)
{
  mpfr_t a, x;

  mpfr_init2 (a, 2);
  mpfr_init2 (x, 2);

  /* check gamma_inc(0,1) = 0.219383934395520 */
  mpfr_set_ui (a, 0, MPFR_RNDN);
  mpfr_set_ui (x, 1, MPFR_RNDN);
  mpfr_gamma_inc (a, a, x, MPFR_RNDN);
  if (mpfr_cmp_ui_2exp (a, 1, -2))
    {
      printf ("Error for gamma_inc(0,1)\n");
      printf ("expected 0.25\n");
      printf ("got      ");
      mpfr_out_str (stdout, 10, 0, a, MPFR_RNDN);
      printf ("\n");
      exit (1);
    }

  mpfr_set_prec (a, 1);
  mpfr_set_prec (x, 1);
  mpfr_set_ui_2exp (a, 1, 32, MPFR_RNDN);
  mpfr_set_ui_2exp (x, 1, -32, MPFR_RNDN);
  mpfr_gamma_inc (a, a, x, MPFR_RNDN);

  mpfr_clear (a);
  mpfr_clear (x);
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

  specials ();

  for (p = MPFR_PREC_MIN; p < 100; p++)
    test_random (p, 10);

  /* FIXME: once the case gamma_inc (-n, x) is implemented, we can activate
     the generic tests below */
  // test_generic (MPFR_PREC_MIN, 100, 100);

  tests_end_mpfr ();
  return 0;
}
