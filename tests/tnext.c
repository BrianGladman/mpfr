/* Test file for mpfr_nextabove, mpfr_nextbelow, mpfr_nexttoward.

Copyright 2003 Free Software Foundation.

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
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"
#include "mpfr-test.h"

/* Generic tests for mpfr_nextabove and mpfr_nextbelow */
static void
generic_abovebelow (void)
{
  int i;

  for (i = 0; i < 20000; i++)
    {
      mpfr_t x, y, z, t;
      mp_prec_t prec;
      int neg, below;

      prec = (randlimb () % 300) + MPFR_PREC_MIN;
      mpfr_inits2 (prec, x, y, z, (void *) 0);
      mpfr_init2 (t, 3);
      do
        mpfr_random (x);
      while (mpfr_cmp_ui (x, 0) == 0);
      neg = randlimb () & 1;
      if (neg)
        mpfr_neg (x, x, GMP_RNDN);
      mpfr_set (y, x, GMP_RNDN);
      below = randlimb () & 1;
      if (below)
        mpfr_nextbelow (y);
      else
        mpfr_nextabove (y);
      mpfr_set_si (t, below ? -5 : 5, GMP_RNDN);
      mpfr_mul_2si (t, t, mpfr_get_exp (x) - prec - 3, GMP_RNDN);
      /* t = (1/2 + 1/8) ulp(x) */
      mpfr_add (z, x, t, GMP_RNDN);
      if (!mpfr_number_p (y) || mpfr_cmp (y, z) != 0)
        {
          printf ("Error in mpfr_next%s for\n",
                  below ? "below" : "above");
          mpfr_out_str (stdout, 2, 0, x, GMP_RNDN);
          printf (", got\n");
          mpfr_out_str (stdout, 2, 0, y, GMP_RNDN);
          printf (" instead of\n");
          mpfr_out_str (stdout, 2, 0, z, GMP_RNDN);
          printf ("\n");
          exit (1);
        }
      mpfr_clears (x, y, z, t, (void *) 0);
    }
}

int
main (void)
{
  tests_start_mpfr ();
  generic_abovebelow ();
  tests_end_mpfr ();
  return 0;
}
