/* Example illustrating how to use mpfr_can_round. */

/*
Copyright 2016 Free Software Foundation, Inc.
Contributed by the AriC and Caramel projects, INRIA.

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
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <stdio.h>
#include <mpfr.h>

int
main (void)
{
  mpfr_t x;
  mpfr_rnd_t r1, r2;
  int ok;

  /* Given an approximation of Pi to 53 bits computed with rounding mode r1,
     we call mpfr_can_round() to see if we can deduced the correct rounding
     of Pi to 50 bits with rounding mode r2.
     The error is at most 1 ulp for directed rounding, and 1/2 ulp for rounding
     to nearest. This translates into err = prec(x) for directed rounding,
     and err = prec(x) + 1 for rounding to nearest. */
  mpfr_init2 (x, 53);
  for (r1 = 0; r1 < 4; r1++)
    {
      mpfr_const_pi (x, r1);
      mpfr_printf ("rnd1=%s approx=%Rb\n", mpfr_print_rnd_mode (r1), x);
      for (r2 = 0; r2 < 4; r2++)
        {
          ok = mpfr_can_round (x, mpfr_get_prec (x) + (r1 == MPFR_RNDN),
                               r1, r2, 50);
          printf ("rnd2=%s can_round=%d\n", mpfr_print_rnd_mode (r2), ok);
        }
    }
  mpfr_clear (x);
  return 0;
}
