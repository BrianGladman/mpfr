/* Test file for mpfr_mul_2exp.

Copyright 1999, 2001, 2002, 2003, 2004 Free Software Foundation.

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

#include "mpfr-test.h"

static const char * const val[] = {
  "1.0001@100","4.0004000000000@102", "4.0004000000000@97",
  "1.ABF012345@-100","6.afc048d140000@-98","6.afc048d140000@-103",
  "F.FFFFFFFFF@10000","3.fffffffffc000@10003","3.fffffffffc000@9998",
  "1.23456789ABCDEF@42","4.8d159e26af37c@44","4.8d159e26af37c@39",
  "17@42","5.c000000000000@45","5.c000000000000@40",
  "42@-17","1.0800000000000@-13","1.0800000000000@-18"
};

int
main (int argc, char *argv[])
{
  mpfr_t w,z;
  unsigned long k;

  tests_start_mpfr ();

  mpfr_inits2(53, w, z, NULL);

  mpfr_set_inf (w, 1);
  mpfr_mul_2exp (w, w, 10, GMP_RNDZ);
  if (!MPFR_IS_INF(w))
    {
      printf ("Inf != Inf");
      exit (1);
    }

  mpfr_set_nan (w);
  mpfr_mul_2exp (w, w, 10, GMP_RNDZ);
  if (!MPFR_IS_NAN(w))
    {
      printf ("NaN != NaN");
      exit (1);
    }

  for( k = 0 ; k < numberof(val) ; k+=3 )
    {
      mpfr_set_str (w, val[k], 16, GMP_RNDN);
      mpfr_mul_2exp (z, w, 10, GMP_RNDZ);
      if (mpfr_cmp_str(z, val[k+1], 16, GMP_RNDN))
	{
	  printf("ERROR for mpfr_mul_2ui for %s\n", val[k]);
	  printf("Expected: %s\n"
		 "Got     : ", val[k+1]);
	  mpfr_out_str(stdout, 16, 0, z, GMP_RNDN);
	  putchar('\n');
	  exit(-1);
	}
      mpfr_div_2exp (z, w, 10, GMP_RNDZ);
      if (mpfr_cmp_str(z, val[k+2], 16, GMP_RNDN))
        {
          printf("ERROR for mpfr_div_2ui for %s\n"
		 "Expected: %s\n"
                 "Got     : ", val[k], val[k+2]);
          mpfr_out_str(stdout, 16, 0, z, GMP_RNDN);
          putchar('\n');
          exit(-1);
        }
    }
  mpfr_set_ui(w, 1, GMP_RNDN);
  mpfr_mul_2ui(w, w, (unsigned long) MPFR_EXP_MAX+12, GMP_RNDN);
  if (!mpfr_inf_p(w))
    {
      printf("Overflow error!\n");
      exit(1);
    }
  mpfr_set_ui(w, 1, GMP_RNDN);
  mpfr_div_2ui(w, w, (unsigned long) MPFR_EXP_MAX+12, GMP_RNDN);
  if (mpfr_cmp_ui(w, 0))
    {
      printf("Underflow error!\n");
      exit(1);
    }

  mpfr_clears(w,z,NULL);

  tests_end_mpfr ();
  return 0;
}
