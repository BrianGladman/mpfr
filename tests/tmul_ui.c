/* Test file for mpfr_mul_ui.

Copyright (C) 1999 Free Software Foundation.

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

#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"
#include "mpfr-impl.h"

int
main (int argc, char *argv[])
{
  mpfr_t x, y;
  unsigned int xprec, yprec;
  
  mpfr_init2(x, 53); mpfr_init2(y, 53);

  /* checks that result is normalized */
  mpfr_set_d(y, 6.93147180559945286227e-01, GMP_RNDZ);
  mpfr_mul_ui(x, y, 1, GMP_RNDZ);
  if (MPFR_MANT(x)[MPFR_PREC(x)/mp_bits_per_limb] >> (mp_bits_per_limb-1) == 0) {
    fprintf(stderr, "Error in mpfr_mul_ui: result not normalized\n");
    exit(1);
  }
  if (mpfr_cmp(x,y)) {
    fprintf(stderr, "Error in mpfr_mul_ui: 1*y != y\n");
    printf("y=  "); mpfr_print_raw(y); putchar('\n');
    printf("1*y="); mpfr_print_raw(x); putchar('\n');
    exit(1);
  }


  mpfr_set_d(x, 1.0/0.0, GMP_RNDZ); 
  mpfr_mul_ui(x, x, 3, GMP_RNDU); 
  if (mpfr_get_d(x) != 1.0/0.0) { 
        fprintf(stderr, "Error in mpfr_mul_ui: Inf*3 does not give Inf\n"); exit(1);
  }

  mpfr_set_d(x, -1.0/0.0, GMP_RNDZ); 
  mpfr_mul_ui(x, x, 3, GMP_RNDU); 
  if (mpfr_get_d(x) != -1.0/0.0) { 
        fprintf(stderr, "Error in mpfr_mul_ui: -Inf*3 does not give -Inf\n"); exit(1);
  }

  mpfr_set_d(x, 0.0/0.0, GMP_RNDZ); 
  mpfr_mul_ui(x, x, 3, GMP_RNDU); 
  if (!isnan(mpfr_get_d(x))) { 
        fprintf(stderr, "Error in mpfr_mul_ui: NaN*3 does not give NaN\n"); exit(1);
  }

  mpfr_set_d(x, 1.0/3.0, GMP_RNDZ); 
  mpfr_mul_ui(x, x, 3, GMP_RNDU); 
  if (mpfr_get_d(x) != 1.0) {
    fprintf(stderr, "Error in mpfr_mul_ui: U(Z(1/3)*3) does not give 1\n"); exit(1);
  }

  /* checks sign is correct */
  mpfr_set_d(x, -2.0, GMP_RNDZ);
  mpfr_set_d(y, 3.0, GMP_RNDZ);
  mpfr_mul_ui(x, y, 4, GMP_RNDZ);
  if (mpfr_cmp_ui(x, 0) <= 0) {
    fprintf(stderr, "Error in mpfr_mul_ui: 4*3.0 does not give a positive result:\n"); 
    mpfr_print_raw(x); putchar('\n');
    printf("mpfr_cmp_ui(x, 0) = %d\n", mpfr_cmp_ui(x, 0));
    exit(1);
  }

  mpfr_set_prec(x, 9);
  mpfr_set_prec(y, 9);
  mpfr_set_str_raw(y,"0.100001111E9");
  mpfr_mul_ui(x, y, 1335, GMP_RNDN);
  mpfr_set_str_raw(y,"0.100111001E19");
  if (mpfr_cmp(x, y)) {
    fprintf(stderr, "Error in mul_ui for 1335*(0.100001111E9)\n"); exit(1);
  }

  mpfr_set_prec(y, 100);
  mpfr_set_prec(x, 100);
  /* y = 1199781142214086656 */
  mpfr_set_str_raw(y, "0.1000010100110011110101001011110010101111000100001E61");
  mpfr_mul_ui(x, y, 121, GMP_RNDD);
  /* 121*y = 145173518207904485376, representable exactly */
  mpfr_set_str_raw(y, "0.1111101111010101111111100011010010111010111110110011001E67");
  if (mpfr_cmp(x, y)) {
    printf("Error for 121*y: expected result is:\n");
    mpfr_print_raw(y); putchar('\n');
  }

  mpfr_set_prec (x, 32);
  mpfr_set_str_raw (x, "0.10000000000000000000000000000000E1");
  mpfr_set_prec (y, 93);
  mpfr_mul_ui (y, x, 1, GMP_RNDN);

  mpfr_set_prec (x, 287);
  mpfr_set_str_raw (x, "0.1111E7");
  mpfr_set_prec (y, 289);
  mpfr_mul_ui (y, x, 6, GMP_RNDN);
  mpfr_set_str_raw (x, "0.101101E10");
  if (mpfr_cmp (x, y))
    {
      printf ("Error for 6 * 120\n");
      exit (1);
    }

  for (xprec = 53; xprec <= 128; xprec++)
    {
      mpfr_set_prec (x, xprec);
      mpfr_set_str_raw (x, "0.1100100100001111110011111000000011011100001100110111E2");
      for (yprec = 53; yprec <= 128; yprec++)
	{
	  mpfr_set_prec (y, yprec);
	  mpfr_mul_ui (y, x, 1, GMP_RNDN);
	  if (mpfr_get_d (x) != mpfr_get_d (y))
	    {
	      fprintf (stderr, "multiplication by 1.0 fails for xprec=%u, yprec=%u\n", xprec, yprec);
	      printf ("expected "); mpfr_print_raw (x); putchar ('\n');
	      printf ("got      "); mpfr_print_raw (y); putchar ('\n');
	      exit (1);
	    }
	}
    }

  mpfr_clear(x); mpfr_clear(y);

  return 0;
}
