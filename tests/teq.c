/* Test file for mpfr_cmp2.

Copyright (C) 1999 PolKA project, Inria Lorraine and Loria

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
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"
#include "mpfr-impl.h"

void teq(mpfr_srcptr x)
{
  mpfr_t y; long k, px, mx; 

  mpfr_init2(y, PREC(x)); 

  mx = (PREC(x) - 1)/BITS_PER_MP_LIMB; 
  px = BITS_PER_MP_LIMB - 2; 

  for (k = 2; k < PREC(x); k++)
    {
      mpfr_set(y, x, GMP_RNDN); 

      MANT(y) [mx] ^= (mp_limb_t) 1 << px; 

      if (mpfr_eq(y, x, k) ||  
	  !mpfr_eq(y, x, k - 1))
	{
	  fprintf(stderr, "Error in eq.\n"); 	  
	  printf("x = "); mpfr_print_raw(x); printf("\n"); 
	  printf("y = "); mpfr_print_raw(y); printf("\n"); 	
	  printf("k = %ld\n", k); 
	  printf("mpfr_eq(y, x, k) = %d\nmpfr_eq(y, x, k - 1) = %d\n", mpfr_eq(y, x, k),mpfr_eq(y, x, k - 1)); 
	  mpfr_clear(x); mpfr_clear(y); 
	  exit(-1); 
	}

      if (px) { --px; } else { --mx; px = BITS_PER_MP_LIMB - 1; }
    }
  mpfr_clear(y); 
}

int main()
{
  int j; mpfr_t x; 

  mpfr_init2(x, 1000); 

  for (j=0;j<1000;j++) {
    mpfr_random(x); 
    teq(x); 
  }
  mpfr_clear (x); 
  return 0;
}

