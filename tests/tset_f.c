/* Test file for mpfr_set_f.

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
#include "gmp.h"
#include "mpfr.h"
#include "time.h"

#if defined(hpux)
#define srandom srand48
#endif

int
main()
{
  mpfr_t x; mpf_t y; mpf_t z; unsigned long k, pr; 
  
  mpfr_init2(x, 100);
  mpf_init(y); 

  srandom(time(NULL)); 
  mpf_random2(y, 10, 0); 
  mpfr_set_f(x, y, 53, rand() & 3); 

  mpf_clear(y); mpfr_clear(x); 

  for (k = 1; k <= 100000; k++)
    {
      pr = 1 + (rand()&255); 
      mpf_init2(z, pr);
      mpf_random2(z, z->_mp_prec, 0);
      mpfr_init2(x, pr);
      mpfr_set_f(x, z, pr, 0);
    }
  return(0);
}
