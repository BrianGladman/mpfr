/* Test file for mpfr_set_q.

Copyright (C) 2000 PolKA project, Inria Lorraine and Loria

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

#include <stdlib.h>
#include <unistd.h>
#include "gmp.h"
#include "mpfr.h"

main()
{
  long int i, n;
  unsigned long int d;
  mpfr_t x;
  mpq_t q;
  double y, z;
  unsigned char rnd;

  mpfr_init2(x, 53); mpq_init(q);
  for (i=0;i<1000000;i++) {
    n = lrand48();
    d = lrand48();
    if (lrand48()%2); n = -n;
    mpq_set_si(q, n, d);
    rnd = lrand48() % 4;
    mpfr_set_machine_rnd_mode(rnd);
    y = (double) n / d;
    mpfr_set_q(x, q, rnd);
    z = mpfr_get_d(x);
    if (y != z) {
      fprintf(stderr, "Error for q=%ld/%lu and rnd=%s\n", n, d, 
	      mpfr_print_rnd_mode(rnd));
      fprintf(stderr, "libm.a gives %1.20e, mpfr_set_q gives %1.20e\n",
	      y, z);
      exit(1);
    }
  }
  mpfr_clear(x); mpq_clear(q);
  return 0;
}
