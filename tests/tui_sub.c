/* Test file for mpfr_ui_sub.

Copyright (C) 2000 Free Software Foundation.

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"
#ifdef TEST
#include "mpfr-test.h"
#endif

void check _PROTO((unsigned long, double, unsigned int, double)); 

/* checks that y/x gives the same results in double
   and with mpfr with 53 bits of precision */
void check (unsigned long y, double x, mp_rnd_t rnd_mode, double z1)
{
  double z2; mpfr_t xx, zz;

  mpfr_init2(xx, 53);
  mpfr_init2(zz, 53);
  mpfr_set_d(xx, x, rnd_mode);
  mpfr_ui_sub(zz, y, xx, rnd_mode);
#ifdef TEST
  mpfr_set_machine_rnd_mode(rnd_mode);
#endif
  if (z1==0.0) z1 = y-x;
  z2 = mpfr_get_d(zz);
  if (z1!=z2 && !(isnan(z1) && isnan(z2))) {
    printf("expected difference is %1.20e, got %1.20e\n",z1,z2);
    printf("mpfr_ui_sub failed for y=%lu x=%1.20e with rnd_mode=%s\n",
	   y, x, mpfr_print_rnd_mode(rnd_mode));
    exit(1);
  }
  mpfr_clear(xx); mpfr_clear(zz);
}

int main(argc,argv) int argc; char *argv[];
{
#ifdef TEST
  double x; unsigned long y, N; int i,rnd_mode,rnd;
#ifdef __mips
    /* to get denormalized numbers on IRIX64 */
    union fpc_csr exp;
    exp.fc_word = get_fpc_csr();
    exp.fc_struct.flush = 0;
    set_fpc_csr(exp.fc_word);
#endif

  srand48(getpid());
  N = (argc<2) ? 1000000 : atoi(argv[1]);
  rnd_mode = (argc<3) ? -1 : atoi(argv[2]);
  for (i=0;i<1000000;i++) {
    x = drand(); 
    y = lrand48();
    if (ABS(x)>2.2e-307) {
      /* avoid denormalized numbers and overflows */
      rnd = (rnd_mode==-1) ? lrand48()%4 : rnd_mode;
      check(y, x, rnd, 0.0);
    }
  }
#endif
  check(1, 1.0/0.0, GMP_RNDN, -1.0/0.0); 
  check(1, -1.0/0.0, GMP_RNDN, 1.0/0.0); 
  check(1, 0.0/0.0, GMP_RNDN, 0.0/0.0); 
  check(1196426492, 1.4218093058435347e-3, GMP_RNDN, 1.1964264919985781e9);
  check(1092583421, -1.0880649218158844e9, GMP_RNDN, 2.1806483428158845901e9);
  check(948002822, 1.22191250737771397120e+20, GMP_RNDN,
	-1.2219125073682338611e20);
  check(832100416, 4.68311314939691330000e-215, GMP_RNDD,
	8.3210041599999988079e8);
  check(1976245324, 1.25296395864546893357e+232, GMP_RNDZ,
	-1.2529639586454686577e232);
  check(2128997392, -1.08496826129284207724e+187, GMP_RNDU,
	1.0849682612928422704e187);
  check(293607738, -1.9967571564050541e-5, GMP_RNDU, 2.9360773800002003e8);
  check(354270183, 2.9469161763489528e3, GMP_RNDN, 3.5426723608382362e8);
  return 0;
}

