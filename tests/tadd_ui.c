/* Test file for mpfr_add_ui

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

/* #define DEBUG */
/* #define VERBOSE */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"
#include "mpfr-impl.h"
#ifdef IRIX64
#include <sys/fpu.h>
#endif

extern int isnan();
extern int getpid();

#define ABS(x) (((x)>0) ? (x) : (-x))

/* checks that x+y gives the same results in double
   and with mpfr with 53 bits of precision */
void check(double x, unsigned long y, unsigned int rnd_mode)
{
  double z1,z2; mpfr_t xx,zz;

  mpfr_init2(xx, 53);
  mpfr_init2(zz, 53);
  mpfr_set_d(xx, x, rnd_mode);
  mpfr_add_ui(zz, xx, y, rnd_mode);
  mpfr_set_machine_rnd_mode(rnd_mode);
  z1 = x+y;
  z2 = mpfr_get_d(zz);
  if (z1!=z2 && !(isnan(z1) && isnan(z2))) {
    printf("expected sum is %1.20e, got %1.20e\n",z1,z2);
    printf("mpfr_add_ui failed for x=%1.20e y=%lu with rnd_mode=%u\n",x,y,rnd_mode);
    exit(1);
  }
  mpfr_clear(xx); mpfr_clear(zz);
}

int main(argc,argv) int argc; char *argv[];
{
  double x; unsigned long y, N; int i,rnd_mode,rnd;
#ifdef IRIX64
    /* to get denormalized numbers on IRIX64 */
    union fpc_csr exp;
    exp.fc_word = get_fpc_csr();
    exp.fc_struct.flush = 0;
    set_fpc_csr(exp.fc_word);
#endif

  check(1.22191250737771397120e+20, 948002822, GMP_RNDN);
  srand(getpid());
  N = (argc<2) ? 1000000 : atoi(argv[1]);
  rnd_mode = (argc<3) ? -1 : atoi(argv[2]);
  for (i=0;i<1000000;i++) {
    x = drand(); 
    y = lrand48();
    if (ABS(x)>2.2e-307 && x+y<1.7e+308 && x+y>-1.7e308) {
      /* avoid denormalized numbers and overflows */
      rnd = (rnd_mode==-1) ? lrand48()%4 : rnd_mode;
      check(x, y, rnd);
    }
  } 
  return 0;
}

