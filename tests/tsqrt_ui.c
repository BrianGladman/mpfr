/* Test file for mpfr_sqrt_ui.

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"
#ifdef IRIX64
#include <sys/fpu.h>
#endif

extern int isnan(), getpid();

int maxulp=0;

void check(a, rnd_mode) unsigned long a; unsigned char rnd_mode;
{
  mpfr_t q; double Q, Q2; int u;

  mpfr_init2(q, 53);
  mpfr_set_machine_rnd_mode(rnd_mode);
  mpfr_sqrt_ui(q, a, rnd_mode);
  Q = sqrt(1.0 * a);
  Q2 = mpfr_get_d(q);
  if (Q!=Q2 && (!isnan(Q) || !isnan(Q2))) {
    u = ulp(Q2,Q);
    printf("mpfr_sqrt_ui failed for a=%lu, rnd_mode=%d\n",a,rnd_mode);
    printf("expected sqrt is %1.20e, got %1.20e (%d ulp)\n",Q,Q2,u);
    exit(1);
  }
  mpfr_clear(q);
}

int main()
{
  int i; unsigned long a;
#ifdef IRIX64
    /* to get denormalized numbers on IRIX64 */
    union fpc_csr exp;
    exp.fc_word = get_fpc_csr();
    exp.fc_struct.flush = 0;
    set_fpc_csr(exp.fc_word);
#endif

  srand(getpid());
  check(0, GMP_RNDN);
  check(2116118, GMP_RNDU);
  for (i=0;i<1000000;i++) {
    a = lrand48();
    /* machine arithmetic must agree if a <= 2.0^53 */
    if (1.0*a < 9007199254872064.0) check(a, rand() % 4);
  }
  return 0;
}
