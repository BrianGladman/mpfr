/* Test file for mpfr_sqrt_ui.

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
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-test.h"

extern int isnan(), getpid();
void check _PROTO((unsigned long, unsigned char, double)); 

int maxulp=0;

void check (unsigned long a, unsigned char rnd_mode, double Q)
{
  mpfr_t q; double Q2; int u, ck;

  mpfr_init2(q, 53);
#ifdef TEST
  mpfr_set_machine_rnd_mode(rnd_mode);
#endif
  mpfr_sqrt_ui(q, a, rnd_mode);
  ck = (Q >= 0.0);
  if (!ck) Q = sqrt(1.0 * a);
  Q2 = mpfr_get_d(q);
  if (Q!=Q2 && (!isnan(Q) || !isnan(Q2))) {
    u = ulp(Q2,Q);
    if (ck) printf("mpfr_sqrt_ui failed");
    else printf("mpfr_sqrt_ui differs from sqrt");
    printf(" for a=%lu, rnd_mode=%s\n",
	   a, mpfr_print_rnd_mode(rnd_mode));
    printf("sqrt gives %1.20e, mpfr_sqrt_ui gives %1.20e (%d ulp)\n",Q,Q2,u);
    exit(1);
  }
  mpfr_clear(q);
}

int main()
{
#ifdef TEST
  int i; unsigned long a;
#ifdef __mips
    /* to get denormalized numbers on IRIX64 */
    union fpc_csr exp;
    exp.fc_word = get_fpc_csr();
    exp.fc_struct.flush = 0;
    set_fpc_csr(exp.fc_word);
#endif

  srand(getpid());
  for (i=0;i<1000000;i++) {
    a = lrand48();
    /* machine arithmetic must agree if a <= 2.0^53 */
    if (1.0*a < 9007199254872064.0) check(a, rand() % 4, -1.0);
  }
#endif
  check(0, GMP_RNDN, 0.0);
  check(2116118, GMP_RNDU, 1.45468828276026215e3);
  return 0;
}
