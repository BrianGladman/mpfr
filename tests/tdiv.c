/* Test file for mpfr_div.

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"
#ifdef __mips
#include <sys/fpu.h>
#endif

extern int isnan();

#define check53(n, d, rnd, res) check4(n, d, rnd, 53, res)

void check4(N, D, rnd_mode, p, Q)
double N, D, Q; unsigned char rnd_mode; int p;
{
  mpfr_t q, n, d; double Q2;

  mpfr_init2(q, p); mpfr_init2(n, p); mpfr_init2(d, p);
  mpfr_set_d(n, N, rnd_mode);
  mpfr_set_d(d, D, rnd_mode);
  mpfr_div(q, n, d, rnd_mode);
#ifdef TEST
  mpfr_set_machine_rnd_mode(rnd_mode);
#endif
  if (Q==0.0) Q = N/D;
  Q2 = mpfr_get_d(q);
  if (p==53 && Q!=Q2 && (!isnan(Q) || !isnan(Q2))) {
    printf("mpfr_div failed for n=%1.20e, d=%1.20e, rnd_mode=%s\n",
	   N, D, mpfr_print_rnd_mode(rnd_mode));
    printf("expected quotient is %1.20e, got %1.20e (%d ulp)\n", Q, Q2,
	   ulp(Q2, Q));
    exit(1);
  }
  mpfr_clear(q); mpfr_clear(n); mpfr_clear(d);  
}

void check24(float N, float D, unsigned char rnd_mode, float Q)
{
  mpfr_t q, n, d; float Q2;

  mpfr_init2(q, 24); mpfr_init2(n, 24); mpfr_init2(d, 24);
  mpfr_set_d(n, N, rnd_mode);
  mpfr_set_d(d, D, rnd_mode);
  mpfr_div(q, n, d, rnd_mode);
  Q2 = mpfr_get_d(q);
  if (Q!=Q2) {
    printf("mpfr_div failed for n=%1.10e, d=%1.10e, prec=24, rnd_mode=%s\n",
	   N, D, mpfr_print_rnd_mode(rnd_mode));
    printf("expected quotient is %1.10e, got %1.10e\n", Q, Q2);
    exit(1);
  }
  mpfr_clear(q); mpfr_clear(n); mpfr_clear(d);  
}

/* the following examples come from the paper "Number-theoretic Test 
   Generation for Directed Rounding" from Michael Parks, Table 2 */
void check_float()
{
  float b=8388608.0; /* 2^23 */

  check24(b*8388610.0, 8388609.0, GMP_RNDN, 8.388609e6);
  check24(b*16777215.0, 16777213.0, GMP_RNDN, 8.388609e6);
  check24(b*8388612.0, 8388611.0, GMP_RNDN, 8.388609e6);
  check24(b*12582914.0, 12582911.0, GMP_RNDN, 8.38861e6);
  check24(12582913.0, 12582910.0, GMP_RNDN, 1.000000238);
  check24(b*16777215.0, 8388609.0, GMP_RNDN, 1.6777213e7);
  check24(b*8388612.0, 8388609.0, GMP_RNDN, 8.388611e6);
  check24(b*12582914.0, 8388610.0, GMP_RNDN, 1.2582911e7);
  check24(b*12582913.0, 8388610.0, GMP_RNDN, 1.258291e7);

  check24(b*8388610.0, 8388609.0, GMP_RNDZ, 8.388608e6);
  check24(b*16777215.0, 16777213.0, GMP_RNDZ, 8.388609e6);
  check24(b*8388612.0, 8388611.0, GMP_RNDZ, 8.388608e6);
  check24(b*12582914.0, 12582911.0, GMP_RNDZ, 8.38861e6);
  check24(12582913.0, 12582910.0, GMP_RNDZ, 1.000000238);
  check24(b*16777215.0, 8388609.0, GMP_RNDZ, 1.6777213e7);
  check24(b*8388612.0, 8388609.0, GMP_RNDZ, 8.38861e6);
  check24(b*12582914.0, 8388610.0, GMP_RNDZ, 1.2582911e7);
  check24(b*12582913.0, 8388610.0, GMP_RNDZ, 1.258291e7);

  check24(b*8388610.0, 8388609.0, GMP_RNDU, 8.388609e6);
  check24(b*16777215.0, 16777213.0, GMP_RNDU, 8.38861e6);
  check24(b*8388612.0, 8388611.0, GMP_RNDU, 8.388609e6);
  check24(b*12582914.0, 12582911.0, GMP_RNDU, 8.388611e6);
  check24(12582913.0, 12582910.0, GMP_RNDU, 1.000000357);
  check24(b*16777215.0, 8388609.0, GMP_RNDU, 1.6777214e7);
  check24(b*8388612.0, 8388609.0, GMP_RNDU, 8.388611e6);
  check24(b*12582914.0, 8388610.0, GMP_RNDU, 1.2582912e7);
  check24(b*12582913.0, 8388610.0, GMP_RNDU, 1.2582911e7);

  check24(b*8388610.0, 8388609.0, GMP_RNDD, 8.388608e6);
  check24(b*16777215.0, 16777213.0, GMP_RNDD, 8.388609e6);
  check24(b*8388612.0, 8388611.0, GMP_RNDD, 8.388608e6);
  check24(b*12582914.0, 12582911.0, GMP_RNDD, 8.38861e6);
  check24(12582913.0, 12582910.0, GMP_RNDD, 1.000000238);
  check24(b*16777215.0, 8388609.0, GMP_RNDD, 1.6777213e7);
  check24(b*8388612.0, 8388609.0, GMP_RNDD, 8.38861e6);
  check24(b*12582914.0, 8388610.0, GMP_RNDD, 1.2582911e7);
  check24(b*12582913.0, 8388610.0, GMP_RNDD, 1.258291e7);
}

void check_convergence()
{
  mpfr_t x, y; int i, j;
  
  mpfr_init2(x, 130); 
  mpfr_set_str_raw(x, "0.1011111101011010101000001010011111101000011100011101010011111011000011001010000000111100100111110011001010110100100001001000111001E6944");
  mpfr_init2(y, 130);
  mpfr_set_ui(y, 5, GMP_RNDN);
  mpfr_div(x, x, y, GMP_RNDD); /* exact division */

  mpfr_set_prec(x, 64);
  mpfr_set_prec(y, 64);
  mpfr_set_str_raw(x, "0.10010010011011010100101001010111100000101110010010101E55");
  mpfr_set_str_raw(y, "0.1E585");
  mpfr_div(x, x, y, GMP_RNDN);
  mpfr_set_str_raw(y, "0.10010010011011010100101001010111100000101110010010101E-529");
  if (mpfr_cmp(x, y)) {
    fprintf(stderr, "Error in mpfr_div for prec=64, rnd=GMP_RNDN\n");
    printf("got        "); mpfr_print_raw(x); putchar('\n');
    printf("instead of "); mpfr_print_raw(y); putchar('\n');
    exit(1);
  }

  for (i=32; i<=64; i+=32) {
      mpfr_set_prec(x, i);
      mpfr_set_prec(y, i);
      mpfr_set_d(x, 1.0, GMP_RNDN);
      for (j=0;j<4;j++) {
	mpfr_set_d(y, 1.0, GMP_RNDN);
	mpfr_div(y, x, y, j);
	if (mpfr_cmp_ui(y, 1)) {
	  fprintf(stderr, "mpfr_div failed for x=1.0, y=1.0, prec=%u rnd=%s\n",
		  i, mpfr_print_rnd_mode(j));
	  printf("got "); mpfr_print_raw(y); putchar('\n');
	  exit(1);
	}
      }
  }

  mpfr_clear(x); mpfr_clear(y);
}

int main(int argc, char *argv[])
{
  int N;
#ifdef TEST
  int i; double n, d, e;
#ifdef __mips
    /* to get denormalized numbers on IRIX64 */
    union fpc_csr exp;
    exp.fc_word = get_fpc_csr();
    exp.fc_struct.flush = 0;
    set_fpc_csr(exp.fc_word);
#endif
#endif

  N = (argc>1) ? atoi(argv[1]) : 100000;
  check_float(); /* checks single precision */
  check_convergence();
  check53(0.0, 1.0, GMP_RNDZ, 0.0);
  check53(-7.4988969224688591e63, 4.8816866450288732e306, GMP_RNDD,
	  -1.5361282826510687291e-243);
  check53(-1.33225773037748601769e+199, 3.63449540676937123913e+79, GMP_RNDZ,
	  -3.6655920045905428978e119);
  check4(4.0, 4.503599627370496e15, GMP_RNDZ, 62, 0.0);
  check4(1.0, 2.10263340267725788209e+187, GMP_RNDU, 65, 0.0);
  check4(2.44394909079968374564e-150, 2.10263340267725788209e+187, GMP_RNDU,
	 65, 0.0);
  check53(9.89438396044940256501e-134, 5.93472984109987421717e-67, GMP_RNDU,
	  1.6672003992376663654e-67);
  check53(9.89438396044940256501e-134, -5.93472984109987421717e-67, GMP_RNDU,
	  -1.6672003992376663654e-67);
  check53(-4.53063926135729747564e-308, 7.02293374921793516813e-84, GMP_RNDD,
	  -6.4512060388748850857e-225);
  check53(6.25089225176473806123e-01, -2.35527154824420243364e-230, GMP_RNDD,
	  -2.6540006635008291192e229);
  check53(6.52308934689126e15, -1.62063546601505417497e273, GMP_RNDN,
	  -4.0250194961676020848e-258);
  check53(1.04636807108079349236e-189, 3.72295730823253012954e-292, GMP_RNDZ,
	  2.810583051186143125e102);
#ifdef TEST
  srand48(getpid());
  for (i=0;i<N;i++) {
    do { n = drand(); d = drand(); e = fabs(n)/fabs(d); }
    /* smallest normalized is 2^(-1022), largest is 2^(1023)*(2-2^(-52)) */
    while (e>=1.7976931348623157081e308 || e<2.225073858507201383e-308);
    check4(n, d, rand() % 4, 53, 0.0);
  }
#endif
  return 0;
}
