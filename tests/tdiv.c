#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

/* #define DEBUG */

double drand()
{
  double d; long int *i;

  i = (long int*) &d;
  i[0] = lrand48();
  i[1] = lrand48();
  if (lrand48()%2) d=-d; /* generates negative numbers */
  return d;
}

/* returns the number of ulp's between a and b */
int ulp(a,b) double a,b;
{
  double eps=1.1102230246251565404e-16; /* 2^(-53) */
  b = (a-b)/a; if (b<0) b = -b;
  return (int) floor(b/eps);
}

#define check(n,d,r) check4(n,d,r,53)

void check4(N, D, rnd_mode, p) double N, D; unsigned char rnd_mode; int p;
{
  mpfr_t q, n, d; double Q,Q2;

#ifdef DEBUG
  printf("N=%1.20e D=%1.20e rnd_mode=%d\n",N,D,rnd_mode);
#endif
  mpfr_init2(q, p); mpfr_init2(n, p); mpfr_init2(d, p);
  mpfr_set_d(n, N, rnd_mode);
  mpfr_set_d(d, D, rnd_mode);
  mpfr_div(q, n, d, rnd_mode);
  mpfr_set_machine_rnd_mode(rnd_mode);
  Q = N/D;
  Q2 = mpfr_get_d(q);
#ifdef DEBUG
    printf("expected quotient is %1.20e, got %1.20e (%d ulp)\n",Q,Q2,
	   ulp(Q2,Q));
    mpfr_print_raw(q); putchar('\n');
#endif
  if (Q!=Q2 && (!isnan(Q) || !isnan(Q2))) {
    printf("mpfr_div failed for n=%1.20e, d=%1.20e, rnd_mode=%d\n",N,D,rnd_mode);
    printf("expected quotient is %1.20e, got %1.20e (%d ulp)\n",Q,Q2,
	   ulp(Q2,Q));
    exit(1);
  }
  mpfr_clear(q); mpfr_clear(n); mpfr_clear(d);  
}

check_convergence()
{
  mpfr_t x, y;
  
  mpfr_init2(x, 130); 
  mpfr_set_str_raw(x, "0.1011111101011010101000001010011111101000011100011101010011111011000011001010000000111100100111110011001010110100100001001000111001E6944");
  mpfr_init_set_ui(y, 5, 130, GMP_RNDN);
  mpfr_div(x, x, y, GMP_RNDD); /* exact division */
  mpfr_clear(x); mpfr_clear(y);
}

void main()
{
  int i; double n, d, e; 

  check_convergence();
  check(0.0, 1.0, 1);
  check(-1.33225773037748601769e+199, 3.63449540676937123913e+79, 1);
  d = 1.0; for (i=0;i<52;i++) d *= 2.0;
  check4(4.0, d, GMP_RNDZ, 62);
  check4(1.0, 2.10263340267725788209e+187, 2, 65);
  check4(2.44394909079968374564e-150, 2.10263340267725788209e+187, 2, 65);
  /* the following tests when d is an exact power of two */
  check(9.89438396044940256501e-134, 5.93472984109987421717e-67, 2);
  check(9.89438396044940256501e-134, -5.93472984109987421717e-67, 2);
  check(-4.53063926135729747564e-308, 7.02293374921793516813e-84, 3);
  check(6.25089225176473806123e-01, -2.35527154824420243364e-230, 3);
  check(6.52308934689126000000e+15, -1.62063546601505417497e+273, 0);
  check(1.04636807108079349236e-189, 3.72295730823253012954e-292, 1);
  for (i=0;i<100000;i++) {
    do { n = drand(); d = drand(); e = fabs(n)/fabs(d); }
    /* smallest normalized is 2^(-1022), largest is 2^(1023)*(2-2^(-52)) */
    while (e>=1.7976931348623157081e308 || e<2.225073858507201383e-308);
    check(n, d, rand() % 4);
  }
}
