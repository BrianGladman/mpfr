#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

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

void check(a, rnd_mode) double a; unsigned char rnd_mode;
{
  mpfr_t q, n; double Q,Q2;

#ifdef DEBUG
  printf("a=%1.20e rnd_mode=%d\n",a,rnd_mode);
#endif
  mpfr_init2(q, 53); mpfr_init2(n, 53);
  mpfr_set_d(n, a, rnd_mode);
  mpfr_set_machine_rnd_mode(rnd_mode);
  mpfr_sqrt(q, n, rnd_mode);
  Q = sqrt(a);
  Q2 = mpfr_get_d(q);
#ifdef DEBUG
    printf("expected sqrt is %1.20e, got %1.20e (%d ulp)\n",Q,Q2,
	   ulp(Q2,Q));
    mpfr_print_raw(q); putchar('\n');
#endif
  if (Q!=Q2 && (!isnan(Q) || !isnan(Q2))) {
    printf("mpfr_sqrt failed for a=%1.20e, rnd_mode=%d\n",a,rnd_mode);
    printf("expected sqrt is %1.20e, got %1.20e (%d ulp)\n",Q,Q2,
	   ulp(Q2,Q));
    exit(1);
  }
  mpfr_clear(q); mpfr_clear(n);
}

void main()
{
  int i; double a;

  srand(getpid());
  check(1.21902794387441766400e+18, 1);
  check(9.89438396044940256501e-134, 2);
  for (i=0;i<100000;i++) {
    do { a = drand(); } while (isnan(a));
    check(a, rand() % 4);
  }
}
