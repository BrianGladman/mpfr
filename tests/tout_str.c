#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gmp.h"
#include "mpfr.h"
#include <time.h>

void print_double(d) double d;
{
  int e, i;

  e = (int) ceil(log(fabs(d))/log(2.0));
  /* d <= 2^e */
  e -= 53;
  if (e>0) for (i=0;i<e;i++) d /= 2.0;
  else for (i=0;i<-e;i++) d *= 2.0;
  printf("%1.0f*2^(%d)",d,e);
}

double drand()
{
  double d; long int *i;

  i = (long int*) &d;
  i[0] = lrand48();
  i[1] = lrand48();
  if (lrand48()%2) d=-d; /* generates negative numbers */
  return d;
}

check(d, rnd, base) double d; unsigned char rnd; int base;
{
  mpfr_t x;

  mpfr_init2(x, 53);
  mpfr_set_d(x, d, rnd);
  mpfr_set_machine_rnd_mode(rnd);
  printf("%1.19e base %d:\n ", d, base);
  mpfr_out_str(stdout, base, (base==2) ? 53 : 0, x, rnd);
  putchar('\n');
  if (base==2) { mpfr_print_raw(x); putchar('\n'); }
  mpfr_clear(x);
}

int
main(int argc, char **argv)
{
  int i; double d;

  srand(getpid());
  /* printf seems to round towards nearest in all cases, at least with gcc */
  check(4.059650008e-83, 0, 10);
  check(-6.606499965302424244461355e233, 0, 10);
  check(-7.4, 0, 10);
  check(0.997, 0, 10);
  check(-4.53063926135729747564e-308, 0, 10);
  check(2.14478198760196000000e+16, 0, 10);
  check(7.02293374921793516813e-84, 0, 10);
  for (i=0;i<100;i++) {
    do { d = drand(); } while (isnan(d));
    check(d, 0, 2 + rand()%35);
  }
}



