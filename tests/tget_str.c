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

check(d, rnd) double d; unsigned char rnd;
{
  mpfr_t x; char *str, str2[30]; int l, l2; mp_exp_t e;

  mpfr_init2(x, 53);
  mpfr_set_d(x, d, rnd);
  str = mpfr_get_str(NULL, &e, 10, 5, x, rnd);
  mpfr_set_machine_rnd_mode(rnd);
  sprintf(str2, "%1.4e", d);
  l2 = strlen(str2);
  l = strlen(str); 
  mpfr_clear(x);
  free(str);
}

int
main(int argc, char **argv)
{
  int i; double d;

  srand(getpid());
  /* printf seems to round towards nearest in all cases, at least with gcc */
  check(4.059650008e-83, 0);
  check(-6.606499965302424244461355e233, 0);
  check(-7.4, 0);
  check(0.997, 0);
  check(-4.53063926135729747564e-308, 0);
  check(2.14478198760196000000e+16, 0);
  check(7.02293374921793516813e-84, 0);
  for (i=0;i<100000;i++) {
    do { d = drand(); } while (isnan(d));
    check(d, 0);
  }
}



