#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "gmp.h"
#include "mpfr.h"
#include "mpfr-impl.h"

#define check(d,r,b) check4(d,r,b,53)

void check4(d, rnd, base, prec) double d; unsigned char rnd; int base, prec;
{
  mpfr_t x;

  mpfr_init2(x, prec);
  mpfr_set_d(x, d, rnd);
mpfr_print_raw(x); printf("\n");
  mpfr_set_machine_rnd_mode(rnd);
  printf("%1.19e base %d:\n ", d, base);
  mpfr_out_str(stdout, base, (base==2) ? prec : 0, x, rnd);
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
  check4(1.0, GMP_RNDN, 10, 120);
  check(1.0, GMP_RNDU, 10);
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
  return 0;
}



