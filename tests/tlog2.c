#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"

/* tlog2 [prec] [rnd] [0 = no print] */

int main(argc, argv) int argc; char *argv[];
{
  mpfr_t x; int p; unsigned char rnd;

  p = (argc>1) ? atoi(argv[1]) : 53;
  rnd = (argc>2) ? atoi(argv[2]) : GMP_RNDZ;
  mpfr_init2(x, p);
  mpfr_log2(x, rnd);
  if (argc>=2) {
    printf("log(2)="); mpfr_out_str(stdout, 10, 0, x, rnd); putchar('\n');
  }
  else if (mpfr_get_d(x) != 6.9314718055994530941e-1)  {
    fprintf(stderr, "mpfr_log2 failed for prec=53\n"); exit(1);
  }
  mpfr_clear(x);
  return 0;
}
