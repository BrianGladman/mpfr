#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"

/* tpi [prec] [rnd] [0 = no print] */

int main(argc, argv) int argc; char *argv[];
{
  mpfr_t x; int p; unsigned char rnd;

  p = (argc>1) ? atoi(argv[1]) : 53;
  rnd = (argc>2) ? atoi(argv[2]) : GMP_RNDZ;
  mpfr_init2(x, p);
  mpfr_pi(x, rnd);
  if (argc>=2) {
    printf("Pi="); mpfr_out_str(stdout, 10, 0, x, rnd); putchar('\n');
  }
  else if (mpfr_get_d(x) != 3.141592653589793116) {
    fprintf(stderr, "mpfr_pi failed for prec=53\n"); exit(1);
  }
  mpfr_clear(x);
  return 0;
}
