#include "gmp.h"
#include "mpfr.h"

/* tlog2 [prec] [rnd] [0 = no print] */

main(argc, argv) int argc; char *argv[];
{
  mpfr_t x; int p; unsigned char rnd;

  p = (argc>1) ? atoi(argv[1]) : 53;
  rnd = (argc>2) ? atoi(argv[2]) : GMP_RNDZ;
  mpfr_init2(x, p);
  mpfr_log2(x, rnd);
  if (argc<=3 || atoi(argv[3])!=0) {
    printf("log(2)="); mpfr_out_str(stdout, 10, 0, x, rnd); putchar('\n');
  }
  mpfr_clear(x);
}
