#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "cputime.h"

extern int mpfr_exp2 (mpfr_ptr, mpfr_srcptr, mp_rnd_t);

int main(int argc,char *argv[])
{
  mpfr_t x, y, z; int prec, rnd, N, i, st;

  srand(getpid());
  if (argc!=4) {
    fprintf(stderr,"Usage: texp2 prec rnd N\n"); exit(1);
  }
  prec = atoi(argv[1]);
  rnd = atoi(argv[2]);
  N = atoi(argv[3]);
  mpfr_init2(x, prec); mpfr_random(x);
  /*  printf("x="); mpfr_print_raw(x); putchar('\n'); */
  mpfr_init2(y, prec);
  mpfr_init2(z, prec);

  mpfr_exp2(z, x, rnd); /* log(2) initialization */

  st=cputime();
  for (i=0; i<N; i++) mpfr_exp2(z, x, rnd);
  printf("mpfr_exp2 takes %dms\n", cputime()-st);

  mpfr_exp(y, x, rnd); /* log(2) initialization */

  st=cputime();
  for (i=0; i<N; i++) mpfr_exp(y, x, rnd);
  printf("mpfr_exp takes %dms\n", cputime()-st);

  if (mpfr_cmp(y,z)) {
    printf("mpfr_exp and mpfr_exp2 disagree for\nx=");
    mpfr_print_raw(x); putchar('\n');
    printf("mpfr_exp gives  "); mpfr_print_raw(y); putchar('\n');
    printf("mpfr_exp2 gives "); mpfr_print_raw(z); putchar('\n');
  }
  return 0;
}


