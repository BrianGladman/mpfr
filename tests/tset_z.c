#include "gmp.h"
#include "mpfr.h"

/* tset_z z rnd prec */

check(long i, unsigned char rnd) {
  mpfr_t f; mpz_t z; 

  mpfr_init2(f, 53); mpz_init(z);
  mpz_set_ui(z, i);
  mpfr_set_z(f, z, rnd);
  if ((long)mpfr_get_d(f) != i) {
    printf("Error in mpfr_set_z for i=%ld rnd_mode=%d\n",i,rnd);
    exit(1);
  }
  mpfr_clear(f); mpz_clear(z);
}

main(argc,argv) int argc; char *argv[];
{
  long i, j; unsigned char rnd;

  srand(getpid());
  for (j=0; j<1000000; j++)
    check(lrand48(), rand()%4);
}
