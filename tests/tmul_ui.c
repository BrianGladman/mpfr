#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "time.h"

int
main(int argc, char **argv)
{
  mpfr_t x, y;
  
  mpfr_init2(x, 53); mpfr_init2(y, 53);

  /* checks that result is normalized */
  mpfr_set_d(y, 6.93147180559945286227e-01, GMP_RNDZ);
  mpfr_mul_ui(x, y, 1, GMP_RNDZ);
  if (MANT(x)[PREC(x)/BITS_PER_MP_LIMB] >> (BITS_PER_MP_LIMB-1) == 0) {
    fprintf(stderr, "Error in mpfr_mul_ui: result not normalized\n");
    exit(1);
  }
  if (mpfr_cmp(x,y)) {
    fprintf(stderr, "Error in mpfr_mul_ui: 1*y != y\n");
    printf("y=  "); mpfr_print_raw(y); putchar('\n');
    printf("1*y="); mpfr_print_raw(x); putchar('\n');
    exit(1);
  }

  mpfr_set_d(x, 1.0/3.0, GMP_RNDZ); 
  mpfr_mul_ui(x, x, 3, GMP_RNDU); 
  if (mpfr_get_d(x) != 1.0) {
    fprintf(stderr, "U(Z(1/3)*3) does not give 1\n"); exit(1);
  }

  mpfr_clear(x); mpfr_clear(y);
  return(0);
}
