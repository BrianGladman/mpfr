#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmp.h"
#include "longlong.h"
#include "mpfr.h"

main()
{
  mpfr_t x; unsigned long i; long s;
  
  mpfr_init(x);

  mpfr_set_ui(x, 3, GMP_RNDZ);
  if (mpfr_cmp_ui(x, i=3)!=0) {
    printf("Error in mpfr_cmp_ui(%1.20f,%d)\n",mpfr_get_d(x), i); exit(1);
  }
  if (mpfr_cmp_ui(x, i=2)<=0) {
    printf("Error in mpfr_cmp_ui(%1.20f,%d)\n",mpfr_get_d(x), i); exit(1);
  }
  if (mpfr_cmp_ui(x, i=4)>=0) {
    printf("Error in mpfr_cmp_ui(%1.20f,%d)\n",mpfr_get_d(x), i); exit(1);
  }

  mpfr_set_si(x, -3, GMP_RNDZ);
  if (mpfr_cmp_si(x, s=-3)!=0) {
    printf("Error in mpfr_cmp_si(%1.20f,%d)\n",mpfr_get_d(x), s); exit(1);
  }
  if (mpfr_cmp_si(x, s=-4)<=0) {
    printf("Error in mpfr_cmp_si(%1.20f,%d)\n",mpfr_get_d(x), s); exit(1);
  }
  if (mpfr_cmp_si(x, s=1)>=0) {
    printf("Error in mpfr_cmp_si(%1.20f,%d)\n",mpfr_get_d(x), s); exit(1);
  }
}
