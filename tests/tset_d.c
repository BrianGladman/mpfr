#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"
#include <time.h>
#include <math.h>

double drand()
{
  double d; long int *i;

  i = (long int*) &d;
  i[0] = lrand48();
  i[1] = lrand48();
  return d;
}

int
main(int argc, char **argv)
{
  mpfr_t x,y,z; unsigned long k,n; double d; 

  mpfr_init2(z, 32);
  mpfr_set_d(z, 1.0, 0);
  if (mpfr_get_d(z) != 1.0) {
    mpfr_print_raw(z); putchar('\n');
    printf("Error: 1.0 != 1.0\n"); exit(1);
  }
  mpfr_init2(x, 53); mpfr_init2(y, 53);
  mpfr_set_d(x, d=-1.08007920352320089721e+150, 0);
  if (mpfr_get_d(x) != d) {
    mpfr_print_raw(x); putchar('\n');
    printf("Error: get_d o set_d <> identity for d = %1.20e %1.20e\n",d,
	   mpfr_get_d(x)); exit(1);
  }
  srand48(time(NULL)); 
  mpfr_set_d(x, 8.06294740693074521573e-310, 0); 
  d = -6.72658901114033715233e-165;
  mpfr_set_d(x, d, 0);
  if (d != mpfr_get_d(x)) {
    mpfr_print_raw(x); putchar('\n');
    printf("Error: get_d o set_d <> identity for d = %1.20e %1.20e\n",d,
	   mpfr_get_d(x)); exit(1);
  }
  n = (argc==1) ? 1000000 : atoi(argv[1]);
  for (k = 1; k <= n; k++)
    {      
      do { d = drand(); } while (isnan(d)); /* does not yet work for NaN */
      mpfr_set_d(x, d, 0); 
      if (d != mpfr_get_d(x)) 
	{ 
	  fprintf(stderr, 
		  "Mismatch on : %1.18g != %1.18g\n", d, mpfr_get_d(x)); 
	  mpfr_print_raw(x); putchar('\n');
	} 
    }

  mpfr_clear(x); 
  return 0; 
}
