#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"
#include <time.h>
#include <math.h>
#ifdef IRIX64
#include <sys/fpu.h>
#endif

extern int isnan();

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
  mpfr_t x,y,z; unsigned long k,n; double d, dd;
#ifdef IRIX64
    /* to get denormalized numbers on IRIX64 */
    union fpc_csr exp;
    exp.fc_word = get_fpc_csr();
    exp.fc_struct.flush = 0;
    set_fpc_csr(exp.fc_word);
#endif

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
      d = drand();
      mpfr_set_d(x, d, 0); 
      dd = mpfr_get_d(x);
      if (d != dd && (!isnan(d) || !isnan(dd)))
	{ 
	  fprintf(stderr, 
		  "Mismatch on : %1.18g != %1.18g\n", d, mpfr_get_d(x)); 
	  mpfr_print_raw(x); putchar('\n');
	} 
    }

  mpfr_clear(x); 
  return 0; 
}
