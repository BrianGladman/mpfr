#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"
#include "time.h"

extern int random();
extern
#ifdef alpha
int
#else
void
#endif
srandom();

int
main(int argc, char **argv)
{
  mpfr_t x; long k, z, d; unsigned long zl, dl, N; 
  
  mpfr_init2(x, 100);

  srandom(time(NULL)); 

  N = (argc==1) ? 1000000 : atoi(argv[1]);

  for (k = 1; k <= N; k++)
    {
      z = random() - (1 << 30);      
      mpfr_set_si(x, z, GMP_RNDZ); 
      d = (int)mpfr_get_d(x);
      if (d != z)
	fprintf(stderr, "Expected %ld got %ld\n", z, d); exit(1);
      
    }

  for (k = 1; k <= N; k++)
    {
      zl = random();
      mpfr_set_ui(x, zl, GMP_RNDZ); 
      dl = (unsigned int) mpfr_get_d(x);
      if (dl != zl)
	fprintf(stderr, "Expected %lu got %lu\n", zl, dl); exit(1);
    }

  mpfr_clear(x); 
  return(0);
}
