#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"

double drand()
{
  double d; long int *i;

  i = (long int*) &d;
  do {
    i[0] = lrand48();
    i[1] = lrand48();
  /*if (lrand48()%2) d=-d; */ /* generates negative numbers */
                              /* useless here */
  } while ((d<1e-153)||(d>1e153));    /* to avoid underflow or overflow
					 in double calculus in sqrt(u*v) */

  return d;
}

/* returns the number of ulp's between a and b */
int ulp(a,b) double a,b;
{
  double eps=1.1102230246251565404e-16; /* 2^(-53) */
  b = (a-b)/a; if (b<0) b = -b;
  return (int) floor(b/eps);
}


void check(double a, unsigned char rnd_mode)
{
  mpfr_t ta, tres;
  double res1,res2 ;

  mpfr_init2(ta, 53);
  mpfr_init2(tres, 53);
  
  mpfr_set_d(ta, a, rnd_mode);

  mpfr_log(tres, ta, rnd_mode);
  mpfr_set_machine_rnd_mode(rnd_mode);
  
  res1=log(a);
  res2=mpfr_get_d(tres);

  if (res1!=res2 && (!isnan(res1) || !isnan(res2))) {
    printf("mpfr_log failed for a=%1.20e, rnd_mode=%d\n",a,rnd_mode);
    printf("expected result is %1.20e, got %1.20e (%d ulp)\n",res1,res2,
	   ulp(res2,res1));
  }
  else
    printf("GOAL !!!\t \t pour le log de %1.20e \n",a);
  mpfr_clear(ta); mpfr_clear(tres); 
  
}

void main() {
  int i;
  double d;
  check(10,GMP_RNDU);
  check(6,GMP_RNDU); 
  check(1,GMP_RNDZ); 
  check(62,GMP_RNDU);
  check(0.5,GMP_RNDZ);  
  check(3,GMP_RNDZ);
  check(234375765,GMP_RNDU);
  check(8,GMP_RNDZ); 
  check(44,GMP_RNDU); 
  check(exp(2),GMP_RNDU);

  for(i=0;i<10;i++) {
    d=drand();
    check(d,GMP_RNDU);
  }
}
