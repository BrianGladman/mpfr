#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
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


double max(double a,double b) {
  if (a>=b)
    return a;
  return b;
}

double min(double a,double b) {
   if (b>=a)
    return a;
  return b;
}



double dagm(double a, double b) { 
  double u,v,tmpu,tmpv;
  
  if ((isnan(a))||(isnan(b)))
    return a+b;

  tmpv=max(a,b);
  tmpu=min(a,b);

  do
    {
      u=tmpu;
      v=tmpv;
      tmpu=sqrt(u*v);
      tmpv=(u+v)/2.0;
    }
  while (!(((tmpu==u)&&(tmpv==v))||(ulp(u,v)==0))); 

  /*  printf("difference : %i ulp\n",ulp(u,v)); */
	 return u;
}



void check(double a, double b, unsigned char rnd_mode)
{
  mpfr_t ta, tb, tres;
  double res1, res2;

  mpfr_init2(ta, 53);
  mpfr_init2(tb, 53);
  mpfr_init2(tres, 53);
  
  mpfr_set_d(ta, a, rnd_mode);
  mpfr_set_d(tb, b, rnd_mode);

  mpfr_agm(tres, ta, tb, rnd_mode);
  mpfr_set_machine_rnd_mode(rnd_mode);
  
  res1=dagm(a,b);
  res2 = mpfr_get_d(tres);

  if (res1!=res2 && (!isnan(res1) || !isnan(res2))) {
    printf("mpfr_agm failed for a=%1.20e, b=%1.20e, rnd_mode=%d\n",a,b,rnd_mode);
    printf("expected result is %1.20e, got %1.20e (%d ulp)\n",res1,res2,
	   ulp(res2,res1));
    exit(1);
  }
  mpfr_clear(ta); mpfr_clear(tb); mpfr_clear(tres); 
  
}

void main() {
   int i;
   double a,b,gd,pt;

   check(2,1,GMP_RNDN);
   check(6,4,GMP_RNDN); 
   check(62,61,GMP_RNDN);
   check(0.5,1,GMP_RNDN);
   check(1,2,GMP_RNDN); 
   check(234375765,234375000,GMP_RNDN);
   check(8,1,GMP_RNDU);
   check(1,44,GMP_RNDU);  
   check(1,3.725290298461914062500000e-9,GMP_RNDU); 

   srand48(getpid()); 

   for (i=0;i<40;i++) {
     a = drand(); 
     b = drand();
     check(a, b, rand() % 4);
   } 
   printf("fin\n");
}

