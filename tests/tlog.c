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
   } while ((d<1e-153)||(d>1e153));    /* to avoid underflow or overflow
					 in double calculus in sqrt(u*v) */
  return d;
}


/* returns the number of ulp's between a and b */
int ulp(a,b) double a,b;
{
  double eps=1.1102230246251565404e-16; /* 2^(-53) */
  b = (a-b)/a;
  return (int) floor(b/eps);
}


void check(double a, unsigned char rnd_mode)
{
  mpfr_t ta, tres;
  double res1,res2 ;

  mpfr_set_machine_rnd_mode(rnd_mode);  
  res1=log(a);

  mpfr_init2(ta, 53);
  mpfr_init2(tres, 53);
  mpfr_set_d(ta, a, GMP_RNDN);

  mpfr_log(tres, ta, rnd_mode);
  res2=mpfr_get_d(tres);


  if (res1!=res2 && (!isnan(res1) || !isnan(res2))) {
    printf("mpfr_log failed for    a=%1.20e, rnd_mode=%d\n",a,rnd_mode);
    printf(" double calculus gives %1.20e\n mpfr_log        gives %1.20e (%d ulp)\n pari            gives \n \n",res1,res2,ulp(res1,res2));
  }
  /*else {
    printf("GOAL !!! for           a=%1.20e, rnd_mode=%d\n",a,rnd_mode);
    printf(" double calculus gives %1.20e\n pari            gives \n \n",res1);
    }*/
  mpfr_clear(ta); mpfr_clear(tres); 
}

check3(double d, unsigned long prec, unsigned char rnd)
{
  mpfr_t x, y;
  
  mpfr_init2(x, prec); mpfr_init2(y, prec);
  mpfr_set_d(x, d, rnd);
  mpfr_log(y, x, rnd);
  mpfr_out_str(stdout, 10, 0, y, rnd); putchar('\n');
  mpfr_clear(x); mpfr_clear(y);
}

void main(int argc, char *argv[]) {
  int i;
  double d;

  if (argc==4) { /* tlog x prec rnd */
    check3(atof(argv[1]), atoi(argv[2]), atoi(argv[3]));
    return;
  }
  printf("SUN Solaris: craffe\n 20000 essais\n");
  printf("GMP_RNDN : %i, GMP_RNDZ : %i,GMP_RNDU : %i,GMP_RNDD : %i\n",GMP_RNDN, GMP_RNDZ,GMP_RNDU, GMP_RNDD); 
   

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
  check(7.53428236571286402512e+34,GMP_RNDZ);
  check(6.18784121531737948160e+19,GMP_RNDZ); 
  check(1.02560267603047283735e+00,GMP_RNDD);
  check(7.53428236571286402512e+34,GMP_RNDZ);
  check(1.42470900831881198052e+49,GMP_RNDZ); 
  
  check(1.08013816255293777466e+11,GMP_RNDN);
  check(6.72783635300509015581e-37,GMP_RNDU);
  check(2.25904918906057891180e-52,GMP_RNDU);
  check(1.48901209246462951085e+00,GMP_RNDD);
  check(1.70322470467612341327e-01,GMP_RNDN);
  check(1.94572026316065240791e+01,GMP_RNDD);
  check(4.01419512207026418764e+04,GMP_RNDD);
  check(9.47077365236487591672e-04,GMP_RNDZ);
  check(3.95906157687589643802e-109,GMP_RNDD);
  check(2.73874914516503004113e-02,GMP_RNDD);
  check(9.18989072589566467669e-17,GMP_RNDZ);
  check(7.70603645360819104595e+54,GMP_RNDZ);
  check(1.74827399630587801934e-23,GMP_RNDZ);
  check(4.35302958401482307665e+22,GMP_RNDD);
  check(9.70791868689332915209e+00,GMP_RNDD);
  check(2.22183639799464011100e-01,GMP_RNDN);
  check(2.27313466156682375540e+00,GMP_RNDD);
  check(6.58057413965851156767e-01,GMP_RNDZ);
  check(7.34302197248998461006e+43,GMP_RNDZ);
  check(6.09969788341579732815e+00,GMP_RNDD);
  check(8.94529798779875738679e+82,GMP_RNDD);
  check(1.68775280934272742250e+00,GMP_RNDZ);





  srand48(getpid());
  for(i=0;i<20000;i++) {
    d=drand();
    check(d,rand() % 4);
  }
}

