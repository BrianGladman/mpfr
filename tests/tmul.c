#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"
#include "mpfr-impl.h"

#define MINNORM 2.2250738585072013831e-308 /* 2^(-1022), smallest normalized */

/* checks that x*y gives the same results in double
   and with mpfr with 53 bits of precision */
void check(double x, double y, unsigned int rnd_mode, unsigned int px, 
unsigned int py, unsigned int pz, double res)
{
  double z1,z2,z3; mpfr_t xx,yy,zz; int i;
  mpf_t xxx,yyy,zzz;

  /* printf("x=%1.20e, y=%1.20e, rnd_mode=%u px=%u py=%u pz=%u\n",x,y,rnd_mode,
     px, py, pz); */
  mpfr_init2(xx, px);
  mpfr_init2(yy, py);
  mpfr_init2(zz, pz);
  mpf_init2(xxx,px); mpf_init2(yyy,py); mpf_init2(zzz,pz);
  mpf_set_d(xxx, x); mpf_set_d(yyy, y);
  mpfr_set_d(xx, x, rnd_mode);
  mpfr_set_d(yy, y, rnd_mode);
for (i=0;i<1;i++)  mpfr_mul(zz, xx, yy, rnd_mode);
  mpf_mul(zzz, xxx, yyy);
  mpfr_set_machine_rnd_mode(rnd_mode);
  z1 = (res==0.0) ? x*y : res;
  z2 = mpfr_get_d(zz);
  z3 = mpf_get_d(zzz);
  if (px==53 && py==53 && pz==53) res=1.0;
  if (res!=0.0 && z1!=z2 && (z1>=MINNORM || z1<=-MINNORM)) {
    printf("expected product is %1.20e, got %1.20e\n",z1,z2);
    printf("mpfr_mul failed for x=%1.20e y=%1.20e with rnd_mode=%u\n",x,y,rnd_mode);
mpfr_print_raw(zz); putchar('\n');
    exit(1);
  }
  mpfr_clear(xx); mpfr_clear(yy); mpfr_clear(zz);
  mpf_clear(xxx); mpf_clear(yyy); mpf_clear(zzz);
}

/* check sign of result */
void check_sign()
{
  mpfr_t a, b;

  mpfr_init2(a, 53); mpfr_init2(b, 53);
  mpfr_set_d(a, -1.0, GMP_RNDN);
  mpfr_set_d(b, 2.0, GMP_RNDN);
  mpfr_mul(a, b, b, GMP_RNDN);
  if (mpfr_get_d(a) != 4.0) {
    fprintf(stderr,"2.0*2.0 gives %1.20e\n", mpfr_get_d(a)); exit(1);
  }
  mpfr_clear(a); mpfr_clear(b);
}

int main(argc,argv) int argc; char *argv[];
{
  double x,y,z; int i,prec,rnd_mode;

  check(6.9314718055994530941514e-1, 0.0, GMP_RNDZ, 53, 53, 53, 0.0);
  check(0.0, 6.9314718055994530941514e-1, GMP_RNDZ, 53, 53, 53, 0.0);
  prec = (argc<2) ? 53 : atoi(argv[1]);
  rnd_mode = (argc<3) ? -1 : atoi(argv[2]);
  check_sign();
  check(-4.165000000e4, -0.00004801920768307322868063274915, GMP_RNDN, 
	53, 53, 53, 2.0); 
  check(2.71331408349172961467e-08, -6.72658901114033715233e-165,
	GMP_RNDZ, 53, 53, 53, 0.0);
  x=0.31869277231188065; y=0.88642843322303122;
  check(x, y, GMP_RNDZ, 53, 53, 53, 0.0);
  x=8.47622108205396074254e-01; y=3.24039313247872939883e-01;
  check(x, y, GMP_RNDU, 28, 45, 1, 0.5);
  x=2.63978122803639081440e-01; 
  y=5736014.0/8388608.0; /* 6.83786096444222835089e-01; */
  check(x, y, GMP_RNDN, 34, 23, 31, 0.180504585267044603);
  x=9.84891017624509146344e-01; /* rounded to 1.0 with prec=6 */
  x=1.0;
  y=1.18351709358762491320e-01;
  check(x, y, GMP_RNDU, 6, 41, 36, 0.1183517093595583);
  /* the following checks that rounding to nearest sets the last
     bit to zero in case of equal distance */
  check(67108865.0, 134217729.0, GMP_RNDN, 53, 53, 53, 0.0);
  x=1.37399642157394197284e-01; y=2.28877275604219221350e-01;
  check(x, y, GMP_RNDN, 49, 15, 32, 0.0314472340833162888);
  x=4.03160720978664954828e-01; y=5.85483042917246621073e-01;
  check(x, y, GMP_RNDZ, 51, 22, 32, 0.2360436821472831);
  x=3.90798504668055102229e-14; y=9.85394674650308388664e-04;
  check(x, y, GMP_RNDN, 46, 22, 12, 0.385027296503914762e-16);
  x=4.58687081072827851358e-01; y=2.20543551472118792844e-01;
  check(x, y, GMP_RNDN, 49, 3, 1, 0.125);
  for (i=0;i<1000000;) {
    x = drand();
    y = drand();
    z = x*y; if (z<0) z=-z;
    if (z<1e+308 && z>1e-308) /* don't test overflow/underflow for now */
      { i++;
      check(x, y, (rnd_mode==-1) ? lrand48()%4 : rnd_mode, 
	    prec, prec, prec, 0.0);
      }
  } 
  return 0;
}

