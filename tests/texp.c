#include <math.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

/* #define DEBUG */

extern void srand48();
extern double drand48();
extern int isnan();
extern int getpid();

int maxu=0;

/* returns the number of ulp's between a and b */
int ulp(a,b) double a,b;
{
  double eps=1.1102230246251565404e-16; /* 2^(-53) */
  b = (a-b)/a; if (b<0) b = -b;
  return (int) floor(b/eps);
}

#define check(d, r) check3(d, r, 0.0)

/* returns the number of ulp of error */
int check3(double d, unsigned char rnd, double e)
{
  mpfr_t x, y; double f; int u=0, ck=0;

#ifdef DEBUG
  printf("d=%1.20e rnd=%d\n",d,rnd);
#endif
  mpfr_init2(x, 53); mpfr_init2(y, 53);
  mpfr_set_machine_rnd_mode(rnd);
  if (e==0.0) e = exp(d); else ck=1; /* really check */
  mpfr_set_d(x, d, rnd); 
  mpfr_exp(y, x, rnd); 
  f = mpfr_get_d(y);
  if (f != e && (!isnan(f) || !isnan(e))) {
    u = ulp(e,f);
    if (u<0) {
      if (u == (mp_limb_t)1<<(BITS_PER_MP_LIMB-1)) u += 1;
      u=-u;
    }
    if (u!=0) {
      if (ck) {
	printf("mpfr_exp failed for x=%1.20e, rnd=%d\n",d,rnd);
	printf("expected result is %1.20e, got %1.20e, dif=%d ulp\n",e,f,u);
	exit(1);
      }
      else if (u>maxu) {
	maxu=u;
	printf("mpfr_exp differs from libm.a for x=%1.20e, rnd=%d\n",d,rnd);
	printf("libm.a gave %1.20e, mpfr_exp got %1.20e, dif=%d ulp\n",e,f,u);
      }
    }
  }
  mpfr_clear(x); mpfr_clear(y);
  return u;
}

/* computes n bits of exp(d) */
int check_large (double d, int n, char rnd)
{
  mpfr_t x; mpfr_t y;
  
  mpfr_init2(x, n); mpfr_init2(y, n);
  if (d==0.0) { /* try exp(Pi*sqrt(163)/3)-640320 */
    mpfr_set_d(x, 163.0, rnd);
    mpfr_sqrt(x, x, rnd);
    mpfr_pi(y, rnd);
    mpfr_mul(x, x, y, rnd);
    mpfr_div_ui(x, x, 3, rnd);
  }
  else mpfr_set_d(x, d, rnd);
  mpfr_exp(x, x, rnd);
  if (d==0.0) {
    mpfr_set_d(y, 640320.0, rnd);
    mpfr_sub(x, x, y, rnd);
    printf("exp(Pi*sqrt(163)/3)-640320=");
  }
  else printf("exp(%1.20e)=",d); 
  mpfr_out_str(stdout, 10, 0, x, rnd);
  putchar('\n');
  if (n==53) printf(" =%1.20e\n", mpfr_get_d(x));
  mpfr_clear(x); mpfr_clear(y);
  return 0;
}

/* expx is the value of exp(X) rounded towards -infinity */
int check_worst_case(double X, double expx)
{
  mpfr_t x, y;

  mpfr_init2(x, 53); mpfr_init2(y, 53);
  mpfr_set_d(x, X, GMP_RNDN);
#ifdef DEBUG
  printf("x="); mpfr_print_raw(x); putchar('\n');
#endif
  mpfr_exp(y, x, GMP_RNDD);
#ifdef DEBUG
  printf("D(exp(x))="); mpfr_print_raw(y); putchar('\n');
#endif
  if (mpfr_get_d(y) != expx) {
    fprintf(stderr, "exp(x) rounded towards -infinity is wrong\n"); exit(1);
  }
  mpfr_exp(x, x, GMP_RNDN);
#ifdef DEBUG
  printf("N(exp(x))="); mpfr_print_raw(x); putchar('\n');
#endif
  mpfr_set_d(x, X, GMP_RNDN);
  mpfr_exp(x, x, GMP_RNDU);
#ifdef DEBUG
  printf("U(exp(x))="); mpfr_print_raw(x); putchar('\n');
#endif
  mpfr_add_one_ulp(y);
  if (mpfr_cmp(x,y)) {
    fprintf(stderr, "exp(x) rounded towards +infinity is wrong\n"); exit(1);
  }
  mpfr_clear(x); mpfr_clear(y);
  return 0;
}

/* worst cases communicated by Jean-Michel Muller and Vincent Lefevre */
int check_worst_cases()
{
  mpfr_t x;

  mpfr_init2(x, 53);
  check_worst_case(4.44089209850062517562e-16, 1.00000000000000022204);
  check_worst_case(6.39488462184069720009e-14, 1.00000000000006372680);
  check_worst_case(1.84741111297455401935e-12, 1.00000000000184718907);
  check_worst_case(1.76177628026265550074e-10, 1.00000000017617751702);
  check3(1.76177628026265550074e-10, GMP_RNDN, 1.00000000017617773906);
  check_worst_case(7.54175277499595900852e-10, 1.00000000075417516676);
  check3(7.54175277499595900852e-10, GMP_RNDN, 1.00000000075417538881);
  mpfr_set_str_raw(x, "1.1001111010011100101110111111110101100000100000001011e-31");
#ifdef DEBUG
  printf("x=%1.20e\n", mpfr_get_d(x));
  printf(" ="); mpfr_print_raw(x); putchar('\n');
#endif
  mpfr_clear(x);
  return 0;
}

int
main(int argc, char **argv)
{
  int i, N, s=0, e, maxe=0; double d, lo, hi;

  if (argc==4) { check_large(atof(argv[1]), atoi(argv[2]), atoi(argv[3])); 
		 exit(1); }
  check_worst_cases();
  check3(-8.88024741073346941839e-17, GMP_RNDU, 1.0);
  check3(8.70772839244701057915e-01, GMP_RNDN, 2.38875626491680437269);
  check3(1.0, GMP_RNDN, 2.71828182845904509080);
  check3(-3.42135637628104173534e-07, GMP_RNDZ, 0.999999657864420798958);
  srand48(getpid());
  N = (argc==1) ? 0 : atoi(argv[1]);
  lo = (argc>=3) ? atof(argv[2]) : -7.083964185e2;
  hi = (argc>=4) ? atof(argv[3]) : 7.097827129e2;
  for (i=0;i<N;i++) {
    /* select d such that exp(d) can be represented as a normalized
       machine double-precision number, 
       i.e. 2^(-1022) <= exp(d) <= 2^(1023)*(2-2^(-52)) */
    d = lo + (hi-lo)*drand48();
    e = check(d, rand() % 4);
    s += e;
    if (e>maxe) maxe=e;
  }
  if (N) printf("mean error=%1.2e max error=%d\n", (double)s/(double)N,maxe);
  check3(2.26523754332090625496e+01, 3, 6.8833785261699581146e9);
  /* errors found in libm.a on PC under Linux */
  check3(1.31478962104089092122e+01, GMP_RNDZ, 5.12930793917860137299e+05);
  check3(4.25637507920002378103e-01, GMP_RNDU, 1.53056585656161181497e+00);
  check3(6.26551618962329307459e-16, GMP_RNDU, 1.00000000000000066613e+00);
  check3(-3.35589513871216568383e-03, GMP_RNDD, 9.96649729583626853291e-01);
  check3(1.95151388850007272424e+01, GMP_RNDU, 2.98756340674767792225e+08);
  check3(2.45045953503350730784e+01, GMP_RNDN, 4.38743344916128387451e+10);
  check3(2.58165606081678085104e+01, GMP_RNDD, 1.62925781879432281494e+11);
  check3(-2.36539020084338638128e+01, GMP_RNDZ, 5.33630792749924762447e-11);
  check3(2.39211946135858077866e+01, GMP_RNDU, 2.44817704330214385986e+10);
  check3(-2.78190533055889162029e+01, 1, 8.2858803483596879512e-13); 
           /* +45 ulp, wrong side */
  check3(2.64028186174889789584e+01, 3, 2.9281844652878973388e11); /* -45 ulp*/
  check3(2.92086338843268329413e+01, 1, 4.8433797301907177734e12); /* -45 ulp*/
  check3(-2.46355324071459982349e+01, 1, 1.9995129297760994791e-11);
           /* +45 ulp, wrong side */
  check3(-2.23509444608605427618e+01, 1, 1.9638492867489702307e-10);
           /* +45 ulp, wrong side */
  check3(-2.41175390197331687148e+01, 3, 3.3564940885530624592e-11);/*-45 ulp*/
  check3(2.46363885231578088053e+01, 2, 5.0055014282693267822e10); /* +45 ulp*/
  d=7819821913254249.0; d /= 70368744177664.0;
  check3(d, GMP_RNDN, 1.8262572323517295459e48); /*-73 ulp*/
  check3(-3.56196340354684821250e+02, 0, 2.0225297096141478156e-155); /*+352 */
  check3(6.59678273772710895173e+02, 2, 3.1234469273830195529e286); /* +459 */
  check3(5.13772529701934331570e+02, 3, 1.3445427121297197752e223); /* -469 */
  check3(3.57430211008718345056e+02, 3, 1.6981197246857298443e155); /* -610 */
  check3(3.82001814471465536371e+02, 2, 7.9667300591087367805e165); /* +705 */
  check3(5.92396038219384422518e+02, 3, 1.880747529554661989e257); /* -707 */
  check3(-5.02678550462488090034e+02, 2, 4.8919201895446217839e-219); /* +708*/
  check3(5.30015757134837031117e+02, 3, 1.5237672861171573939e230); /* -709 */
  check3(5.16239362447650933063e+02, 1, 1.5845518406744492105e224); /* -710 */
  /* between 1/2 and 1 */
  check3(6.00812634798592370977e-01, 0, 1.823600119339019443); /* +1 ulp */
  return 0;
}
