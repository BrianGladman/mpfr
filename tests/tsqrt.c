#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"
#ifdef IRIX64
#include <sys/fpu.h>
#endif

extern int isnan(), getpid();

#define check(a,r) check3(a,r,-1.0)

int maxulp=0;

void check3(a, rnd_mode, Q) double a; unsigned char rnd_mode; double Q;
{
  mpfr_t q; double Q2; int ck,u;

  ck = (Q!=-1.0); /* if ck=1, then Q is certified correct */
  mpfr_init2(q, 53);
  mpfr_set_d(q, a, rnd_mode);
  mpfr_set_machine_rnd_mode(rnd_mode);
  mpfr_sqrt(q, q, rnd_mode);
  if (ck==0) Q = sqrt(a);
  Q2 = mpfr_get_d(q);
  if (Q!=Q2 && (!isnan(Q) || !isnan(Q2))) {
    u = ulp(Q2,Q);
    if (ck) {
      printf("mpfr_sqrt failed for a=%1.20e, rnd_mode=%d\n",a,rnd_mode);
      printf("expected sqrt is %1.20e, got %1.20e (%d ulp)\n",Q,Q2,u);
      exit(1);
    }
    else if (u>maxulp || u<-maxulp) {
      maxulp = (u>maxulp) ? u : -u;
      printf("libm.a differs from mpfr_sqrt for a=%1.20e, rnd_mode=%d\n",a,rnd_mode);
      printf("libm.a gives %1.20e, mpfr_sqrt gives %1.20e (%d ulp)\n",Q,Q2,u);
    }
  }
  mpfr_clear(q);
}

int main()
{
  int i; double a;
#ifdef IRIX64
    /* to get denormalized numbers on IRIX64 */
    union fpc_csr exp;
    exp.fc_word = get_fpc_csr();
    exp.fc_struct.flush = 0;
    set_fpc_csr(exp.fc_word);
#endif

  srand(getpid());
  check(6.37983013646045901440e+32, GMP_RNDN);
  check(1.0, GMP_RNDN);
  check(1.0, GMP_RNDZ);
  check(3.725290298461914062500000e-9, GMP_RNDN);
  check(3.725290298461914062500000e-9, GMP_RNDZ);
  a=1190456976439861.0;
  check3(a, GMP_RNDZ, dbl(4630914205854029.0,-27));
  check3(1024.0*a, GMP_RNDZ, dbl(4630914205854029.0,-22));
  check(9.89438396044940256501e-134, GMP_RNDU);
  for (i=0;i<100000;i++) {
    a = drand();
    check(a, rand() % 4);
  }
  return 0;
}
