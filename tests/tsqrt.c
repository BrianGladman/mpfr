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

void check(a, rnd_mode) double a; unsigned char rnd_mode;
{
  mpfr_t q, n; double Q,Q2;

#ifdef DEBUG
  printf("a=%1.20e rnd_mode=%d\n",a,rnd_mode);
#endif
  mpfr_init2(q, 53); mpfr_init2(n, 53);
  mpfr_set_d(n, a, rnd_mode);
  mpfr_set_machine_rnd_mode(rnd_mode);
  mpfr_sqrt(q, n, rnd_mode);
  Q = sqrt(a);
  Q2 = mpfr_get_d(q);
#ifdef DEBUG
    printf("expected sqrt is %1.20e, got %1.20e (%d ulp)\n",Q,Q2,
	   ulp(Q2,Q));
    mpfr_print_raw(q); putchar('\n');
#endif
  if (Q!=Q2 && (!isnan(Q) || !isnan(Q2))) {
    printf("mpfr_sqrt failed for a=%1.20e, rnd_mode=%d\n",a,rnd_mode);
    printf("expected sqrt is %1.20e, got %1.20e (%d ulp)\n",Q,Q2,
	   ulp(Q2,Q));
    exit(1);
  }
  mpfr_clear(q); mpfr_clear(n);
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
  check(1.0, 0);
  check(1.0, 1);
  check(3.725290298461914062500000e-9, 0);
  check(3.725290298461914062500000e-9, 1);
  check(1.21902794387441766400e+18, 1);
  check(9.89438396044940256501e-134, 2);
  for (i=0;i<100000;i++) {
    a = drand();
    check(a, rand() % 4);
  }
  return 0;
}
