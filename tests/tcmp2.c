#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmp.h"
#include "longlong.h"
#include "mpfr.h"
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

void tcmp2(x, y, i) double x, y; int i;
{
  mpfr_t xx,yy; int j;

  if (i==-1) i = (int) floor(log(x)/log(2.0)) - (int) floor(log(x-y)/log(2.0));
  mpfr_init2(xx, 53); mpfr_init2(yy, 53);
  mpfr_set_d(xx, x, 0);
  mpfr_set_d(yy, y, 0);
  if ((j=mpfr_cmp2(xx, yy)) != i) {
    printf("Error in mpfr_cmp2: x=%1.16e y=%1.16e mpfr_cmp2(x,y)=%d instead of %d\n",x,y,j,i); 
    exit(1);
  }
  mpfr_clear(xx); mpfr_clear(yy);
}

int main()
{
  int i,j; double x=1.0, y, z;
#ifdef IRIX64
    /* to get denormalized numbers on IRIX64 */
    union fpc_csr exp;
    exp.fc_word = get_fpc_csr();
    exp.fc_struct.flush = 0;
    set_fpc_csr(exp.fc_word);
#endif

  tcmp2(1.06022698059744327881e+71, 1.05824655795525779205e+71, -1);
  tcmp2(1.0, 1.0, 53);
  for (i=0;i<54;i++) {
    tcmp2(1.0, 1.0-x, i);
    x /= 2.0;
  }
  for (j=0;j<1000000;j++) {
    x = drand(); if (x<0) x = -x;
    y = drand(); if (y<0) y = -y;
    if (!isnan(x) && !isnan(y)) {
      if (x<y) { z=x; x=y; y=z; }
      tcmp2(x, y, -1);
    }
  }
  return 0;
}

