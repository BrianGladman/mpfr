#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmp.h"
#include "longlong.h"
#include "mpfr.h"
#include "mpfr-impl.h"

void main()
{
  double x,y; mpfr_t xx,yy; int i,c;

  mpfr_init2(xx, 65); mpfr_init2(yy, 65);
  mpfr_set_str_raw(xx, "0.10011010101000110101010000000011001001001110001011101011111011101E623");
  mpfr_set_str_raw(yy, "0.10011010101000110101010000000011001001001110001011101011111011100E623");
  if (mpfr_cmp2(xx,yy)!=64) { printf("Error (1) in mpfr_cmp\n"); exit(1); }
  mpfr_set_str_raw(xx, "0.10100010001110110111000010001000010011111101000100011101000011100");
  mpfr_set_str_raw(yy, "0.10100010001110110111000010001000010011111101000100011101000011011");
  if (mpfr_cmp2(xx,yy)!=64) { printf("Error (1) in mpfr_cmp\n"); exit(1); }
  mpfr_set_prec(xx,53); mpfr_set_prec(yy,200);
  mpfr_set_d(xx, 1.0, 0);
  mpfr_set_d(yy, 1.0, 0);
  if (mpfr_cmp(xx,yy)!=0) {
    printf("Error in mpfr_cmp: 1.0 != 1.0\n"); exit(1);
  }
  mpfr_set_prec(yy, 31);
  mpfr_set_d(xx, 1.0000000002, 0);
  mpfr_set_d(yy, 1.0, 0);
  if (!(mpfr_cmp(xx,yy)>0)) {
    printf("Error in mpfr_cmp: not 1.0000000002 > 1.0\n"); exit(1);
  }
  mpfr_set_prec(yy, 53);
  for (i=0;i<1000000;) {    x=drand(); y=drand();
    if (!isnan(x) && !isnan(y)) {
      i++;
      mpfr_set_d(xx, x, 0);
      mpfr_set_d(yy, y, 0);
      c = mpfr_cmp(xx,yy);
      if ((c>0 && x<=y) || (c==0 && x!=y) || (c<0 && x>=y)) {
	printf("Error in mpfr_cmp with x=%1.20e, y=%1.20e mpfr_cmp(x,y)=%d\n",
	       x,y,c); exit(1);
      }
    }
  }
  mpfr_clear(xx); mpfr_clear(yy);
}
