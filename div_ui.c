#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

/* #define DEBUG */

/* returns 0 if result exact, non-zero otherwise */
int
#ifdef __STDC__
mpfr_div_ui(mpfr_ptr y, mpfr_srcptr x, unsigned long u, unsigned char rnd_mode)
#else
mpfr_div_ui(y, x, u, rnd_mode)
     mpfr_ptr y;  
     mpfr_srcptr x;
     unsigned long u;
     unsigned char rnd_mode; 
#endif
{
  int xn, yn, dif, sh, i; mp_limb_t *xp, *yp, c, d;

  if (FLAG_NAN(x)) { SET_NAN(y); return 1; }
  if (u==0) { printf("infinity\n"); return 1; }

  xn = (PREC(x)-1)/BITS_PER_MP_LIMB + 1;
  yn = (PREC(y)-1)/BITS_PER_MP_LIMB + 1;

  xp = MANT(x);
  yp = MANT(y);
  EXP(y) = EXP(x);
  if (SIGN(x)!=SIGN(y)) CHANGE_SIGN(y);

  /* save limb yp[-1] that will be used to store an extra limb of
     the quotient */
  d = yp[-1];
  dif = yn+1-xn;
#ifdef DEBUG
  printf("dif=%d u=%lu\n",dif,u);
  printf("x="); mpfr_print_raw(x); putchar('\n');
#endif
  if (dif>=0)
    c = mpn_divrem_1(yp-1, dif, xp, xn, u);
  else /* dif < 0 i.e. xn > yn */
    c = mpn_divrem_1(yp-1, 0, xp-dif, yn, u);
#ifdef DEBUG
printf("y="); mpfr_print_raw(y); putchar('\n');
#endif

  /* shift left to normalize */
  count_leading_zeros(sh, yp[yn-1]);
  if (sh) { mpn_lshift(yp-1, yp-1, yn+1, sh); EXP(y) -= sh; }
  yp[-1] = d; /* restore value of yp[-1] */
#ifdef DEBUG
printf("y="); mpfr_print_raw(y); putchar('\n');
#endif

  sh = yn*BITS_PER_MP_LIMB - PREC(y);
  /* it remains sh bits in less significant limb of y */

  d = *yp & (((mp_limb_t)1 << sh) - 1);
  *yp ^= d; /* set to zero lowest sh bits */

  if ((c | d)==0) {
    for (i=0; i<-dif && xp[i]==0; i++);
    if (i>=-dif) return 0; /* result is exact */
  }

  switch (rnd_mode) {
  case GMP_RNDZ:
    return 1; /* result is inexact */
  case GMP_RNDU:
    if (SIGN(y)>0) mpfr_add_one_ulp(y);
    return 1; /* result is inexact */
  case GMP_RNDD:
    if (SIGN(y)<0) mpfr_add_one_ulp(y);
    return 1; /* result is inexact */
  case GMP_RNDN:
    if (d < ((mp_limb_t)1 << (sh-1))) return 1;
    else if (d > ((mp_limb_t)1 << (sh-1))) {
      mpfr_add_one_ulp(y);
    }
    else { /* d = (mp_limb_t)1 << (sh-1) */
      if (c) mpfr_add_one_ulp(y);
      else {
	for (i=0; i<-dif && xp[i]==0; i++);
	if (i<-dif) mpfr_add_one_ulp(y);
	else { /* exactly in the middle */
	  if (*yp & ((mp_limb_t)1 << sh)) mpfr_add_one_ulp(y);
	}
      }
    }
    return 1;
  }
  return 0; /* to prevent warning from gcc */
}
