#include <stdio.h>
#include "gmp.h"
#include "mpfr.h"

/* sets x to y^n */
void 
#if __STDC__
mpfr_pow_ui (mpfr_ptr x, mpfr_srcptr y, unsigned int n, unsigned char rnd)
#else
mpfr_pow_ui (x, y, n, rnd)
     mpfr_ptr x;
     mpfr_srcptr y; 
     unsigned int n; 
     unsigned char rnd; 
#endif
{
  int i;
  
  if (n==0) { mpfr_set_ui(x, 1, rnd); return; }
  mpfr_set(x, y, rnd);
  for (i=0;(1<<i)<=n;i++);
  /* now 2^(i-1) <= n < 2^i */
  for (i-=2; i>=0; i--) {
    mpfr_mul(x, x, x, rnd);
    if (n & (1<<i)) mpfr_mul(x, x, y, rnd);
  }
  return;
}

/* sets x to y^n */
void mpfr_ui_pow_ui (mpfr_ptr x, unsigned int y, unsigned int n, 
		     unsigned char rnd)
{
  int i;

  if (n==0) { mpfr_set_ui(x, 1, rnd); return; }
  mpfr_set_ui(x, y, rnd);
  for (i=0;(1<<i)<=n;i++);
  /* now 2^(i-1) <= n < 2^i */
  for (i-=2; i>=0; i--) {
    mpfr_mul(x, x, x, rnd);
    if (n & (1<<i)) mpfr_mul_ui(x, x, y, rnd);
  }
  return;
}
