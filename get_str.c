#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

/*
  Convert op to a string in base 'base' with 'n' digits and writes the 
  mantissa in 'str', the exponent in 'expptr'.
  The format is 0.xxxxxxxxEyyyy.
  The result is rounded wrt 'rnd_mode'.
 */

/* #define DEBUG */

char *mpfr_get_str(char *str, char *expptr, int base, size_t n,
		  mpfr_srcptr op, unsigned char rnd_mode)
{
  double d; long e, f, q, i, neg, p, err, prec, sh; mpfr_t a, b; mpz_t bz;
  char *str0; unsigned char rnd1;

#ifdef DEBUG
  printf("op="); mpfr_print_raw(op); printf(" rnd_mode=%d\n",rnd_mode);
  printf("  =%1.20e\n",mpfr_get_d(op));
#endif
  /* first determines the exponent */
  e = EXP(op); 
  EXP(op)=0; d=fabs(mpfr_get_d(op)); EXP(op)=e;
  /* the absolute value of op is between 1/2*2^e and 2^e */
  /* the output exponent f is such that base^(f-1) <= |op| < base^f
     i.e. f = 1 + floor(log(|op|)/log(base))
     = 1 + floor((log(|m|)+e*log(2))/log(base)) */
  f = 1 + (int) floor((log(d)+(double)e*log(2.0))/log((double)base));
#ifdef DEBUG
  printf("exponent = %d\n",f);
#endif
  /* now the first n digits of the mantissa are obtained from
     rnd(op*base^(n-f)) */
  prec = (long) ceil((double)n*log((double)base)/log(2.0));
  err = 5;
  q = prec+err;
  /* one has to use at least q bits */
  q = ((q-1)/mp_bits_per_limb)*mp_bits_per_limb;
  mpfr_init(a); mpfr_init(b);
  p = n-f; if ((neg=(p<0))) p=-p;
  rnd1 = (neg) ? GMP_RNDU : GMP_RNDZ; /* if neg we divide by base^p */
  do {
    q += mp_bits_per_limb;
    /* compute base^p with q bits and rounding towards zero */
    mpfr_set_prec(b, q, GMP_RNDZ);
    if (p==0) {
      mpfr_set(b, op, GMP_RNDZ);
    } 
    else {
      mpfr_set_prec(a, q, rnd1); 
      mpfr_set_ui(a, base, rnd1);
      for (i=0;(1<<i)<=p;i++);
      /* now 2^(i-1) <= p < 2^i */
       for (i-=2; i>=0; i--) {
	 mpfr_mul(b, a, a, rnd1);
	 if (p & (1<<i)) mpfr_mul_ui(a, b, base, rnd1);
	 else mpfr_set(a, b, rnd1);
       }
       /* now a is an approximation by default of base^p */
       if (neg) mpfr_div(b, op, a, GMP_RNDZ);
       else mpfr_mul(b, op, a, GMP_RNDZ);
    }
    if (SIGN(op)<0) CHANGE_SIGN(b);
  } while (mpfr_can_round(b, q-err, GMP_RNDZ, rnd_mode, prec)==0);
  if (SIGN(op)<0)
    switch (rnd_mode) {
    case GMP_RNDU: rnd_mode=GMP_RNDZ; break;
    case GMP_RNDD: rnd_mode=GMP_RNDU; break;
  }
#ifdef DEBUG
printf("rnd=%d\n",rnd_mode);
printf("b="); mpfr_print_raw(b); putchar('\n');
printf("=%1.20e\n",mpfr_get_d(b));
#endif
  prec=EXP(b); 
  mpfr_round(b, rnd_mode, prec);
  prec=EXP(b); /* may have chnaged due to rounding */
#ifdef DEBUG
printf("b="); mpfr_print_raw(b); putchar('\n');
printf("prec=%d q=%d b=",prec,q); mpfr_print_raw(b); putchar('\n');
printf("=%1.20e\n",mpfr_get_d(b));
#endif
  /* now the mantissa is the integer part of b */
  mpz_init(bz); q=1+(prec-1)/mp_bits_per_limb; _mpz_realloc(bz, q);
  sh = prec%mp_bits_per_limb;
  if (sh) mpn_rshift(PTR(bz), MANT(b), q, mp_bits_per_limb-sh);
  else MPN_COPY(PTR(bz), MANT(b), q);
  bz->_mp_size=q;
#ifdef DEBUG
printf("bz="); mpz_out_str(stdout,10,bz); putchar('\n');
printf("b="); mpfr_print_raw(b); putchar('\n');
#endif
  /* computes the number of characters needed */
  q = ((SIGN(op)<0) ? 1 : 0) + 2 + n + 2 +
    + (int) ceil(log((double)fabs(f))/log(10.0));
  if (f<10 || f>-10) q++;
  if (str==NULL) str0=str=malloc(q);
  if (SIGN(op)<0) *str++='-';
  if (n>1) *str++ = '.';
  mpz_get_str(str, base, bz); /* n digits of mantissa */
  if (strlen(str)==n+1) f++; /* possible due to rounding */
  str[n++] = 'e'; 
  f--; /* replaces 0.xxx*b^f by x.xx*b^(f-1) */
  str[n++] = (f>=0) ? '+' : '-'; /* is there a rule for f=0 ? */
  if (f<10 && f>-10) str[n++]='0';
  sprintf(str+n, "%ld", (f<0) ? -f : f);
  if (str[-1]=='.') { str[-1]=str[0]; str[0]='.'; }
  mpfr_clear(a); mpfr_clear(b); mpz_clear(bz);
  return str0;
}
