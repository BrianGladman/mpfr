#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"
#include <math.h>

/*
  Convert op to a string in base 'base' with 'n' digits and writes the 
  mantissa in 'str', the exponent in 'expptr'.
  The result is rounded wrt 'rnd_mode'.

  For op = 3.1416 we get str = "31416" and expptr=1.
 */

/* #define DEBUG */

char *mpfr_get_str(char *str, mp_exp_t *expptr, int base, size_t n,
		  mpfr_srcptr op, unsigned char rnd_mode)
{
  double d; long e, q, neg, p, err, prec, sh; mpfr_t a, b; mpz_t bz;
  char *str0; unsigned char rnd1; int f, pow2;

  if (base<2 || 36<base) {
    fprintf(stderr, "Error: too small or too large base in mpfr_get_str: %d\n",
	    base);
    exit(1);
  }
  count_leading_zeros(pow2, (mp_limb_t)base); 
  pow2 = mp_bits_per_limb - pow2 - 1;
  if (base != (1<<pow2)) pow2=0; 
  /* if pow2 <> 0, then base = 2^pow2 */

#ifdef DEBUG
  printf("op="); mpfr_print_raw(op); printf(" rnd_mode=%d\n",rnd_mode);
  printf("  =%1.20e\n",mpfr_get_d(op));
#endif
  /* first determines the exponent */
  e = EXP(op); 
  d = fabs(mpfr_get_d2(op, 0));
  /* the absolute value of op is between 1/2*2^e and 2^e */
  /* the output exponent f is such that base^(f-1) <= |op| < base^f
     i.e. f = 1 + floor(log(|op|)/log(base))
     = 1 + floor((log(|m|)+e*log(2))/log(base)) */
  f = 1 + (int) floor((log(d)+((double)e)*log(2.0))/log((double)base));
#ifdef DEBUG
  printf("exponent = %d\n",f);
#endif
  if (n==0)
    n = (int) ceil((double)PREC(op)*log(2.0)/log((double)base)+
		   log(4.0*fabs((double)((f==0) ? 1 : f)))/log(2.0));
  /* now the first n digits of the mantissa are obtained from
     rnd(op*base^(n-f)) */
  prec = (long) ceil((double)n*log((double)base)/log(2.0));
  err = 5;
  q = prec+err;
  /* one has to use at least q bits */
  q = ((q-1)/mp_bits_per_limb)*mp_bits_per_limb;
  mpfr_init(a); mpfr_init(b);
  p = n-f; if ((neg=(p<0))) p=-p;
#ifdef DEBUG
  printf("n=%d prec=%d p=%d\n",n,prec,p);
#endif
  rnd1 = rnd_mode;
  if (neg) {
    /* if neg we divide by base^p so we have to invert the rounding mode */
    switch (rnd1) {
    case GMP_RNDN: rnd1=GMP_RNDN; break;
    case GMP_RNDZ: rnd1=GMP_RNDU; break;
    case GMP_RNDU: rnd1=GMP_RNDZ; break;
    case GMP_RNDD: rnd1=GMP_RNDZ; break;
    }
  }
  do {
    q += mp_bits_per_limb;
    if (pow2) {
      if (neg) mpfr_div_2exp(b, op, pow2*p, rnd_mode);
      else mpfr_mul_2exp(b, op, pow2*p, rnd_mode);
    } 
    else {
       /* compute base^p with q bits and rounding towards zero */
       mpfr_set_prec(b, q);
       if (p==0) mpfr_set(b, op, rnd_mode);
       else {
	 mpfr_set_prec(a, q); 
	 mpfr_ui_pow_ui(a, base, p, rnd1);
	 /* now a is an approximation by default of base^p */
	 if (neg) mpfr_div(b, op, a, rnd_mode);
	 else mpfr_mul(b, op, a, rnd_mode);
       }
    }
    if (SIGN(op)<0) CHANGE_SIGN(b);
    if (q>2*prec+mp_bits_per_limb) {
      printf("no convergence in mpfr_get_str\n"); exit(1);
    }
  } while (pow2==0 && mpfr_can_round(b, q-err, rnd_mode, rnd_mode, prec)==0);
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
  q = ((SIGN(op)<0) ? 1 : 0) + n + 1;
  if (str==NULL) str0=str=malloc(q); 
  if (SIGN(op)<0) *str++='-';
  /*  if (n>1) *str++ = '.'; */
  mpz_get_str(str, base, bz); /* n digits of mantissa */
  if (strlen(str)==n+1) f++; /* possible due to rounding */
  /*  str[n++] = 'e'; */
  /* str[n++] = (f>=0) ? '+' : '-'; */ /* is there a rule for f=0 ? */
  /*  if (f<10 && f>-10) str[n++]='0'; */
  *expptr = f;
  /*  if (str[-1]=='.') { str[-1]=str[0]; str[0]='.'; } */
  mpfr_clear(a); mpfr_clear(b); mpz_clear(bz);
  return str0;
}

