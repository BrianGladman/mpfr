#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

/* #define DEBUG2 */

#ifdef DEBUG2
mp_limb_t *ap0, *ap1;
#define ck(x, dx, i) /* mp_limb_t *x; int dx; */ \
{ \
  if ((x)<ap0 || (x)+(dx)-1>ap1) { \
    printf("error %d\n",(i)); \
    printf("%1.20e,%d,%1.20e,%d,%d,%d\n",mpfr_get_d(b),PREC(b),mpfr_get_d(c), \
	   PREC(c),PREC(a),rnd_mode); \
    exit(1); } \
}
#else
#define ck(x, dx, i) {}
#endif

void mpfr_sub(a, b, c, rnd_mode) 
mpfr_ptr a; mpfr_srcptr b, c; unsigned char rnd_mode;
{
  CHANGE_SIGN(c);
  mpfr_add(a, b, c, rnd_mode);
  CHANGE_SIGN(c);
}

/* put in ap[0]..ap[an-1] the value of bp[0]..bp[n-1] shifted by sh bits
   to the left minus ap[0]..ap[n-1], with 0 <= sh < mp_bits_per_limb, and
   returns the borrow.
*/
mp_limb_t mpn_sub_lshift_n (ap, bp, n, sh, an) mp_limb_t *ap, *bp; int n,sh,an;
{
  mp_limb_t c, bh;

  /* shift b to the left */
  if (sh) bh = mpn_lshift(bp, bp, n, sh);
  c = (n<an) ? mpn_sub_n(ap, bp, ap, n) : mpn_sub_n(ap, bp+(n-an), ap, an);
  /* shift back b to the right */
  if (sh) {
    mpn_rshift(bp, bp, n, sh);
    bp[n-1] += bh<<(mp_bits_per_limb-sh);
  }
  return c;
}

/* signs of b and c differ, abs(b)>=abs(c), diff_exp>=0 */
void mpfr_sub1(a, b, c, rnd_mode, diff_exp) 
mpfr_t a, b, c; unsigned char rnd_mode; int diff_exp;
{
  mp_limb_t *ap, *bp, *cp, cc, c2; unsigned int an,bn,cn; 
  int sh,dif,k,cancel,cancel1,cancel2;

#ifdef DEBUG2
  printf("b=  "); if (SIGN(b)>=0) putchar(' ');
  mpfr_print_raw(b); putchar('\n');
  printf("c=  "); if (SIGN(c)>=0) putchar(' ');
  for (k=0;k<diff_exp;k++) putchar(' '); mpfr_print_raw(c);
  putchar('\n');
  printf("b=%1.20e c=%1.20e\n",mpfr_get_d(b),mpfr_get_d(c));
#endif
  cancel = mpfr_cmp2(b, c);
  /* we have to take into account the first (PREC(a)+cancel) bits from b */
  cancel1 = cancel/mp_bits_per_limb; cancel2 = cancel%mp_bits_per_limb;
  ap = MANT(a);
  bp = MANT(b);
  cp = MANT(c);
  an = (PREC(a)-1)/mp_bits_per_limb+1; /* number of significant limbs of a */
#ifdef DEBUG2
ap0=ap; ap1=ap+an-1;
#endif
  sh = an*mp_bits_per_limb-PREC(a); /* non-significant bits in low limb */
  bn = (PREC(b)-1)/mp_bits_per_limb+1; /* number of significant limbs of b */
  cn = (PREC(c)-1)/mp_bits_per_limb + 1;
  EXP(a) = EXP(b)-cancel;
  /* adjust sign to that of b */
  if (SIGN(a)*SIGN(b)<0) CHANGE_SIGN(a);
  /* case 1: diff_exp>=prec(a), i.e. c only affects the last bit
     through rounding */
  dif = PREC(a)-diff_exp;
#ifdef DEBUG2
  printf("PREC(a)=%d an=%u PREC(b)=%d bn=%u PREC(c)=%d diff_exp=%u dif=%d\n",
	 PREC(a),an,PREC(b),bn,PREC(c),diff_exp,dif);
#endif
  if (dif<=0) { /* diff_exp>=PREC(a): c does not overlap with a */
    /* either PREC(b)<=PREC(a), and we can copy the mantissa of b directly 
       into that of a, or PREC(b)>PREC(a) and we have to round b-c */
    if (PREC(b)<=PREC(a)+cancel) {
      if (cancel2) { ck(ap+(an-bn+cancel1), bn-cancel1, 1);
	mpn_lshift(ap+(an-bn+cancel1), bp, bn-cancel1, cancel2); }
      else { ck(ap+(an-bn+cancel1), bn-cancel1, 2);
	MPN_COPY(ap+(an-bn+cancel1), bp, bn-cancel1); }
      /* fill low significant limbs with zero */
ck(ap, an-bn+cancel1, 3);
      MPN_ZERO(ap, an-bn+cancel1);
      /* now take c into account */
      if (rnd_mode==GMP_RNDN) { /* to nearest */
	/* if diff_exp > PREC(a), no change */
	if (diff_exp==PREC(a)) {
	  /* if c is not zero, then as it is normalized, we have to subtract
	     one to the lsb of a if c>1/2, or c=1/2 and lsb(a)=1 (round to
	     even) */
	  if (NOTZERO(c)) { /* c is not zero */
	    /* check whether mant(c)=1/2 or not */
	    cc = *cp - (1<<(mp_bits_per_limb-1));
	    if (cc==0) {
	      bp = cp+(PREC(c)-1)/mp_bits_per_limb;
	      while (cp<bp && cc==0) cc = *++cp;
	    }
ck(ap+an-1, 1, 4);
	    if (cc || (ap[an-1] & 1<<sh)) goto sub_one_ulp;
	      /* mant(c) > 1/2 or mant(c) = 1/2: subtract 1 iff lsb(a)=1 */
	  }
	}
      }
      else if ((ISNONNEG(b) && rnd_mode==GMP_RNDU) || 
	       (ISNEG(b) && rnd_mode==GMP_RNDD)) {
	/* round up: nothing to do */
      }
      else {
	/* round down: subtract 1 ulp iff c<>0 */
	if (NOTZERO(c)) goto sub_one_ulp;
      }
    }
    else { /* PREC(b)>PREC(a) : we have to round b-c */
      k=bn-an;
      /* first copy the 'an' most significant limbs of b to a */
ck(ap, an, 5);
      MPN_COPY(ap, bp+k, an);
      if (rnd_mode==GMP_RNDN) { /* to nearest */
	/* first check whether the truncated bits from b are 1/2*lsb(a) */
	if (sh) {
ck(ap, 1, 6);
	  cc = *ap & (((mp_limb_t)1<<sh)-1);
	  *ap &= ~cc; /* truncate last bits */
	  cc -= (mp_limb_t)1<<(sh-1);
	}
	else /* no bit to truncate */
	  cc = bp[--k] - (1<<(mp_bits_per_limb-1));
	if ((long)cc>0) { /* suppose sizeof(long)=sizeof(mp_limb_t) */
	  goto add_one_ulp; /* trunc(b)>1/2*lsb(a) -> round up */
	}
	else if (cc==0) {
	  while (k>1 && cc==0) cc=bp[--k];
	  /* now if the truncated part of b = 1/2*lsb(a), check whether c=0 */
ck(ap, 1, 7);
	  if (NOTZERO(c) || (*ap & (1<<sh))) goto sub_one_ulp;
	  /* if trunc(b)-c is exactly 1/2*lsb(a) : round to even lsb */
	}
	/* if cc<0 : trunc(b) < 1/2*lsb(a) -> round down, i.e. do nothing */
      }
      else { /* round towards infinity or zero */
	if (sh) {
ck(ap, 1, 8);
	  cc = *ap & ((1<<sh)-1);
	  *ap &= ~cc; /* truncate last bits */
	}
	else cc=0;
	cn--;
        c2 = (dif>-sh) ? cp[cn]>>(mp_bits_per_limb-dif-sh) : 0;
	while (cc==c2 && (k || cn)) {
	  if (k) cc = bp[--k];
	  if (cn) {
	    c2 = cp[cn]<<(dif+sh);
	    if (--cn) c2 += cp[cn]>>(mp_bits_per_limb-dif-sh);
	  }
	}
	dif = ((ISNONNEG(b) && rnd_mode==GMP_RNDU) || 
	       (ISNEG(b) && rnd_mode==GMP_RNDD));
	/* round towards infinity if dif=1, towards zero otherwise */
	if (dif && cc>c2) goto add_one_ulp;
	else if (dif==0 && cc<c2) goto sub_one_ulp;
      }
    }
  }
  else { /* case 2: diff_exp < PREC(a) : c overlaps with a by dif bits */
    /* first copy upper part of c into a (after shift) */
    int overlap;
    dif += cancel;
    k = (dif-1)/mp_bits_per_limb + 1; /* only the highest k limbs from c
					 have to be considered */
    if (k<an) {
ck(ap+k, an-k, 9);
      MPN_ZERO(ap+k, an-k); /* do it now otherwise ap[k] may be 
				       destroyed in case dif<0 */
    }
#ifdef DEBUG2
    printf("cancel=%d dif=%d k=%d cn=%d sh=%d\n",cancel,dif,k,cn,sh);
#endif
    if (dif<=PREC(c)) { /* c has to be truncated */
      dif = dif % mp_bits_per_limb;
      dif = (dif) ? mp_bits_per_limb-dif-sh : -sh;
      /* we have to shift by dif bits to the right */
      if (dif>0) {
ck(ap, (k<=an) ? k : an, 10);
	mpn_rshift(ap, cp+(cn-k), (k<=an) ? k : an, dif);
        if (k>an) ap[an-1] += cp[cn-k+an]<<(mp_bits_per_limb-dif);
      }
      else if (dif<0) {
ck(ap, k, 11);
	cc = mpn_lshift(ap, cp+(cn-k), k, -dif);
	if (k<an) { ck(ap+k, 1, 12); ap[k]=cc; }
	/* put the non-significant bits in low limb for further rounding */
	ap[0] += cp[cn-k-1]>>(mp_bits_per_limb+dif);
      }
      else { ck(ap, k, 13); MPN_COPY(ap, cp+(cn-k), k); }
      overlap=1;
    }
    else { /* c is not truncated, but we have to fill low limbs with 0 */
ck(ap, k-cn, 14);
      MPN_ZERO(ap, k-cn);
      overlap = cancel-diff_exp;
#ifdef DEBUG2
      printf("0:a="); mpfr_print_raw(a); putchar('\n');
      printf("overlap=%d\n",overlap);
#endif
      if (overlap>=0) {
	cn -= overlap/mp_bits_per_limb;
	/* warning: a shift of zero with mpn_lshift is not allowed */
	if (overlap) {
	  if (an<cn) {
ck(ap, an, 15);
	    mpn_lshift(ap, cp+(cn-an), an, overlap%mp_bits_per_limb);
	    ap[0] += cp[cn-an-1]>>(mp_bits_per_limb-overlap);
	  }
	  else {
ck(ap+(an-cn), cn, 15);
	    mpn_lshift(ap+(an-cn), cp, cn, overlap%mp_bits_per_limb);
	  }
	}
	else { ck(ap+(an-cn), cn, 15); MPN_COPY(ap+(an-cn), cp, cn); }
      }
      else { /* shift to the right by -overlap bits */
	overlap = -overlap;
	k = overlap/mp_bits_per_limb;
	cc=mpn_rshift(ap+(an-k-cn),cp,cn,overlap%mp_bits_per_limb);
	if (an>k+cn) ap[an-k-cn-1]=cc;
      }
      overlap=0;
    }
#ifdef DEBUG2
      printf("1:a="); mpfr_print_raw(a); putchar('\n');
#endif
    /* here overlap=1 iff ulp(c)<ulp(a) */
    /* then put high limbs to zero */
    /* now add 'an' upper limbs of b in place */
    if (PREC(b)<=PREC(a)+cancel) { int i, imax;
      overlap += 2;
      /* invert low limbs */
      imax = (int)an-(int)bn+cancel1;
      for (i=0;i<imax;i++) { ck(ap+i, 1, 17); ap[i] = ~ap[i]; }
if (i) ck(ap, i, 18);
      cc = (i) ? mpn_add_1(ap, ap, i, 1) : 1;
ck(ap+i, ((bn-cancel1)<an) ? (bn-cancel1) : an, 19);
      mpn_sub_lshift_n(ap+i, bp, bn-cancel1, cancel2, an);
ck(ap+i, an-i, 20);
      mpn_sub_1(ap+i, ap+i, an-i, 1-cc);
    }
    else { /* PREC(b) > PREC(a): we have to truncate b */
ck(ap, an, 21);
      mpn_sub_lshift_n(ap, bp+(bn-an-cancel1), an, cancel2, an);
    }
    /* remains to do the rounding */
#ifdef DEBUG2
      printf("2:a="); mpfr_print_raw(a); putchar('\n');
      printf("overlap=%d\n",overlap);
#endif
    if (rnd_mode==GMP_RNDN) { /* to nearest */
      int kc;
      /* four cases: overlap =
         (0) PREC(b) > PREC(a) and diff_exp+PREC(c) <= PREC(a)
         (1) PREC(b) > PREC(a) and diff_exp+PREC(c) > PREC(a)
         (2) PREC(b) <= PREC(a) and diff_exp+PREC(c) <= PREC(a)
         (3)  PREC(b) <= PREC(a) and diff_exp+PREC(c) > PREC(a) */
      switch (overlap)
	{
        case 1: /* both b and c to round */
	  kc = cn-k; /* remains kc limbs from c */
	  k = bn-an; /* remains k limbs from b */
	  /* truncate last bits and store the difference with 1/2*ulp in cc */
	  cc = *ap & ((1<<sh)-1);
	  *ap &= ~cc; /* truncate last bits */
	  cc -= 1<<(sh-1);
	  while ((cc==0 || cc==-1) && k!=0 && kc!=0) {
	    kc--;
	    cc -= mpn_sub_1(&c2, bp+(--k), 1, (cp[kc]>>dif) +
			    (cp[kc+1]<<(mp_bits_per_limb-dif)));
	    if (cc==0 || cc==-1) cc=c2;
	  }
	  if ((long)cc>0) goto add_one_ulp;
	  else if ((long)cc<-1) return; /* the carry can be at most 1 */
	  else if (kc==0) goto round_b;
	  /* else round c: go through */
	case 3: /* only c to round */
	  bp=cp; k=cn-k; goto to_nearest;
	case 0: /* only b to round */
        round_b:
	  k=bn-an; dif=0; goto to_nearest;
        /* otherwise the result is exact: nothing to do */
	}
    }
    else if ((ISNONNEG(b) && rnd_mode==GMP_RNDU) || 
	     (ISNEG(b) && rnd_mode==GMP_RNDD)) {
      cc = *ap & ((1<<sh)-1);
      *ap &= ~cc; /* truncate last bits */
      if (cc) goto add_one_ulp; /* will happen most of the time */
      else { /* same four cases too */
	int kc = cn-k; /* remains kc limbs from c */
	switch (overlap)
	{
        case 1: /* both b and c to round */
	  k = bn-an; /* remains k limbs from b */
          dif = diff_exp % mp_bits_per_limb;
	  while (cc==0 && k!=0 && kc!=0) {
	    kc--;
	    cc = bp[--k] - (cp[kc]>>dif);
	    if (dif) cc -= (cp[kc+1]<<(mp_bits_per_limb-dif));
	  }
	  if (cc) goto add_one_ulp;
	  else if (kc==0) goto round_b2;
	  /* else round c: go through */
	case 3: /* only c to round: nothing to do */
	  /* while (kc) if (cp[--kc]) goto add_one_ulp; */
	  /* if dif>0 : remains to check last dif bits from c */
	  /* if (dif>0 && (cp[0]<<(mp_bits_per_limb-dif))) goto add_one_ulp; */
	  break;
	case 0: /* only b to round */
        round_b2:
	  k=bn-an;
	  while (k) if (bp[--k]) goto add_one_ulp;
        /* otherwise the result is exact: nothing to do */
	}
      }
    }
    /* else round to zero: remove 1 ulp if neglected bits from b-c are < 0 */
    else {
      cc = *ap & ((1<<sh)-1); 
      *ap &= ~cc;
      if (cc==0) { /* otherwise neglected difference cannot be < 0 */
	/* take into account bp[0]..bp[bn-cancel1-1] shifted by cancel2 to left
	   and cp[0]..cp[cn-k-1] shifted by dif bits to right */
	switch (overlap) {
	case 0:
	case 2:
	  break; /* c is not truncated ==> no borrow */
	case 1: /* both b and c are truncated */
	  break;
	case 3: /* only c is truncated */
	  cn -= k; /* take into account cp[0]..cp[cn-1] shifted by dif bits
		      to the right */
	  cc = (dif>0) ? cp[cn]<<(mp_bits_per_limb-dif) : 0;
	  while (cc==0 && cn>0) cc = cp[--cn];
          if (cc) goto sub_one_ulp;
	  break;
	}
      }
    }
  }
  goto end_of_sub;
    
  to_nearest: /* 0 <= sh < mp_bits_per_limb : number of bits of a to truncate
                 bp[k] : last significant limb from b */
#ifdef DEBUG2
mpfr_print_raw(a); putchar('\n');
#endif
        if (sh) {
	  cc = *ap & ((1<<sh)-1);
	  *ap &= ~cc; /* truncate last bits */
	  c2 = 1<<(sh-1);
	}
	else /* no bit to truncate */
	  { cc = bp[--k]; c2 = 1<<(mp_bits_per_limb-1); }
#ifdef DEBUG2
	printf("cc=%lu c2=%lu k=%u\n",cc,c2,k);
#endif
	if (cc>c2) goto add_one_ulp; /* trunc(b)>1/2*lsb(a) -> round up */
	else if (cc==c2) {
	  cc=0; while (k && cc==0) cc=bp[--k];
#ifdef DEBUG2
	  printf("cc=%lu\n",cc);
#endif
	  /* special case of rouding c shifted to the right */
	  if (cc==0 && dif>0) cc=bp[0]<<(mp_bits_per_limb-dif);
	  /* now if the truncated part of b = 1/2*lsb(a), check whether c=0 */
          if (bp!=cp) { 
	    if (cc || (*ap & (1<<sh))) goto add_one_ulp;
	  }
	  else {
	    /* subtract: if cc>0, do nothing */
	    if (cc==0 && (*ap & (1<<sh))) goto add_one_ulp;
	  }
	}
        goto end_of_sub;

 sub_one_ulp:
    cc = mpn_sub_1(ap, ap, an, 1<<sh);
    if (cc) { printf("carry(2) in mpfr_sub1\n"); exit(1); }
  goto end_of_sub;

  add_one_ulp: /* add one unit in last place to a */
    cc = mpn_add_1(ap, ap, an, 1<<sh);
    if (cc) { printf("carry(3) in mpfr_sub1\n"); exit(1); }

 end_of_sub:
#ifdef DEBUG2
printf("b-c="); if (SIGN(a)>0) putchar(' '); mpfr_print_raw(a); putchar('\n');
#endif
  return;
}

