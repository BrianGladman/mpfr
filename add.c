/* (c) PolKA project, Inria Lorraine */
/* written by Paul Zimmermann, February 1999 */
/* to do: - destination may be identical to operands */


#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

extern void mpfr_sub1 _PROTO((mpfr_ptr, mpfr_srcptr, mpfr_srcptr, 
			      unsigned char, int));

/* #define DEBUG2 */

#ifdef DEBUG2
mp_limb_t *ap0, *ap1;
#define ck(x, dx, i) /* mp_limb_t *x; int dx; */ \
{ \
  if (x<ap0 || x+dx-1>ap1) { \
    printf("error %d\n",i); \
    printf("%1.20e,%d,%1.20e,%d,%d,%d\n",mpfr_get_d(b),PREC(b),mpfr_get_d(c), \
	   PREC(c),PREC(a),rnd_mode); \
    exit(1); } \
}
#else
#define ck(x, dx, i) {}
#endif

/* signs of b and c are supposed equal,
   diff_exp is the difference between the exponents of b and c,
   which is supposed >= 0 */

void 
#if __STDC__
mpfr_add1(mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, 
	  unsigned char rnd_mode, int diff_exp) 
#else
mpfr_add1(a, b, c, rnd_mode, diff_exp) 
     mpfr_ptr a; 
     mpfr_srcptr b; 
     mpfr_srcptr c; 
     unsigned char rnd_mode; 
     int diff_exp;
#endif
{
  mp_limb_t *ap, *bp, *cp, cc, c2, c3=0; unsigned int an,bn,cn; int sh,dif,k;
  TMP_DECL(marker); 

#ifdef DEBUG2
  printf("b=  "); mpfr_print_raw(b); putchar('\n');
  printf(" =  %1.20e\n",mpfr_get_d(b));
  printf("c=  "); for (k=0;k<diff_exp;k++) putchar(' '); mpfr_print_raw(c);
  putchar('\n');
#endif
  
  TMP_MARK(marker); 
  ap = MANT(a);
  bp = MANT(b);
  cp = MANT(c); 
  if (ap == bp) {
    bp = (mp_ptr) TMP_ALLOC(ABSSIZE(b) * BYTES_PER_MP_LIMB); 
    MPN_COPY (bp, ap, ABSSIZE(b));
    if (ap == cp) { cp = bp; }
  }
  else if (ap == cp)
    {
      cp = (mp_ptr) TMP_ALLOC (ABSSIZE(c) * BYTES_PER_MP_LIMB);
      MPN_COPY(cp, ap, ABSSIZE(c)); 
    }

  an = (PREC(a)-1)/mp_bits_per_limb+1; /* number of significant limbs of a */

#ifdef DEBUG2
  ap0=ap; ap1=ap+an-1;
#endif

  sh = an*mp_bits_per_limb-PREC(a); /* non-significant bits in low limb */
  bn = (PREC(b)-1)/mp_bits_per_limb+1; /* number of significant limbs of b */
  EXP(a) = EXP(b);

  if (SIGN(a)!=SIGN(b)) CHANGE_SIGN(a);

  /* case 1: diff_exp>=prec(a), i.e. c only affects the last bit
     through rounding */
  dif = PREC(a)-diff_exp;

#ifdef DEBUG2
  if (bn & ((mp_limb_t)1<<(mp_bits_per_limb-1))==0) {
    printf("Error in mpfr_add: b not msb-normalized\n"); exit(1);
  }
printf("an=%u prec(a)=%d bn=%u prec(b)=%d prec(c)=%d diff_exp=%u dif=%d\n",
       an,PREC(a),bn,PREC(b),PREC(c),diff_exp,dif);
#endif

  if (dif<=0) { 
    
    /* diff_exp>=PREC(a): c does not overlap with a */
    /* either PREC(b)<=PREC(a), and we can copy the mantissa of b directly 
       into that of a, or PREC(b)>PREC(a) and we have to round b+c */

    if (PREC(b)<=PREC(a)) {

      ck(ap+an-bn,bn,1);

      MPN_COPY(ap+(an-bn), bp, bn);
      /* fill low significant limbs with zero */

      for (bp=ap;bn<an;bn++) { ck(bp,1,2); *bp++=0; }

      /* now take c into account */
      if (rnd_mode==GMP_RNDN) { 
	
	/* to nearest */
	/* if diff_exp > PREC(a), no change */

	if (diff_exp==PREC(a)) {

	  /* if c is not zero, then as it is normalized, we have to add
	     one to the lsb of a if c>1/2, or c=1/2 and lsb(a)=1 (round to
	     even) */
	  
	  if (NOTZERO(c)) { 

	    /* c is not zero */
	    /* check whether mant(c)=1/2 or not */

	    cc = *cp - (1<<(mp_bits_per_limb-1));
	    if (cc==0) {
	      bp = cp+(PREC(c)-1)/mp_bits_per_limb;
	      while (cp<bp && cc==0) cc = *++cp;
	    }

	    ck(ap+an-1, 1, 3);
	    if (cc || (ap[an-1] & 1<<sh)) goto add_one_ulp;
	    /* mant(c) != 1/2 or mant(c) = 1/2: add 1 iff lsb(a)=1 */
	  }
	}
      }
      else if ((ISNONNEG(b) && rnd_mode==GMP_RNDU) || 
	       (ISNEG(b) && rnd_mode==GMP_RNDD)) {
	
	/* round up */
	if (NOTZERO(c)) goto add_one_ulp;
      }
      /* in the other cases (round to zero, or up/down with sign -/+),
         nothing to do */
    }
    else { 

      /* PREC(b)>PREC(a) : we have to round b+c */      
      k=bn-an;

      /* first copy the 'an' most significant limbs of b to a */
      ck(ap, an, 4);
      MPN_COPY(ap, bp+k, an);
      if (rnd_mode==GMP_RNDN) { 
	
	/* to nearest */
	/* first check whether the truncated bits from b are 1/2*lsb(a) */
	
	if (sh) {
	  ck(ap, 1, 5);
	  cc = *ap & ((1<<sh)-1);
	  *ap &= ~cc; /* truncate last bits */
	  cc -= 1<<(sh-1);
	}
	else /* no bit to truncate */
	  cc = bp[--k] - (1<<(mp_bits_per_limb-1));

	if ((long)cc>0) goto add_one_ulp; /* trunc(b)>1/2*lsb(a) -> round up */
	else if (cc==0) {

	  while (k>1 && cc==0) cc=bp[--k];

	  /* now if the truncated part of b = 1/2*lsb(a), check whether c=0 */
	  ck(ap, 1, 6);
	  if (NOTZERO(c) || (*ap & (1<<sh))) goto add_one_ulp;
	  /* if trunc(b)+c is exactly 1/2*lsb(a) : round to even lsb */
	}

	/* if cc<0 : trunc(b) < 1/2*lsb(a) -> round down, i.e. do nothing */
      }
      else if ((ISNONNEG(b) && rnd_mode==GMP_RNDU) || 
	       (ISNEG(b) && rnd_mode==GMP_RNDD)) {

	/* first check whether trunc(b)+c is zero or not */
	if (sh) {
	  ck(ap, 1, 7);
	  cc = *ap & ((1<<sh)-1); *ap &= ~cc; /* truncate last bits */
	}
	else cc = bp[--k] - (1<<(mp_bits_per_limb-1));
	while (cc==0 && k>1) cc=bp[--k];
	if (cc || NOTZERO(c)) goto add_one_ulp;
      }

      /* in the other cases (round to zero, or up/down with sign -/+),
         nothing to do, since b and c don't overlap, there can't be any
	 carry */

    }
  }
  else { 
    /* diff_exp < PREC(a) : c overlaps with a by dif bits */
    /* first copy upper part of c into a (after shift) */
    unsigned char overlap;
    
    k = (dif-1)/mp_bits_per_limb + 1; /* only the highest k limbs from c
					 have to be considered */
    cn = (PREC(c)-1)/mp_bits_per_limb + 1;
ck(ap+k, an-k, 8);
    MPN_ZERO(ap+k, an-k); /* do it now otherwise ap[k] may be destroyed
			     in case dif<0 */

#ifdef DEBUG2
    printf("k=%d cn=%d\n",k,cn);
#endif

    if (dif<=PREC(c)) { 
      /* c has to be truncated */
      dif = dif % mp_bits_per_limb;
      dif = (dif) ? mp_bits_per_limb-dif-sh : -sh;

      /* we have to shift by dif bits to the right */

      if (dif>0) { ck(ap, k, 9); mpn_rshift(ap, cp+(cn-k), k, dif); }
      else if (dif<0) {
	ck(ap, k+1, 10);
	ap[k] = mpn_lshift(ap, cp+(cn-k), k, -dif);

	/* put the non-significant bits in low limb for further rounding */

	if (cn >= k+1)
	  ap[0] += cp[cn-k-1]>>(mp_bits_per_limb+dif);
      }
      else { ck(ap, k, 11); MPN_COPY(ap, cp+(cn-k), k); }
      overlap=1;
    }
    else { 

      /* c is not truncated, but we have to fill low limbs with 0 */

      k = diff_exp/mp_bits_per_limb;
      overlap = diff_exp%mp_bits_per_limb;

      /* warning: a shift of zero bit is not allowed */
      MPN_ZERO(ap, an-k-cn); 
      if (overlap) { 
	cc=mpn_rshift(ap+(an-k-cn), cp, cn, overlap); 
	if (an-k-cn>0) ap[an-k-cn-1]=cc;
      }
      else { ck(ap+(an-k-cn), cn, 13); MPN_COPY(ap+(an-k-cn), cp, cn); }
      overlap=0;
    }

    /* here overlap=1 iff ulp(c)<ulp(a) */
    /* then put high limbs to zero */
    /* now add 'an' upper limbs of b in place */

    if (PREC(b)<=PREC(a)) {
      overlap += 2;
      ck(ap+(an-bn), bn, 14);
      cc = mpn_add_n(ap+(an-bn), ap+(an-bn), bp, bn);
    }
    else { 
      
      /* PREC(b) > PREC(a): we have to truncate b */
      
      ck(ap, an, 15);      
      cc = mpn_add_n(ap, ap, bp+(bn-an), an); 
    }

    if (cc) { 

      /* shift one bit to the right */

      ck(ap, an, 16);
      c3 = (ap[0]&1) && (PREC(a)%mp_bits_per_limb==0);
      mpn_rshift(ap, ap, an, 1);
      ap[an-1] += (mp_limb_t)1<<(mp_bits_per_limb-1);
      EXP(a)++;
    }
    
    /* remains to do the rounding */

    if (rnd_mode==GMP_RNDN) { 

      /* to nearest */

      int kc;
      
      /* four cases: overlap =
         (0) PREC(b) > PREC(a) and diff_exp+PREC(c) <= PREC(a)
         (1) PREC(b) > PREC(a) and diff_exp+PREC(c) > PREC(a)
         (2) PREC(b) <= PREC(a) and diff_exp+PREC(c) <= PREC(a)
         (3)  PREC(b) <= PREC(a) and diff_exp+PREC(c) > PREC(a) */
#ifdef DEBUG2
printf("overlap=%d\n",overlap);
#endif
      
      switch (overlap)
	{
        case 1: /* both b and c to round */
	  kc = cn-k; /* remains kc limbs from c */
	  k = bn-an; /* remains k limbs from b */
	
	  /* truncate last bits and store the difference with 1/2*ulp in cc */

	  ck(ap, 1, 17);
	  cc = *ap & ((1<<sh)-1);
	  *ap &= ~cc; /* truncate last bits */
	  cc -= 1<<(sh-1);
	  while ((cc==0 || cc==-1) && k!=0 && kc!=0) {
	    kc--;
	    cc += mpn_add_1(&c2, bp+(--k), 1,(cp[kc+1]<<(mp_bits_per_limb-dif))
			    +(cp[kc]>>dif));
	    if (cc==0 || cc==-1) cc=c2;
	  }
	  if ((long)cc>0) goto add_one_ulp;
	  else if ((long)cc<-1) 
	    { TMP_FREE(marker); return; /* the carry can be at most 1 */ }
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
      ck(ap, 1, 18);
      cc = *ap & ((1<<sh)-1);
      *ap &= ~cc; /* truncate last bits */
      if (cc || c3) goto add_one_ulp; /* will happen most of the time */
      else { 
	
	/* same four cases too */
	
	int kc = cn-k; /* remains kc limbs from c */
	switch (overlap)
	{
        case 1: /* both b and c to round */
	  k = bn-an; /* remains k limbs from b */
	  while (cc==0 && k!=0 && kc!=0) {
	    kc--;
	    cc = mpn_add_1(&c2, bp+(--k), 1,(cp[kc+1]<<(mp_bits_per_limb-dif))
			   + (cp[kc]>>dif));
	  }
	  if (cc) goto add_one_ulp;
	  else if (kc==0) goto round_b2;
	  /* else round c: go through */
	case 3: /* only c to round */
	  while (kc) if (cp[--kc]) goto add_one_ulp;
	  /* if dif>0 : remains to check last dif bits from c */
	  if (dif>0 && (cp[0]<<(mp_bits_per_limb-dif))) goto add_one_ulp;
	  break;
	case 0: /* only b to round */
        round_b2:
	  k=bn-an;
	  while (k) if (bp[--k]) goto add_one_ulp;
        /* otherwise the result is exact: nothing to do */
	}
      }
    }
    /* else nothing to do: round towards zero, i.e. truncate last sh bits */
    else {
ck(ap, 1, 19);
      *ap &= ~((1<<sh)-1);
    }
  }
  goto end_of_add;
    
  to_nearest: /* 0 <= sh < mp_bits_per_limb : number of bits of a to truncate
                 bp[k] : last significant limb from b */
        if (sh) {
ck(ap, 1, 20);
	  cc = *ap & ((1<<sh)-1);
	  *ap &= ~cc; /* truncate last bits */
	  c2 = 1<<(sh-1);
	}
	else /* no bit to truncate */
	  { if (k) cc = bp[--k]; else cc = 0; c2 = 1<<(mp_bits_per_limb-1); }
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
	  if (cc==0 && dif>0) cc=cp[0]<<(mp_bits_per_limb-dif);
	  /* now if the truncated part of b = 1/2*lsb(a), check whether c=0 */
ck(ap, 1, 21);
	  if (cc || (*ap & (1<<sh))) goto add_one_ulp;
	}
        goto end_of_add;

  add_one_ulp: /* add one unit in last place to a */
ck(ap, an, 22);
    cc = mpn_add_1(ap, ap, an, 1<<sh);
    if (cc) { printf("carry(3) in mpfr_add\n"); exit(1); }

 end_of_add:
#ifdef DEBUG2
printf("b+c="); mpfr_print_raw(a); putchar('\n');
#endif
  TMP_FREE(marker); 
  return;
}

void
#if __STDC__
mpfr_add(mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, 
	      unsigned char rnd_mode) 
#else
mpfr_add(a, b, c, rnd_mode)
     mpfr_ptr a; 
     mpfr_srcptr b; 
     mpfr_srcptr c; 
     unsigned char rnd_mode; 
#endif
{
  int diff_exp;

  if (FLAG_NAN(b) || FLAG_NAN(c)) {
    SET_NAN(a); return;
  }

  if (!NOTZERO(b)) { mpfr_set(a, c, rnd_mode); return; }
  if (!NOTZERO(c)) { mpfr_set(a, b, rnd_mode); return; }

  diff_exp = EXP(b)-EXP(c);
  if (SIGN(b) != SIGN(c)) { /* signs differ, it's a subtraction */
    if (diff_exp<0) {
      mpfr_sub1(a, c, b, rnd_mode, -diff_exp);
    }
    else if (diff_exp>0) mpfr_sub1(a, b, c, rnd_mode, diff_exp);
    else { /* diff_exp=0 */
      CHANGE_SIGN(c);
      diff_exp = mpfr_cmp(b,c);
      CHANGE_SIGN(c);
      /* if b>0 and diff_exp>0 or b<0 and diff_exp<0: abs(b) > abs(c) */
      if (diff_exp==0) SET_ZERO(a);
      else if (diff_exp*SIGN(b)>0) mpfr_sub1(a, b, c, rnd_mode, 0);
      else mpfr_sub1(a, c, b, rnd_mode, 0);
    }
  }
  else /* signs are equal, it's an addition */
    if (diff_exp<0) mpfr_add1(a, c, b, rnd_mode, -diff_exp);
    else mpfr_add1(a, b, c, rnd_mode, diff_exp);
}

