/* written by Paul Zimmermann, November 1998-January 1999 */

#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

/* #define DEBUG */ /* check possible errors */
/* #define DEBUG2 */ /* print informations */

/* avoids two divisions */
#define GETBIT(a, n) (a)[(n)>>LOG_MP_BITS_PER_LIMB]>>((n)&((1<<LOG_MP_BITS_PER_LIMB)-1))

/* and here comes the code. To do:
- take sign(b) and sign(c) into account, and write sign(a)
- does not use all bits of b and c when PREC(b)>PREC(a) or PREC(c)>PREC(a)
- use fast mpn_mul multiplication
*/

void mpfr_mul(a, b, c, rnd_mode) 
mpfr_ptr a; mpfr_srcptr b, c; unsigned char rnd_mode;
{
  unsigned int bn, cn, an, k; int sh; unsigned char b1;
  mp_limb_t *ap, *bp, *cp, cc;
  long int sign_product;
  TMP_DECL(marker); 

  /* the multiplication works as follows:
     1. multiply all significant limbs from the mantissa of b and c in a 
        temporary space (ap)
     2. rounds the result up to PREC(a) bits
     3. copies the result into the destination
  */

#ifdef DEBUG2
  printf("b="); mpfr_print_raw(b); putchar('\n');
  printf("c="); mpfr_print_raw(c); putchar('\n');
#endif
  sign_product = SIGN(b) * SIGN(c);
  bn = (PREC(b)-1)/mp_bits_per_limb+1; /* number of significant limbs of b */
  cn = (PREC(c)-1)/mp_bits_per_limb+1; /* number of significant limbs of c */
  k = bn+cn; /* effective nb of limbs used by b*c */
  TMP_MARK(marker); 
  ap = (mp_limb_t*) TMP_ALLOC(k*BYTES_PER_MP_LIMB);
  bp = MANT(b); cp = MANT(c);

  /* step 1: multiplies two mantissa */
  if (bn>=cn) mpn_mul(ap, bp, bn, cp, cn);
  else mpn_mul(ap, cp, cn, bp, bn);

  /* now ap[0]..ap[k-1] contains the product of both mantissa,
     with ap[k-1]>=2^(mp_bits_per_limb-2) */
  an = (PREC(a)-1)/mp_bits_per_limb+1; /* number of significant limbs of a */
  sh = k*mp_bits_per_limb-PREC(a); /* number of bits to truncate */
  b1 = GETBIT(ap+k-1,mp_bits_per_limb-1); /* msb from the product */
  sh = sh+b1-1; /* sh-- if b1=0 */
#ifdef DEBUG2
printf("an=%u bn=%u cn=%u b1=%u sh=%u bp[0]=%lu cp[0]=%lu\n",an,bn,cn,b1,sh,*bp,*cp);
#endif
  /* if sh<=0, then no rounding is needed */
  if (sh<=0) { /* copy most significant k*mp_bits_per_limb bits */
    bp = a->_mp_d; cp = bp+an-k;
    if (b1) {
      MPN_COPY(cp, ap, k); 
      EXP(a) = EXP(b) + EXP(c);
    }
    else { 
      mpn_lshift(cp, ap, k, 1); 
      EXP(a) = EXP(b) + EXP(c) - 1;
    }
    /* set to zero low bits */
    while (cp>bp) *--cp=0;
  }
  else { /* k*mp_bits_per_limb + (b1-1) > PREC(a) */
     /* sh>0 is the number of low significant bits to truncate.
        (i) down mode: just put them to zero
        (ii) nearest mode: add 2^(sh-1) and truncate
        (iii) up mode: add 2^sh-1 and truncate
     */
     if (rnd_mode==GMP_RNDN) { /* nearest mode */
       /* if middle of two representable numbers, we must
          set the lower bit to zero according to the IEEE norm */
       sh--; /* to keep one guard bit */
       cc=0; /* 0 iff middle */
       while (cc==0 && sh>=mp_bits_per_limb) {
	 cc = *ap++; k--; sh -= mp_bits_per_limb;
       }
       ap += (sh/mp_bits_per_limb);
       k -= (sh/mp_bits_per_limb);
       sh %= mp_bits_per_limb; /* now 0<=sh<mp_bits_per_limb */
       if (cc==0)
	 /* middle iff last sh bits are zero and bit (sh+1) is 1 */
	 cc = (*ap<<(mp_bits_per_limb-sh-1)) - ((mp_limb_t)1<<(mp_bits_per_limb-1));
       if (cc)
	 cc = mpn_add_1(ap, ap, k, (mp_limb_t)1<<sh);
       else /* put lower (sh+1) bits to zero */
	 *ap = *ap & ~(((mp_limb_t)2<<sh)-1);
       if (cc!=0) { /* may happen especially when PREC(a) is small */
	 mpn_rshift(ap, ap, k, 1);
	 ap[k-1] += (mp_limb_t)1<<(mp_bits_per_limb-1);
	 b1++; /* so that EXP(a) will be incremented by 1 below */
       }
       else sh++;
       /* 0 <= sh <= mp_bits_per_limb */
     }
     else if ((sign_product>=0 && rnd_mode==GMP_RNDU) ||
         (sign_product<0 && rnd_mode==GMP_RNDD)) { /* up mode */
       cc = mpn_sub_1(ap, ap, k, 1); /* subtract 1 */
       ap += (sh/mp_bits_per_limb);
       k -= (sh/mp_bits_per_limb);
       sh %= mp_bits_per_limb; /* now 0<=sh<mp_bits_per_limb */
       cc = mpn_add_1(ap, ap, k, (mp_limb_t)1<<sh)-cc;
       if (cc!=0) { /* may happen especially when PREC(a) is small */
	 mpn_rshift(ap, ap, k, 1);
	 ap[k-1] += (mp_limb_t)1<<(mp_bits_per_limb-1);
	 b1++; /* so that EXP(a) will be incremented by 1 below */
	 if (sh==0) { sh=mp_bits_per_limb; ap--; k++; }
	 sh--;
       }
     }
     else { /* down mode */
       ap += (sh/mp_bits_per_limb);
       k -= (sh/mp_bits_per_limb);
       sh %= mp_bits_per_limb;
     }
     /* now truncate last sh bits in ap[0] */
#ifdef DEBUG2
printf("sh=%u *ap=%u\n",sh,*ap);
#endif
     if (sh==mp_bits_per_limb) { ap++; k--; sh=0; }
     else *ap = *ap & (~(mp_limb_t)0 << sh);
#ifdef DEBUG2
printf("*ap=%u\n",*ap);
#endif
     /* step 3: copies the result into the destination */
     cp=a->_mp_d; 
     EXP(a) = EXP(b) + EXP(c) + b1 - 1;
     if (b1) {
#ifdef DEBUG
       if (k!=an) /* normally this should not happen */
	 { printf("k=%u != an=%u\n",k,an); exit(1); }
#endif
       MPN_COPY(cp, ap, k);
     }
     else { 
       if (k!=an) { /* only when k=an+1 and sh=mp_bits_per_limb-1 */
#ifdef DEBUG
	 if (k!=an+1 || sh!=mp_bits_per_limb-1) 
	   { printf("k=%u != an+1=%u or sh=%u != mp_bits_per_limb-1\n",k,an+1,
		    sh); exit(1); }
#endif
	 mpn_lshift(cp, ap+1, an, 1);
	 *cp += ap[0]>>sh;
       }
       else mpn_lshift(cp, ap, an, 1); 
     }
   }
  if (sign_product<0) CHANGE_SIGN(a);
  TMP_FREE(marker); 
#ifdef DEBUG2
printf("b*c="); mpfr_print_raw(a); putchar('\n');
#endif
}

     
