/* (c) PolKA project at Inria Lorraine

written by Paul Zimmermann, February 1999

    returns 0 iff b = c
            a positive value iff b > c
            a negative value iff b < c

More precisely, in case b and c are of same sign, the absolute value 
of the result is one plus the absolute difference between the exponents 
of b and c, i.e. one plus the number of bits shifts to align b and c
(this value is useful in mpfr_sub).

*/

#include <stdio.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

/* #define DEBUG */

int 
#if __STDC__
mpfr_cmp ( mpfr_srcptr b, mpfr_srcptr c )
#else
mpfr_cmp(b, c)
     mpfr_srcptr b;
     mpfr_srcptr c; 
#endif
{
   long int s, diff_exp;
   unsigned long bn, cn;
   mp_limb_t *bp, *cp;

   s = SIGN(b) * SIGN(c);
   if (s<0) return(SIGN(b));

   /* now signs are equal */

   diff_exp = EXP(b)-EXP(c);
   s = (SIGN(b)>0) ? 1 : -1;

   if (diff_exp>0) return(s*(1+diff_exp));
   else if (diff_exp<0) return(s*(-1+diff_exp));
   /* both signs and exponents are equal */

   bn = (PREC(b)-1)/mp_bits_per_limb+1;
   cn = (PREC(c)-1)/mp_bits_per_limb+1;
   bp = MANT(b); cp = MANT(c);

   while (bn && cn) {
     if (bp[--bn] != cp[--cn])
       return((bp[bn]>cp[cn]) ? s : -s);
   }

   if (bn) { while (bn) if (bp[--bn]) return(s); }
   else if (cn) while (cn) if (cp[--cn]) return(-s);

   return 0;
}

/* returns the number of cancelled bits when one subtracts abs(c) from abs(b). 
   Assumes b>=c, which implies EXP(b)>=EXP(c).
   if b=c, returns prec(b).
*/
int 
#if __STDC__
mpfr_cmp2 ( mpfr_srcptr b, mpfr_srcptr c )
#else
mpfr_cmp2(b, c)
     mpfr_srcptr b; 
     mpfr_srcptr c; 
#endif
{
  long int d, bn, cn, k;
  mp_limb_t *bp, *cp, t, u=0, cc=0;

#ifdef DEBUG
  printf("b="); mpfr_print_raw(b); putchar('\n');
  printf("c="); mpfr_print_raw(c); putchar('\n');
#endif  
  if (NOTZERO(c)==0) return 0;
  d = EXP(b)-EXP(c);
  k = 0; /* result can be d or d+1 if d>1, or >= d otherwise */
#ifdef DEBUG
  if (d<0) { printf("assumption EXP(b)<EXP(c) violated\n"); exit(1); }
#endif
  bn = (PREC(b)-1)/mp_bits_per_limb;
  cn = (PREC(c)-1)/mp_bits_per_limb;
  bp = MANT(b); cp = MANT(c);
  /* subtracts c from b from most significant to less significant limbs,
     and first determines first non zero limb difference */
  if (d)
    {
      cc = bp[bn--];
      if (d<mp_bits_per_limb)
	cc -= cp[cn]>>d; 
    }
  else { /* d=0 */
    while (bn>=0 && cn>=0 && (cc=(bp[bn--]-cp[cn--]))==0) {
      k+=mp_bits_per_limb;
    }

    if (cc==0) { /* either bn<0 or cn<0 */
      while (bn>=0 && (cc=bp[bn--])==0) k+=mp_bits_per_limb;
    }
    /* now bn<0 or cc<>0 */
    if (cc==0 && bn<0) return(PREC(b));
  }

  /* the first non-zero limb difference is cc, and the number
     of cancelled bits in the upper limbs is k */

  count_leading_zeros(u, cc);
  k += u;

  if (cc != (1<<(mp_bits_per_limb-u-1))) return k;
  /* now cc is an exact power of two */
  
  if (cc != 1) 
    /* We just need to compare the following limbs */
    /* until two of them differ. The result is either k or k+1. */
    {
      /* First flush all the unmatched limbs of b ; they all have to
	 be 0 in order for the process to go on */
      while (bn >= 0)
	{
	  if (cn < 0) { return k; }
	  t = bp[bn--]; 
	  if (d < mp_bits_per_limb)
	    {
	      if (d) 
		{ 
		  u = cp[cn--] << (mp_bits_per_limb - d); 
		  if (cn >= 0) u+=(cp[cn]>>d); 
		}
	      else u = cp[cn--]; 
	      
	      if (t > u || (t == u && cn < 0)) return k; 
	      if (t < u) return k+1; 
	    }
	  else
	    if (t) return k; else d -= mp_bits_per_limb; 
	}
	  
      /* bn < 0; if some limb of c is nonzero, return k+1, otherwise return k*/

      if (cn>=0 && (cp[cn--] << (mp_bits_per_limb - d))) { return k+1; }

      while (cn >= 0) 
	if (cp[cn--]) return k+1; 
      return k; 
    }

  /* cc = 1. Too bad. */

  if (d > mp_bits_per_limb) { return k; }

  while (bn >= 0)
    { 
      if (cn < 0) { return k; }

      if (d) 
	 { 
	   u = cp[cn--] << (mp_bits_per_limb - d); 
	   if (cn >= 0) u+=(cp[cn]>>d);
	 }
       else u = cp[cn--]; 
       
       if ((cc = (bp[bn--] | ~u)) != 0)
	 { count_leading_zeros(u, cc); return k + u; }
       else k += mp_bits_per_limb; 
    }

  if (cn >= 0)
    count_leading_zeros(cc, ~(cp[cn--] << (mp_bits_per_limb - d))); 
  else { cc = 0; }

  k += cc; 
  if (cc < d) return k;
  
  while (cn >= 0 && !~cp[cn--]) { k += mp_bits_per_limb; } 
  if (cn >= 0) { count_leading_zeros(cc, ~cp[cn--]); return (k + cc); }

  return k; /* We **need** that the nonsignificant limbs of c are set
	       to zero there */	       
}
