#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

/* set f to the integer z */
int 
#if __STDC__
mpfr_set_z (mpfr_ptr f, mpz_srcptr z, unsigned char rnd)
#else
mpfr_set_z (f, z, rnd) 
     mpfr_ptr f; 
     mpz_srcptr z; 
     unsigned char rnd;
#endif
{
  int fn, zn, k, dif, sign_z, sh; mp_limb_t *fp = MANT(f), *zp, cc, c2;

  sign_z = mpz_cmp_ui(z,0);
  if (sign_z==0) return (SIZE(f)=0);
  fn = 1 + (PREC(f)-1)/BITS_PER_MP_LIMB;
  zn = SIZ(z);
  dif = zn-fn;
  zp = PTR(z);
  count_leading_zeros(k, zp[zn-1]);
  EXP(f) = zn*BITS_PER_MP_LIMB-k;
  if (SIGN(f)*sign_z<0) CHANGE_SIGN(f);
  if (dif>=0) { /* number has to be truncated */
    if (k) {
      mpn_lshift(fp, zp + dif, fn, k);
      if (dif) *fp += zp[dif-1] >> (BITS_PER_MP_LIMB-k);
    }
    else MPN_COPY(fp, zp + dif, fn);
    sh = fn*BITS_PER_MP_LIMB-PREC(f);
    cc = *fp & (((mp_limb_t)1 << sh) - 1);
    *fp = *fp & ~cc;
    if (rnd==GMP_RNDN) {
      if (sh) c2 = (mp_limb_t)1 << (sh-1);
      else { /* sh=0 */
	c2 = ((mp_limb_t)1) << (BITS_PER_MP_LIMB-1);
	dif--;
	cc = (dif>=0) ? ((zp[dif])<<k) : 0;
	if (dif>0 && k) cc += zp[dif-1] >> (BITS_PER_MP_LIMB-k);
      }
      /* now compares cc to c2 */
      if (cc>c2) { mpfr_add_one_ulp(f); return cc; }
      else if (cc<c2) goto towards_zero;
      else {
	cc=0;
	while (dif>0 && (cc=zp[dif-1])==0) dif--;
	if (cc) { mpfr_add_one_ulp(f); return cc; }
	else /* exactly in middle: inexact in both cases */
	  if (*fp & ((mp_limb_t)1<<sh)) { mpfr_add_one_ulp(f); return 1; }
	  else return 1;
      }
    }
    else if ((sign_z>0 && rnd==GMP_RNDU) || (sign_z<0 && rnd==GMP_RNDD)) {
      /* round towards infinity */
      /* result is exact iff all remaining bits are zero */
      if (dif>0 && cc==0) cc=zp[--dif]<<k;
      while (cc==0 && dif>0) cc=zp[--dif];
      if (cc) { mpfr_add_one_ulp(f); return 1; }
      else return 0;
    }
    else { /* round towards zero */
      /* result is exact iff all remaining bits are zero */
    towards_zero:
      if (cc==0 && dif>0) cc=zp[--dif]<<k;
      while (cc==0 && dif>0) cc=zp[--dif];
      return cc;
    }
  }
  else {
    if (k) mpn_lshift(fp-dif, zp, zn, k);
    else MPN_COPY(fp-dif, zp, zn);
    /* fill with zeroes */
    MPN_ZERO(fp, -dif);
    return 0; /* result is exact */
  }
}



