#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

/* set f to the integer z */
int mpfr_set_z (f, z, rnd) mpfr_ptr f; mpz_srcptr z; 
unsigned char rnd;
{
  int fn, zn, k, dif, sign_z, sh; mp_limb_t *fp = MANT(f), cc, c2;

  sign_z = mpz_cmp_ui(z,0);
  if (sign_z==0) return (SIZE(f)=0);
  fn = 1 + (PREC(f)-1)/mp_bits_per_limb;
  zn = SIZ(z);
  dif = zn-fn;
  count_leading_zeros(k, PTR(z)[zn-1]);
  EXP(f) = zn*mp_bits_per_limb-k;
  if (SIGN(f)*sign_z<0) CHANGE_SIGN(f);
  if (dif>=0) { /* number has to be truncated */
    if (k) {
      mpn_lshift(fp, PTR(z) + dif, fn, k);
      if (dif) *fp += PTR(z)[dif-1] >> (mp_bits_per_limb-k);
    }
    else MPN_COPY(fp, PTR(z) + dif, fn);
    sh = fn*mp_bits_per_limb-PREC(f);
    cc = *fp & (((mp_limb_t)1 << sh) - 1);
    *fp = *fp & ~cc;
    if (rnd==GMP_RNDN) {
      if (sh) c2 = (mp_limb_t)1 << (sh-1);
      else { /* sh=0 */
	c2 = (mp_limb_t)1 << (mp_bits_per_limb-1);
	dif--;
	cc = (dif>=0) ? PTR(z)[dif] : 0;
      }
      /* now compares cc to c2 */
      if (cc>c2) { mpfr_add_one_ulp(f); return cc; }
      else if (cc<c2) goto towards_zero;
      else {
	cc=0;
	while (dif>0 && (cc=PTR(z)[dif-1])==0) dif--;
	if (cc) { mpfr_add_one_ulp(f); return cc; }
	else /* exactly in middle: inexact in both cases */
	  if (*fp & ((mp_limb_t)1<<sh)) { mpfr_add_one_ulp(f); return 1; }
	  else return 1;
      }
    }
    else if ((sign_z>0 && rnd==GMP_RNDU) || (sign_z<0 && rnd==GMP_RNDD)) {
      /* round towards infinity */
      /* result is exact iff all remaining bits are zero */
      if (dif>0 && cc==0) cc=PTR(z)[--dif]<<k;
      while (cc==0 && dif>0) cc=PTR(z)[--dif];
      if (cc) { mpfr_add_one_ulp(f); return 1; }
      else return 0;
    }
    else { /* round towards zero */
      /* result is exact iff all remaining bits are zero */
    towards_zero:
      if (cc==0 && dif>0) cc=PTR(z)[--dif]<<k;
      while (cc==0 && dif>0) cc=PTR(z)[--dif];
      return cc;
    }
  }
  else {
    if (k) mpn_lshift(fp-dif, PTR(z), zn, k);
    else MPN_COPY(fp-dif, PTR(z), zn);
    /* fill with zeroes */
    MPN_ZERO(fp, -dif);
    return 0; /* result is exact */
  }
}



