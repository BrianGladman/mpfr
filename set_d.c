#include <math.h> /* for isnan */
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

#define NaN sqrt(-1) /* ensures a machine-independent NaN */

/* Included from gmp-2.0.2, patched to support denorms */

#ifdef XDEBUG
#undef _GMP_IEEE_FLOATS
#endif

#ifndef _GMP_IEEE_FLOATS
#define _GMP_IEEE_FLOATS 0
#endif

#define MP_BASE_AS_DOUBLE (2.0 * ((mp_limb_t) 1 << (BITS_PER_MP_LIMB - 1)))
int
#if __STDC__
__mpfr_extract_double (mp_ptr rp, double d, int e)
#else
__mpfr_extract_double (rp, d)
     mp_ptr rp;
     double d;
     int e;
#endif
     /* e=0 iff rp has only one limb */
{
  long exp;
  mp_limb_t manh, manl;

  /* BUGS

     1. Should handle Inf and NaN in IEEE specific code.
     2. Handle Inf and NaN also in default code, to avoid hangs.
     3. Generalize to handle all BITS_PER_MP_LIMB >= 32.
     4. This lits is incomplete and misspelled.
   */

  if (d == 0.0)
    {
      rp[0] = 0;
#if BITS_PER_MP_LIMB == 32
      if (e) rp[1] = 0;
#endif
      return 0;
    }

#if _GMP_IEEE_FLOATS
  {
    union ieee_double_extract x;
    x.d = d;

    exp = x.s.exp;
    if (exp) 
      {
#if BITS_PER_MP_LIMB == 64
	manl = (((mp_limb_t) 1 << 63)
		| ((mp_limb_t) x.s.manh << 43) | ((mp_limb_t) x.s.manl << 11));
#else
	manh = ((mp_limb_t) 1 << 31) | (x.s.manh << 11) | (x.s.manl >> 21);
	manl = x.s.manl << 11;      
#endif
      }
    else
      {
#if BITS_PER_MP_LIMB == 64
	manl = ((mp_limb_t) x.s.manh << 43) | ((mp_limb_t) x.s.manl << 11);
#else
    manh = (x.s.manh << 11) | (x.s.manl >> 21);
	manl = x.s.manl << 11;      
#endif
      }
  }
#else
  {
    /* Unknown (or known to be non-IEEE) double format.  */
    exp = 0;
    if (d >= 1.0)
      {
        if (d * 0.5 == d)
          abort ();

        while (d >= 32768.0)
          {
            d *= (1.0 / 65536.0);
            exp += 16;
          }
        while (d >= 1.0)
          {
            d *= 0.5;
            exp += 1;
          }
      }
    else if (d < 0.5)
      {
        while (d < (1.0 / 65536.0))
          {
            d *=  65536.0;
            exp -= 16;
          }
        while (d < 0.5)
          {
            d *= 2.0;
            exp -= 1;
          }
      }

    d *= MP_BASE_AS_DOUBLE;
#if BITS_PER_MP_LIMB == 64
    manl = d;
#else
    manh = d;
    manl = (d - manh) * MP_BASE_AS_DOUBLE;
#endif

    exp += 1022;
  }
#endif

  if (exp) exp = (unsigned) exp - 1022; else exp = -1021; 

#if BITS_PER_MP_LIMB == 64
      rp[0] = manl;
#else
      if (e) {
	rp[1] = manh;
	rp[0] = manl;
      }
      else {
	rp[0] = manh;
      }
#endif

  return exp;
}

/* End of part included from gmp-2.0.2 */
/* Part included from gmp temporary releases */
double
#if __STDC__
__mpfr_scale2 (double d, int exp)
#else
__mpfr_scale2 (d, exp)
     double d;
     int exp;
#endif
{
#if _GMP_IEEE_FLOATS
  {
    union ieee_double_extract x;
    x.d = d;
    exp += x.s.exp;
    x.s.exp = exp;
    if (exp >= 2047)
      {
        /* Return +-infinity */
        x.s.exp = 2047;
        x.s.manl = x.s.manh = 0;
      }
    else if (exp < 1)
      {
        x.s.exp = 1;            /* smallest exponent (biased) */
        /* Divide result by 2 until we have scaled it to the right IEEE
           denormalized number, but stop if it becomes zero.  */
        while (exp < 1 && x.d != 0)
          {
            x.d *= 0.5;
            exp++;
          }
      }
    return x.d;
  }
#else
  {
    double factor, r;

    factor = 2.0;
    if (exp < 0)
      {
        factor = 0.5;
        exp = -exp;
      }
    r = d;
    if (exp != 0)
      {
        if ((exp & 1) != 0)
          r *= factor;
        exp >>= 1;
        while (exp != 0)
          {
            factor *= factor;
            if ((exp & 1) != 0)
              r *= factor;
            exp >>= 1;
          }
      }
    return r;
  }
#endif
}


/* End of part included from gmp */

void
mpfr_set_d(mpfr_t r, double d, unsigned char rnd_mode)
{
  int negative, sizer; unsigned int cnt;

  if (d == 0)
    {      
      EXP(r) = 0;
      return;
    }
  else if (isnan(d)) { SET_NAN(r); return; }

  negative = d < 0;
  d = ABS (d);

  sizer = MPFR_LIMBS_PER_DOUBLE; if (ABSSIZE(r)<sizer) sizer=ABSSIZE(r);
  /* warning: __mpfr_extract_double requires at least two limbs */
  EXP(r) = __mpfr_extract_double (MANT(r), d, (sizer>=2) );
  
  count_leading_zeros(cnt, MANT(r)[sizer-1]);
  if (cnt) mpn_lshift(MANT(r), MANT(r), sizer, cnt); 
  
  EXP(r) -= cnt; 
  SIZE(r) = sizer; if (negative) CHANGE_SIGN(r);

  mpfr_round(r, rnd_mode, PREC(r)); 
  return; 
}

double
mpfr_get_d2(mpfr_srcptr src, long e)
{
  double res;
  mp_size_t size, i, n_limbs_to_use;
  mp_ptr qp;
  int negative;

  if (FLAG_NAN(src)) { 
#ifdef DEBUG
    printf("recognized NaN\n");
#endif
    return NaN; }
  if (NOTZERO(src)==0) return 0.0;
  size = 1+(PREC(src)-1)/BITS_PER_MP_LIMB;
  qp = MANT(src);
  negative = (SIGN(src)==-1);

  /* Warning: don't compute the abs(res) and set the sign afterwards,
     otherwise the current machine rounding mode will not be taken
     correctly into account. */
  /* res = (negative) ? -(double)qp[size - 1] : qp[size - 1]; */
  res = 0.0;
  /* Warning: an arbitrary number of limbs may be required for an exact 
     rounding. The following code is correct but not optimal since one
     may be able to decide without considering all limbs. */
  /* n_limbs_to_use = MIN (MPFR_LIMBS_PER_DOUBLE, size); */
  n_limbs_to_use = size;
  /* Accumulate the limbs from less significant to most significant
     otherwise due to rounding we may accumulate several ulps,
     especially in rounding towards -/+infinity. */
  for (i = n_limbs_to_use; i>=1; i--)
    res = res / MP_BASE_AS_DOUBLE +
      ((negative) ? -(double)qp[size - i] : qp[size - i]);
  res = __mpfr_scale2 (res, e - BITS_PER_MP_LIMB); 

  return res;
}

double 
mpfr_get_d(mpfr_srcptr src)
{
  mpfr_get_d2(src, EXP(src));
}

