#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"

#include "timming.h"

#define MAX_BASECASE_THREEHOLD 500
#define TOLERANCE 102/100

/*

               <----------- un ----------->
               ___________________________
              /            .             /  ^
            / .            .           /    |
          /   .            .         /      |
        /     .            .       /        vn
      /       .            .     /          |
    /         .            .   /            |
  /           .            . /              |
/__________________________/                v

<----- vn ----><-- un-vn --> <---- vn--- ->

rn ------>            ignore low part of v and u for all partial products
rn --------------->   ignore low part of u for all partial products
rn -----------------------------> ignore low part of u for some partial products

Algorithm:
  1. Use common mpn_mul_1/mpn_addmul_1 loops, but exclude results more than one
     limbs under result part.
  2. If (lowest limb) > (number of rows of first ignored column), compute
     another column, propagate carry.
*/

#if 0
void
mpn_mulhigh_n_basecase (mp_ptr rp, mp_ptr up, mp_ptr vp, mp_size_t n)
{
  mp_size_t i;

  {
    mp_limb_t tmp1, tmp2, tmp3;
    umul_ppmm (tmp1, rp[0], up[n-2], vp[0]);
    umul_ppmm (tmp3, tmp2, up[n-1], vp[0]);
    add_ssaaaa (rp[2], rp[1], tmp3, tmp2, 0, tmp1);
  }
  for (i = 3 ; i <= n ; i++)
    rp[i] = mpn_addmul_1 (rp, up + n - i, i, vp[i - 2]);
  rp[n+1] = mpn_addmul_1 (rp + 1, up, n, vp[n - 1]);
}
#endif

/* Put in  rp[n-1..2n-1] an approximation of the n+1 high limbs
   of {mp, n} * {np, n}. 
   The error is at worst of ln(n) for rp[n] and rp[n-1] is totally wrong */
void
mpn_mulhigh_n_basecase (mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n)
{
  mp_size_t i;

  rp += n-1;
  umul_ppmm (rp[1], rp[0], up[n-1], vp[0]);
  for (i = 1 ; i < n ; i++)
    rp[i+1] = mpn_addmul_1 (rp, up + (n - i - 1), i+1, vp[i]);
}

mp_size_t mpn_mulhigh_ktab[1024];

void
mpn_mulhigh_n (mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
  mp_size_t k;

  if (n >= 1024)
    abort ();

  k = mpn_mulhigh_ktab[n];
  if (k < 0)
    mpn_mul_basecase (rp, np, n, mp, n);
  else if (k == 0)
    mpn_mulhigh_n_basecase (rp, np, mp, n);
  else
    {
      mp_size_t l = n - k;
      mp_limb_t cy;

      mpn_mul_n (rp + 2 * l, np + l, mp + l, k); /* fills rp[2l..2n-1] */
      mpn_mulhigh_n (rp, np + k, mp, l);          /* fills rp[l-1..2l-1] */
      cy = mpn_add_n (rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
      mpn_mulhigh_n (rp, np, mp + k, l);          /* fills rp[l-1..2l-1] */
      cy += mpn_add_n (rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
      mpn_add_1 (rp + n + l, rp + n + l, k, cy); /* propagate carry */
    }
}

mp_size_t
find_best_k (mp_size_t n)
{
  mp_limb_t a[n], b[n], c[2*n];
  mp_size_t k, kbest;
  unsigned long long t, tbest;
  
  if (n == 0)
    return -1;

  CALCUL_OVERHEAD;
  mpn_random (a, n);
  mpn_random (b, n);

  /* Check k == -1, mpn_mul_base_basecase */
  mpn_mulhigh_ktab[n] = -1;
  kbest = -1;
  tbest = MEASURE (mpn_mulhigh_n (c, a, b, n) );

  /* Check k == 0, mpn_mulhigh_basecase */
  mpn_mulhigh_ktab[n] = 0;
  t = MEASURE (mpn_mulhigh_n (c, a, b, n) );
  if (t * TOLERANCE < tbest)
    kbest = 0, tbest = t;

  /* Check Mulder */
  for (k = (n+1)/2 ; k < n ; k++) {
    mpn_mulhigh_ktab[n] = k;
    t = MEASURE (mpn_mulhigh_n (c, a, b, n));
    if (t *TOLERANCE < tbest)
      kbest = k, tbest = t;
  }
  return kbest;
}

void
tune (mp_size_t n)
{
  /* Find ThreeHold */
  mp_size_t k;
  for (k = 0 ; k <= n ; k++) {
    mpn_mulhigh_ktab[k] = find_best_k (k);
    printf ("%d, ", mpn_mulhigh_ktab[k]);
    fflush (stdout);
  }
  putchar ('\n');
}

#ifndef SIZE  
#define SIZE 15
#endif


int 
main (int argc, const char *argv[])
{
  mp_limb_t *a, *b, *c, *r;
  mp_size_t an, bn, cn, size, size_end, size_step;
  int i;
  unsigned long long t1, t2, t3;   

  printf ("%s [SIZE_START [SIZE_END [SIZE_STEP]]]\n", argv[0]);

  size = SIZE;
  if (argc >= 2)
    size = atol (argv[1]);
  size_end = size;
  if (argc >= 3)
    size_end = atol (argv[2]);
  size_step = 1;
  if (argc >= 4)
    size_step = atol (argv[3]);

  printf ("Tune in progress...\n");
  tune (size_end);
  
  for ( ; size <= size_end ; size += size_step) {
    printf ("Size= %4u ", size);
    
    a = malloc (sizeof(mp_limb_t)*size);
    b = malloc (sizeof(mp_limb_t)*size);
    c = malloc (sizeof(mp_limb_t)*size*3);
    r = malloc (sizeof(mp_limb_t)*size*3);
    
    CALCUL_OVERHEAD;
    mpn_random (a, size);
    mpn_random (b, size);
    
    t1 = MEASURE (mpn_mulhigh_n (c, a, b, size));
    printf ("mulhigh_n:  %7Lu ", t1);
    t2 = MEASURE (mpn_mul_n (r, a, b, size));
    printf ("mpn_mul_n:  %7Lu ", t2);
    t3 = MEASURE (mpn_mulhigh_n_basecase (c, a, b, size));
    printf ("mulhigh_bc: %7Lu ", t3);

    printf ("Ratio: %f %c %c\n", (double) t2 / t1,
	    t1 < t2 ? '*' : ' ',
	    t1 < t3 ? '$' : ' ' );

    if (size < 50 && size == size_end)
      {
	printf ("High: ");
	for (i = 2*size-1 ; i>=size-1 ; i--)
	  printf ("%08X ", c[i]);
	printf("\nmuln: ");
	for (i = 2*size-1 ; i>=size-1 ; i--)
	  printf ("%08X ", r[i]);
	printf("\n");
      }
    free (a);
    free (b);
    free (c);
    free (r);
  }
  
  return 0;
}

#if 0
int mul () 
{
  /* multiplies two mantissa in temporary allocated space */
  b1 = MPFR_LIKELY (bn >= cn)
    ? mpn_mul (tmp, MPFR_MANT (b), bn, MPFR_MANT (c), cn)
    : mpn_mul (tmp, MPFR_MANT (c), cn, MPFR_MANT (b), bn);

  -->;
  
  if (MPFR_LIKELY (bn == cn)) {
    mp_prec_t prec_bn;
    if (MPFR_LIKELY (bn < MPFR_MUL_BASECASE_THREEHOLD))
      goto mul_normal;
    prec_bn = bn*BITS_PER_MP_LIMB-MPFR_INT_CEIL_LOG2 (bn);
    if (MPFR_UNLIKELY (MPFR_PREC (a) > prec_bn-4))
      goto mul_normal;
    else {
      /* mp_size_t offset;
	 offset = (prec_bn - 4 - MPFR_PREC (a)) / BITS_PER_MP_LIMB; */
      mpfr_mpn_mulhigh_n (tmp, MPFR_MANT (b), MPFR_MANT (c), bn);
      if (MPFR_LIKELY (mpfr_can_round_raw (tmp, bn+tn, sign, prec_bn,
					   GMP_RNDN, GMP_RNDZ,
					   MPFR_PREC(a)+(rnd_mode==GMP_RNDN))))
	b1 = tmp[2*bn-1];
      else
	goto mul_normal;
    }
  } else if (bn > cn) {
  mul_normal:
    b1 = mpn_mul (tmp, MPFR_MANT (b), bn, MPFR_MANT (c), cn);
  } else
    b1 = mpn_mul (tmp, MPFR_MANT (c), cn, MPFR_MANT (b), bn);
}
#endif
