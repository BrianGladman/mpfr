#include <stdio.h>
#include <stdlib.h>

#if 0

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

#define MPFR_MUL_BASECASE_THREEHOLD 5
#define MPFR_MULHIGH_TAB_SIZE (sizeof(mulhigh_ktab) / sizeof(mulhigh_ktab[0]))
static short mulhigh_ktab[] = {MPFR_MULHIGH_TAB};

#else

#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"

#define MPFR_MUL_BASECASE_THREEHOLD 5
#define MPFR_MULHIGH_TAB_SIZE 1024
static short mulhigh_ktab[MPFR_MULHIGH_TAB_SIZE];

#endif


/* Put in  rp[n-1..2n-1] an approximation of the n+1 high limbs
   of {mp, n} * {np, n}. 
   The error is at worst of ln(n) for rp[n] and rp[n-1] is totally wrong */
static void
mpfr_mulhigh_n_basecase (mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n)
{
  mp_size_t i;

  rp += n-1;
  umul_ppmm (rp[1], rp[0], up[n-1], vp[0]);
  for (i = 1 ; i < n ; i++)
    rp[i+1] = mpn_addmul_1 (rp, up + (n - i - 1), i+1, vp[i]);
}

void
mpfr_mulhigh_n (mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
{
  mp_size_t k;

  k = __builtin_expect (n < MPFR_MULHIGH_TAB_SIZE, 1) 
    ?  mulhigh_ktab[n] : 2*n/3;
  if (k < 0)
    mpn_mul_basecase (rp, np, n, mp, n);
  else if (k == 0)
    mpfr_mulhigh_n_basecase (rp, np, mp, n);
  else
    {
      mp_size_t l = n - k;
      mp_limb_t cy;

      mpn_mul_n (rp + 2 * l, np + l, mp + l, k); /* fills rp[2l..2n-1] */
      mpfr_mulhigh_n (rp, np + k, mp, l);          /* fills rp[l-1..2l-1] */
      cy = mpn_add_n (rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
      mpfr_mulhigh_n (rp, np, mp + k, l);          /* fills rp[l-1..2l-1] */
      cy += mpn_add_n (rp + n - 1, rp + n - 1, rp + l - 1, l + 1);
      mpn_add_1 (rp + n + l, rp + n + l, k, cy); /* propagate carry */
    }
}


#if 1 /* Tune program */

#include "timming.h"

#define MAX_BASECASE_THREEHOLD 500
#define TOLERANCE 102/100

/* Tune: to improve */
mp_size_t
find_best_k (mp_size_t n)
{
  mp_limb_t a[n], b[n], c[2*n];
  mp_size_t k, kbest;
  unsigned long long t, tbest;
  
  if (n == 0)
    return -1;

  /* Amelioration:
      Si n > 32
      Recherche min, max dans [n-33,n-1]
      max = MAX(BITS_PER_MP_LIMB+max,n)
      Recherche entre min et max 
      Marche pas pour P4 (Trop chaotique!) */

  CALCUL_OVERHEAD;
  mpn_random (a, n);
  mpn_random (b, n);

  /* Check k == -1, mpn_mul_base_basecase */
  mulhigh_ktab[n] = -1;
  kbest = -1;
  tbest = MEASURE (mpfr_mulhigh_n (c, a, b, n) );

  /* Check k == 0, mpn_mulhigh_basecase */
  mulhigh_ktab[n] = 0;
  t = MEASURE (mpfr_mulhigh_n (c, a, b, n) );
  if (t * TOLERANCE < tbest)
    kbest = 0, tbest = t;

  /* Check Mulder */
  for (k = (n+1)/2 ; k < n ; k++) {
    mulhigh_ktab[n] = k;
    t = MEASURE (mpfr_mulhigh_n (c, a, b, n));
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
    mulhigh_ktab[k] = find_best_k (k);
    printf ("%d, ", mulhigh_ktab[k]);
    fflush (stdout);
  }
  putchar ('\n');
}


int 
main (int argc, const char *argv[])
{
  mp_limb_t *a, *b, *c, *r;
  mp_size_t an, bn, cn, size, size_end, size_step;
  int i;
  unsigned long long t1, t2, t3;   

  printf ("%s [SIZE_START [SIZE_END [SIZE_STEP]]]\n", argv[0]);

  size = 15;
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
    
    t1 = MEASURE (mpfr_mulhigh_n (c, a, b, size));
    printf ("mulhigh_n:  %7Lu ", t1);
    t2 = MEASURE (mpn_mul_n (r, a, b, size));
    printf ("mpn_mul_n:  %7Lu ", t2);
    t3 = MEASURE (mpfr_mulhigh_n_basecase (c, a, b, size));
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
#endif

#if 0
int mpfr_mul () 
{
  /* multiplies two mantissa in temporary allocated space */
  b1 = MPFR_LIKELY (bn >= cn)
    ? mpn_mul (tmp, MPFR_MANT (b), bn, MPFR_MANT (c), cn)
    : mpn_mul (tmp, MPFR_MANT (c), cn, MPFR_MANT (b), bn);


  /* TO REPLACE WITH -->; */
  
  if (MPFR_UNLIKELY (bn < cn)) 
    {
      mpfr_ptr  tmp = b;
      mp_size_t tn  = bn;
      b = c;
      c = tmp;
      bn = cn;
      cn = tn;
    }
  if (MPFR_UNLIKELY (bn > MPFR_MUL_BASECASE_THREEHOLD))
    {
      mp_size_t cancel;
      mp_prec_t prec_cn;

      prec_cn = cn*BITS_PER_MP_LIMB-MPFR_INT_CEIL_LOG2 (cn);
      /* prec_cn is the expected precision of mulhigh */

      MPFR_ASSERTD (bn >= cn);
      /* FIXME: Find best guard bits to add */
      if (MPFR_UNLIKELY (MPFR_PREC (a) > prec_cn - 4))
	/* MulHigh can't produce a roundable result.
	   Do the full multiply instead. */
	goto full_multiply;
      cancel = 0;
      if (MPFR_UNLIKELY (MPFR_PREC (a) < prec_cn - 4 -  BITS_PER_MP_LIMB))
	{
	  /* MulHigh will computes too much bits */
	  cancel = (prec_cn - 4 - MPFR_PREC (a)) / BITS_PER_MP_LIMB;
	  MPFR_ASSERTD (cancel >= 1);
	}
      mpfr_mulhigh_n (tmp+2*cancel, MPFR_MANT (b) + cancel, 
		      MPFR_MANT (c) + cancel, cn-cancel);
      /* FIXME: tn or k? */
      if (MPFR_LIKELY (mpfr_can_round_raw (tmp, k, sign, prec_cn,
					   GMP_RNDN, GMP_RNDZ,
					   MPFR_PREC(a)+(rnd_mode==GMP_RNDN))))
	b1 = tmp[2*bn-1];
      else
	goto full_multiply;
      }
    }
  else
    {
    full_multiply:
      b1 = mpn_mul (tmp, MPFR_MANT (b), bn, MPFR_MANT (c), cn);      
    }

#endif
