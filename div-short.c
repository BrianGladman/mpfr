#include <stdio.h>
#include <stdlib.h>

#include "gmp.h"
#include "gmp-impl.h"
#include "cputime.h"

#define DIVREM(qp,np,dp,nn) \
        (nn < DIV_DC_THRESHOLD) \
           ? mpn_sb_divrem_mn (qp, np, nn + nn, dp, nn) \
           : mpn_dc_divrem_n_new (qp, np, dp, nn)

static mp_limb_t
mpn_dc_divrem_n_new (mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t n)
{
  mp_size_t l, m;
  mp_limb_t qh, cc, ql;
  mp_ptr tmp;
  TMP_DECL (marker);

  l = n / 2;
  m = n - l; /* m >= l */

  TMP_MARK (marker);

  qh = DIVREM (qp + l, np + 2 * l, dp + l, m);
  /* partial remainder is in {np, 2l+m} = {np, n+l} */
  /* subtract Q1 * D0, where Q1 = {qp + l, m}, D0 = {d, l} */
  tmp = TMP_ALLOC_LIMBS (n);
  mpn_mul (tmp, qp + l, m, dp, l);
  cc = mpn_sub_n (np + l, np + l, tmp, n);
  if (qh) cc += mpn_sub_n (np + n, np + n, dp, l);
  /* have to subtract cc at np[n+l] */
  while (cc)
    {
      qh -= mpn_sub_1 (qp + l, qp + l, m, (mp_limb_t) 1);
      cc -= mpn_add_n (np + l, np + l, dp, n);
    }

  ql = DIVREM (qp, np + m, dp + m, l);
  /* partial remainder is in {np, m+l} = {np, n} */
  /* subtract Q0 * D0', where Q0 = {qp, l}, D0' = {d, m} */
  mpn_mul (tmp, dp, m, qp, l);
  cc = mpn_sub_n (np, np, tmp, n);
  TMP_FREE (marker);
  if (ql) cc += mpn_sub_n (np + l, np + l, dp, m);
  while (cc)
    {
      ql -= mpn_sub_1 (qp, qp, l, (mp_limb_t) 1);
      cc -= mpn_add_n (np, np, dp, n);
    }
  
  /* propagate ql */
  qh += mpn_add_1 (qp + l, qp + l, m, ql);
  
  return qh;
}

#define DIVREM_HIGH(qp,np,dp,nn) \
        (nn < DIV_DC_THRESHOLD) \
           ? mpn_sb_divrem_mn (qp, np, nn + nn, dp, nn) \
           : mpn_dc_divrem_n_high (qp, np, dp, nn)

static mp_limb_t
mpn_dc_divrem_n_high (mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t n)
{
  mp_size_t l, m;
  mp_limb_t qh, cc, ql;
  mp_ptr tmp;
  TMP_DECL (marker);

  l = n / 2;
  m = n - l; /* m >= l */

  TMP_MARK (marker);

  qh = DIVREM (qp + l, np + 2 * l, dp + l, m);
  /* partial remainder is in {np, 2l+m} = {np, n+l} */
  /* subtract Q1 * D0, where Q1 = {qp + l, m}, D0 = {d, l} */
  tmp = TMP_ALLOC_LIMBS (n);
  mpn_mul (tmp, qp + l, m, dp, l); /* FIXME: use a short product */
  //mpfr_mulhigh_n (tmp+m+l-2*l, qp+l+m-l, dp, l); 
  cc = mpn_sub_n (np + l, np + l, tmp, n);
  TMP_FREE (marker);
  if (qh) cc += mpn_sub_n (np + n, np + n, dp, l);
  /* have to subtract cc at np[n+l] */
  while (cc)
    {
      qh -= mpn_sub_1 (qp + l, qp + l, m, (mp_limb_t) 1);
      cc -= mpn_add_n (np + l, np + l, dp, n);
    }

  ql = DIVREM_HIGH (qp, np + m, dp + m, l);
  
  /* propagate ql */
  qh += mpn_add_1 (qp + l, qp + l, m, ql);
  
  return qh;
}

int
main (int argc, char *argv[])
{
  int n = (argc > 1) ? atoi (argv[1]) : 10000;
  int k = (argc > 2) ? atoi (argv[2]) : 1;
  mp_limb_t *n0p, *np, *n2p, *qp, *q2p, *dp;
  int st;
  int i;

  n0p = malloc (2 * n * sizeof (mp_limb_t));
  np = malloc (2 * n * sizeof (mp_limb_t));
  n2p = malloc (2 * n * sizeof (mp_limb_t));
  dp = malloc (n * sizeof (mp_limb_t));
  qp = malloc (n * sizeof (mp_limb_t));
  q2p = malloc (n * sizeof (mp_limb_t));

  mpn_random (n0p, 2 * n);
  mpn_random (dp, n);
  dp[n - 1] |= GMP_LIMB_HIGHBIT;

  printf ("DIV_DC_THRESHOLD=%u\n", DIV_DC_THRESHOLD);

  st = cputime ();
  for (i = 0; i < k; i++)
    {
      MPN_COPY (np, n0p, 2 * n);
      mpn_divrem (qp, 0, np, 2 * n, dp, n);
    }
  printf ("mpn_divrem took %dms\n", cputime () - st);

#if 0
  st = cputime ();
  for (i = 0; i < k; i++)
    {
      MPN_COPY (n2p, n0p, 2 * n);
      mpn_sb_divrem_mn (q2p, n2p, n + n, dp, n);
    }
  printf ("mpn_sb_divrem_mn took %dms\n", cputime () - st);

  if (mpn_cmp (np, n2p, n) || mpn_cmp (qp, q2p, n))
    abort ();
#endif

  st = cputime ();
  for (i = 0; i < k; i++)
    {
      MPN_COPY (n2p, n0p, 2 * n);
      mpn_dc_divrem_n_new (q2p, n2p, dp, n);
    }
  printf ("mpn_dc_divrem_n_new took %dms\n", cputime () - st);

  if (mpn_cmp (np, n2p, n) || mpn_cmp (qp, q2p, n))
    abort ();

  st = cputime ();
  for (i = 0; i < k; i++)
    {
      MPN_COPY (n2p, n0p, 2 * n);
      mpn_dc_divrem_n_high (q2p, n2p, dp, n);
    }
  printf ("mpn_dc_divrem_n_high took %dms\n", cputime () - st);

  for (i = n - 1; i >= 0 && q2p[i] == qp[i]; i--);

  if (i >= 0)
    printf ("limbs %d differ: %lu %lu\n", i, qp[i], q2p[i]);

  return 0;
}
