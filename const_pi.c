/* mpfr_const_pi -- compute Pi

Copyright 1999, 2000, 2001, 2002, 2003, 2004, 2005 Free Software Foundation.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include "mpfr-impl.h"

/* Time: mpfr_const_pi2 is this new version.
   norroy% time ./pi2 1000000
   mpfr_const_pi took 14110ms
   mpfr_const_pi2 took 4730ms

   norroy% time ./pi2 2000000
   mpfr_const_pi took 39390ms
   mpfr_const_pi2 took 13140ms

   norroy% time ./pi2 5000000
   mpfr_const_pi took 252120ms
   mpfr_const_pi2 took 51780ms */

int
(mpfr_const_pi) (mpfr_ptr x, mp_rnd_t rnd_mode)
{
  mpfr_t a, A, B, D, S;
  mp_prec_t px, p, cancel, k, kmax;
  int inex, ok;

  px = MPFR_PREC (x);

  /* we need 9*2^kmax - 4 >= px+2*kmax+8 */
  for (kmax = 2; ((px + 2 * kmax + 12) / 9) >> kmax; kmax ++);

  p = px + 2 * kmax + 14; /* guarantees no recomputation for px <= 10000 */

  do {
    p += kmax;
    
    mpfr_init2 (a, p);
    mpfr_init2 (A, p);
    mpfr_init2 (B, p);
    mpfr_init2 (D, p);
    mpfr_init2 (S, p);
    
    mpfr_set_ui (a, 1, GMP_RNDN);          /* a = 1 */
    mpfr_set_ui (A, 1, GMP_RNDN);          /* A = a^2 = 1 */
    mpfr_set_ui_2exp (B, 1, -1, GMP_RNDN); /* B = b^2 = 1/2 */
    mpfr_set_ui_2exp (D, 1, -2, GMP_RNDN); /* D = 1/4 */
    
#define b B
#define ap a
#define Ap A
#define Bp B
    for (k = 0, cancel = 0; ; k++)
      {
	mpfr_add (S, A, B, GMP_RNDN);
	mpfr_div_2exp (S, S, 2, GMP_RNDN); /* exact */
	mpfr_sqrt (b, B, GMP_RNDN);
	mpfr_add (ap, a, b, GMP_RNDN);
	mpfr_div_2exp (ap, ap, 1, GMP_RNDN); /* exact */
	mpfr_mul (Ap, ap, ap, GMP_RNDN);
	mpfr_sub (Bp, Ap, S, GMP_RNDN);
	mpfr_mul_2exp (Bp, Bp, 1, GMP_RNDN);
	mpfr_sub (S, Ap, Bp, GMP_RNDN);
	cancel = mpfr_cmp_ui (S, 0) ? -mpfr_get_exp(S) : p;
	/* MPFR_ASSERTN (cancel >= px || cancel >= 9 * (1 << k) - 4); */
	mpfr_mul_2exp (S, S, k, GMP_RNDN);
	mpfr_sub (D, D, S, GMP_RNDN);
	/* stop when |A_k - B_k| <= 2^(k-p) i.e. cancel >= p-k */
	if (cancel + k >= p)
	  break;
      }
#undef b
#undef ap
#undef Ap
#undef Bp

      mpfr_div (A, B, D, GMP_RNDN);

      /* MPFR_ASSERTN(p >= 2 * k + 8); */
      if ((ok = mpfr_can_round (A, p - 2 * k - 8, GMP_RNDN, GMP_RNDZ,
                                px + (rnd_mode == GMP_RNDN))))
        inex = mpfr_set (x, A, rnd_mode);
      
      mpfr_clear (a);
      mpfr_clear (A);
      mpfr_clear (B);
      mpfr_clear (D);
      mpfr_clear (S);
  } while (ok == 0);

  return inex;
}

/*
 int
main (int argc, char *argv[])
{
  int n = atoi (argv[1]);
  int k = (argc > 2) ? atoi (argv[2]) : 1;
  int i;
  mpfr_t x, y;
  int st;
  FILE *fp;

  mpfr_init2 (x, n);
  mpfr_init2 (y, n);

#if 1
  st = cputime ();
  for (i = 0; i < k; i++)
    (mpfr_const_pi) (x, GMP_RNDN);
  printf ("mpfr_const_pi took %dms\n", cputime () - st);
#endif

  st = cputime ();
  for (i = 0; i < k; i++)
    mpfr_const_pi2 (y, GMP_RNDN);
  printf ("mpfr_const_pi2 took %dms\n", cputime () - st);

  fp = fopen ("pi", "w");
  st = cputime ();
  fprintf (fp, "pi2:=");
  mpfr_out_str (fp, 10, 0, y, GMP_RNDZ);
  fprintf (fp, ":\n");
  printf ("output took %dms\n", cputime () - st);
  fclose (fp);
  
#if 1
  if (mpfr_cmp (x, y))
    {
      fprintf (stderr, "results differ\n");
      exit (1);
    }
#endif

  mpfr_clear (x);
  mpfr_clear (y);
  
  return 0;
}
*/
