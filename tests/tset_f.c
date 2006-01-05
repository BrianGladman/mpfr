/* Test file for mpfr_set_f.

Copyright 1999, 2001, 2002, 2003, 2004, 2005, 2006 Free Software Foundation, Inc.

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
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h> /* for ULONG_MAX */

#include "mpfr-test.h"

int
main (void)
{
  mpfr_t x, u;
  mpf_t y, z;
  unsigned long k, pr;
  int r, inexact;

  MPFR_TEST_USE_RANDS ();
  tests_start_mpfr ();

  mpf_init (y);
  mpf_init (z);

  mpf_set_d (y, 0.0);

  /* check prototype of mpfr_init_set_f */
  mpfr_init_set_f (x, y, GMP_RNDN);
  mpfr_set_prec (x, 100);
  mpfr_set_f (x, y, GMP_RNDN);

  mpf_random2 (y, 10, 0);
  mpfr_set_f (x, y, (mp_rnd_t) RND_RAND());

  /* bug found by Jean-Pierre Merlet */
  mpfr_set_prec (x, 256);
  mpf_set_prec (y, 256);
  mpfr_init2 (u, 256);
  mpfr_set_str (u,
     "7.f10872b020c49ba5e353f7ced916872b020c49ba5e353f7ced916872b020c498@2",
     16, GMP_RNDN);
  mpf_set_str (y, "2033033E-3", 10); /* avoid 2033.033 which is
                                        locale-sensitive */
  mpfr_set_f (x, y, GMP_RNDN);
  if (mpfr_cmp (x, u))
    {
      printf ("mpfr_set_f failed for y=2033033E-3\n");
      exit (1);
    }
  mpf_set_str (y, "-2033033E-3", 10); /* avoid -2033.033 which is
                                         locale-sensitive */
  mpfr_set_f (x, y, GMP_RNDN);
  mpfr_neg (u, u, GMP_RNDN);
  if (mpfr_cmp (x, u))
    {
      printf ("mpfr_set_f failed for y=-2033033E-3\n");
      exit (1);
    }

  mpf_set_prec (y, 300);
  mpf_set_str (y, "1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111", -2);
  mpf_mul_2exp (y, y, 600);
  mpfr_set_prec (x, 300);
  mpfr_set_f (x, y, GMP_RNDN);
  if (mpfr_check (x) == 0)
    {
      printf ("Error in mpfr_set_f: corrupted result\n");
      mpfr_dump (x);
      exit (1);
    }
  MPFR_ASSERTN(mpfr_cmp_ui_2exp (x, 1, 901) == 0);

  mpfr_clear (u);

  for (k = 1; k <= 100000; k++)
    {
      pr = 2 + (randlimb () & 255);
      mpf_set_prec (z, pr);
      mpf_random2 (z, z->_mp_prec, 0);
      mpfr_set_prec (x, pr);
      mpfr_set_f (x, z, (mp_rnd_t) 0);
    }

  /* Check for +0 */
  mpfr_set_prec (x, 53);
  mpf_set_prec (y, 53);
  mpf_set_ui (y, 0);
  for (r = 0 ; r < GMP_RND_MAX ; r++)
    {
      int i;
      for (i = -1; i <= 1; i++)
        {
          if (i)
            mpfr_set_si (x, i, GMP_RNDN);
          inexact = mpfr_set_f (x, y, (mp_rnd_t) r);
          if (!MPFR_IS_ZERO(x) || !MPFR_IS_POS(x) || inexact)
            {
              printf ("mpfr_set_f(x,0) failed for %s, i = %d\n",
                      mpfr_print_rnd_mode ((mp_rnd_t) r), i);
              exit (1);
            }
        }
    }

  /* coverage test */
  mpf_set_prec (y, 2);
  mpfr_set_prec (x, 3 * mp_bits_per_limb);
  mpf_set_ui (y, 1);
  for (r = 0; r < mp_bits_per_limb; r++)
    {
      mpfr_random (x); /* to fill low limbs with random data */
      inexact = mpfr_set_f (x, y, GMP_RNDN);
      MPFR_ASSERTN(inexact == 0 && mpfr_cmp_ui_2exp (x, 1, r) == 0);
      mpf_mul_2exp (y, y, 1);
    }

  mpf_set_ui (y, 1);
  mpf_mul_2exp (y, y, ULONG_MAX);
  mpfr_set_f (x, y, GMP_RNDN);
  if (mpfr_inf_p (x) == 0 || mpfr_cmp_ui (x, 0) < 0)
    {
      printf ("Error: mpfr_set_f (x, y, GMP_RNDN) for y=2^ULONG_MAX\n");
      exit (1);
    }

  mpf_set_ui (y, 1);
  mpf_mul_2exp (y, y, __mpfr_emax);
  mpfr_set_f (x, y, GMP_RNDN);
  if (mpfr_inf_p (x) == 0 || mpfr_cmp_ui (x, 0) < 0)
    {
      printf ("Error: mpfr_set_f (x, y, GMP_RNDN) for y=2^__mpfr_emax\n");
      exit (1);
    }

  mpf_set_ui (y, 1);
  mpf_mul_2exp (y, y, __mpfr_emax - 1);
  mpfr_set_f (x, y, GMP_RNDN);
  if (mpfr_cmp_ui_2exp (x, 1, __mpfr_emax - 1) != 0)
    {
      printf ("Error: mpfr_set_f (x, y, GMP_RNDN) for y=2^(__mpfr_emax-1)\n");
      exit (1);
    }

  mpfr_clear (x);
  mpf_clear (y);
  mpf_clear (z);

  tests_end_mpfr ();
  return 0;
}
