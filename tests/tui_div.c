/* Test file for mpfr_ui_div.

Copyright 2000, 2001, 2002, 2003 Free Software Foundation, Inc.

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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "mpfr-impl.h"
#include "mpfr-test.h"

/* checks that y/x gives the same results in double
   and with mpfr with 53 bits of precision */
static void
check (unsigned long y, double x, mp_rnd_t rnd_mode, double z1)
{
  double z2;
  mpfr_t xx, zz;

  mpfr_init2 (xx, 53);
  mpfr_init2 (zz, 53);
  mpfr_set_d (xx, x, rnd_mode);
  mpfr_ui_div (zz, y, xx, rnd_mode);
  z2 = mpfr_get_d1 (zz);
  if (z1 != z2 && !(Isnan(z1) && Isnan(z2)))
    {
      printf ("expected quotient is %1.20e, got %1.20e\n", z1, z2);
      printf ("mpfr_ui_div failed for y=%lu x=%1.20e with rnd_mode=%s\n",
              y, x, mpfr_print_rnd_mode (rnd_mode));
      exit (1);
    }
  mpfr_clear (xx);
  mpfr_clear (zz);
}

static void
check_inexact (void)
{
  mpfr_t x, y, z;
  mp_prec_t px, py;
  int inexact, cmp;
  unsigned long int u;
  mp_rnd_t rnd;

  mpfr_init (x);
  mpfr_init (y);
  mpfr_init (z);

  for (px=2; px<300; px++)
    {
      mpfr_set_prec (x, px);
      do
        {
          mpfr_random (x);
        }
      while (mpfr_cmp_ui (x, 0) == 0);
      u = randlimb ();
      for (py=2; py<300; py++)
        {
          mpfr_set_prec (y, py);
          mpfr_set_prec (z, py + px);
          for (rnd=0; rnd<4; rnd++)
            {
              inexact = mpfr_ui_div (y, u, x, rnd);
              if (mpfr_mul (z, y, x, rnd))
                {
                  printf ("z <- y * x should be exact\n");
                  exit (1);
                }
              cmp = mpfr_cmp_ui (z, u);
              if (((inexact == 0) && (cmp != 0)) ||
                  ((inexact > 0) && (cmp <= 0)) ||
                  ((inexact < 0) && (cmp >= 0)))
                {
                  printf ("Wrong inexact flag for u=%lu, rnd=%s\n",
                          u, mpfr_print_rnd_mode (rnd));
                  printf ("expected %d, got %d\n", cmp, inexact);
                  printf ("x="); mpfr_print_binary (x); puts ("");
                  printf ("y="); mpfr_print_binary (y); puts ("");
                  printf ("y*x="); mpfr_print_binary (z); puts ("");
                  exit (1);
                }
            }
        }
    }

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
}

#if 0
static void
check_nan (void)
{
  mpfr_t  d, q;

  mpfr_init2 (d, 100L);
  mpfr_init2 (q, 100L);

  /* 1/+inf == 0 */
  MPFR_CLEAR_FLAGS (d);
  MPFR_SET_INF (d);
  MPFR_SET_POS (d);
  ASSERT_ALWAYS (mpfr_ui_div (q, 1L, d, GMP_RNDZ) == 0); /* exact */
  ASSERT_ALWAYS (mpfr_number_p (q));
  ASSERT_ALWAYS (mpfr_sgn (q) == 0);

  /* 1/-inf == -0 */
  MPFR_CLEAR_FLAGS (d);
  MPFR_SET_INF (d);
  MPFR_SET_NEG (d);
  ASSERT_ALWAYS (mpfr_ui_div (q, 1L, d, GMP_RNDZ) == 0); /* exact */
  ASSERT_ALWAYS (mpfr_number_p (q));
  ASSERT_ALWAYS (mpfr_sgn (q) == 0);

  /* 1/nan == nan */
  MPFR_SET_NAN (d);
  ASSERT_ALWAYS (mpfr_ui_div (q, 1L, d, GMP_RNDZ) == 0); /* exact */
  ASSERT_ALWAYS (mpfr_nan_p (q));

  /* 0/0 == nan */
  mpfr_set_ui (d, 0L, GMP_RNDN);
  ASSERT_ALWAYS (mpfr_ui_div (q, 0L, d, GMP_RNDZ) == 0); /* exact */
  ASSERT_ALWAYS (mpfr_nan_p (q));

  mpfr_clear (d);
  mpfr_clear (q);
}
#endif

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  check_inexact ();
  check(948002822, 1.22191250737771397120e+20, GMP_RNDN,
	7.758352715731357946e-12);
  check(1976245324, 1.25296395864546893357e+232, GMP_RNDZ,
	1.5772563211925444801e-223);
  check(740454110, 2.11496253355831863313e+183, GMP_RNDZ,
	3.5010270784996976041e-175);
  check(1690540942, 1.28278599852446657468e-276, GMP_RNDU,
	1.3178666932321966062e285);
  check(1476599377, -2.14191393656148625995e+305, GMP_RNDD,
	-6.8938315017943889615e-297);

  tests_end_mpfr ();
  return 0;
}
