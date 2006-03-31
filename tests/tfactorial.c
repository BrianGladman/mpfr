/* Test file for mpfr_factorial.

Copyright 2001, 2002, 2003, 2004, 2005, 2006 Free Software Foundation, Inc.

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
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>

#include "mpfr-test.h"

#define TEST_FUNCTION mpfr_fac_ui

static void
special (void)
{
  mpfr_t x, y;
  int inex;

  mpfr_init (x);
  mpfr_init (y);

  mpfr_set_prec (x, 21);
  mpfr_set_prec (y, 21);
  mpfr_fac_ui (x, 119, GMP_RNDZ);
  mpfr_set_str_binary (y, "0.101111101110100110110E654");
  if (mpfr_cmp (x, y))
    {
      printf ("Error in mpfr_fac_ui (119)\n");
      exit (1);
    }

  mpfr_set_prec (y, 206);
  inex = mpfr_fac_ui (y, 767, GMP_RNDN);
  mpfr_set_prec (x, 206);
  mpfr_set_str_binary (x, "0.110111100001000001101010010001000111000100000100111000010011100011011111001100011110101000111101101100110001001100110100001001111110000101010000100100011100010011101110000001000010001100010000101001111E6250");
  if (mpfr_cmp (x, y))
    {
      printf ("Error in mpfr_fac_ui (767)\n");
      exit (1);
    }
  if (inex <= 0)
    {
      printf ("Wrong flag for mpfr_fac_ui (767)\n");
      exit (1);
    }

  mpfr_set_prec (y, 202);
  mpfr_fac_ui (y, 69, GMP_RNDU);

  mpfr_clear (x);
  mpfr_clear (y);
}

static void
test_int (void)
{
  unsigned long n0 = 1, n1 = 80, n;
  mpz_t f;
  mpfr_t x, y;
  mp_prec_t prec_f, p;
  int r;
  int inex1, inex2;

  mpz_init (f);
  mpfr_init (x);
  mpfr_init (y);

  mpz_fac_ui (f, n0 - 1);
  for (n = n0; n <= n1; n++)
    {
      mpz_mul_ui (f, f, n); /* f = n! */
      prec_f = mpz_sizeinbase (f, 2) - mpz_scan1 (f, 0);
      for (p = MPFR_PREC_MIN; p <= prec_f; p++)
        {
          mpfr_set_prec (x, p);
          mpfr_set_prec (y, p);
          for (r = 0; r < GMP_RND_MAX; r++)
            {
              inex1 = mpfr_fac_ui (x, n, (mp_rnd_t) r);
              inex2 = mpfr_set_z (y, f, (mp_rnd_t) r);
              if (mpfr_cmp (x, y))
                {
                  printf ("Error for n=%lu prec=%lu rnd=%s\n",
                          n, (unsigned long) p, mpfr_print_rnd_mode ((mp_rnd_t) r));
                  exit (1);
                }
              if ((inex1 < 0 && inex2 >= 0) || (inex1 == 0 && inex2 != 0)
                  || (inex1 > 0 && inex2 <= 0))
                {
                  printf ("Wrong inexact flag for n=%lu prec=%lu rnd=%s\n",
                          n, (unsigned long) p, mpfr_print_rnd_mode ((mp_rnd_t) r));
                  exit (1);
                }
            }
        }
    }

  mpz_clear (f);
  mpfr_clear (x);
  mpfr_clear (y);
}

int
main (int argc, char *argv[])
{
  unsigned int prec, err, yprec, n, k, zeros;
  int rnd;
  mpfr_t x, y, z, t;
  int inexact;

  tests_start_mpfr ();

  special ();

  test_int ();

  mpfr_init (x);
  mpfr_init (y);
  mpfr_init (z);
  mpfr_init (t);

  mpfr_fac_ui (y, 0, GMP_RNDN);

  if (mpfr_cmp_ui (y, 1))
    {
      printf ("mpfr_fac_ui(0) does not give 1\n");
      exit (1);
    }

  for (prec = 2; prec <= 100; prec++)
    {
      mpfr_set_prec (x, prec);
      mpfr_set_prec (z, prec);
      mpfr_set_prec (t, prec);
      yprec = prec + 10;
      mpfr_set_prec (y, yprec);

      for (n = 0; n < 50; n++)
        for (rnd = 0; rnd < GMP_RND_MAX; rnd++)
          {
            inexact = mpfr_fac_ui (y, n, (mp_rnd_t) rnd);
            err = (rnd == GMP_RNDN) ? yprec + 1 : yprec;
            if (mpfr_can_round (y, err, (mp_rnd_t) rnd, (mp_rnd_t) rnd, prec))
              {
                mpfr_set (t, y, (mp_rnd_t) rnd);
                inexact = mpfr_fac_ui (z, n, (mp_rnd_t) rnd);
                /* fact(n) ends with floor(n/2)+floor(n/4)+... zeros */
                for (k=n/2, zeros=0; k; k >>= 1)
                  zeros += k;
                if (MPFR_EXP(y) <= (mp_exp_t) (prec + zeros))
                  /* result should be exact */
                  {
                    if (inexact)
                      {
                        printf ("Wrong inexact flag: expected exact\n");
                        exit (1);
                      }
                  }
                else /* result is inexact */
                  {
                    if (!inexact)
                      {
                        printf ("Wrong inexact flag: expected inexact\n");
                        printf ("n=%u prec=%u\n", n, prec);
                        mpfr_print_binary(z); puts ("");
                        exit (1);
                      }
                  }
                if (mpfr_cmp (t, z))
                  {
                    printf ("results differ for x=");
                    mpfr_out_str (stdout, 2, prec, x, GMP_RNDN);
                    printf (" prec=%u rnd_mode=%s\n", prec,
                            mpfr_print_rnd_mode ((mp_rnd_t) rnd));
                    printf ("   got ");
                    mpfr_out_str (stdout, 2, prec, z, GMP_RNDN);
                    puts ("");
                    printf ("   expected ");
                    mpfr_out_str (stdout, 2, prec, t, GMP_RNDN);
                    puts ("");
                    printf ("   approximation was ");
                    mpfr_print_binary (y);
                    puts ("");
                    exit (1);
                  }
              }
          }
    }

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
  mpfr_clear (t);

  tests_end_mpfr ();
  return 0;
}
