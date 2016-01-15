/* Test file for mpfr_log_ui.

Copyright 2016 Free Software Foundation, Inc.
Contributed by the AriC and Caramel projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include "mpfr-test.h"

#define TEST_FUNCTION mpfr_log_ui

int
main (int argc, char *argv[])
{
  unsigned int prec, yprec;
  int rnd;
  mpfr_t x, y, z, t, v;
  unsigned long m, n;
  int inex;
  mpfr_exp_t emin, emax;
  int i;

  tests_start_mpfr ();

  emin = mpfr_get_emin ();
  emax = mpfr_get_emax ();

  mpfr_inits2 (53, x, y, z, t, (mpfr_ptr) 0);
  mpfr_init2 (v, sizeof (unsigned long) * CHAR_BIT);

  if (argc >= 3) /* tlog_ui n prec [rnd] */
    {
      mpfr_set_prec (x, atoi (argv[2]));
      TEST_FUNCTION (x, atoi (argv[1]),
                    argc > 3 ? (mpfr_rnd_t) atoi (argv[3]) : MPFR_RNDN);
      mpfr_out_str (stdout, 10, 0, x, MPFR_RNDN);
      printf ("\n");
      goto clear_and_exit;
    }

  mpfr_set_prec (x, 33);
  mpfr_set_prec (y, 33);
  mpfr_log_ui (x, 3, MPFR_RNDZ);
  mpfr_set_str_binary (y, "1.0001100100111110101001111010101");
  if (mpfr_cmp (x, y))
    {
      printf ("Error for log(3), prec=33, MPFR_RNDZ\n");
      printf ("expected "); mpfr_dump (y);
      printf ("got      "); mpfr_dump (x);
      exit (1);
    }

  mpfr_set_prec (x, 60);
  mpfr_set_prec (y, 60);
  mpfr_log_ui (x, 19, MPFR_RNDU);
  mpfr_set_str_binary (y, "10.1111000111000110110000001100000010010110011001011000111010");
  if (mpfr_cmp (x, y))
    {
      printf ("Error for log(19), prec=60, MPFR_RNDU\n");
      printf ("expected "); mpfr_dump (y);
      printf ("got      "); mpfr_dump (x);
      exit (1);
    }

  mpfr_clear_flags ();
  inex = mpfr_log_ui (x, 0, MPFR_RNDN);
  MPFR_ASSERTN (inex == 0);
  MPFR_ASSERTN (mpfr_inf_p (x));
  MPFR_ASSERTN (mpfr_sgn (x) < 0);
  MPFR_ASSERTN (__gmpfr_flags == MPFR_FLAGS_DIVBY0);

  mpfr_clear_flags ();
  inex = mpfr_log_ui (x, 1, MPFR_RNDN);
  MPFR_ASSERTN (inex == 0);
  MPFR_ASSERTN (mpfr_zero_p (x));
  MPFR_ASSERTN (mpfr_signbit (x) == 0);
  MPFR_ASSERTN (__gmpfr_flags == 0);

  for (prec = MPFR_PREC_MIN; prec <= 100; prec++)
    {
      mpfr_set_prec (x, prec);
      mpfr_set_prec (z, prec);
      mpfr_set_prec (t, prec);
      yprec = prec + 20;
      mpfr_set_prec (y, yprec);

      for (m = 2; m < 70; m++)
        RND_LOOP (rnd)
          {
            /* Start with n = 2 to 49 (mpfr_can_round would fail for n < 2),
               then ULONG_MAX down to ULONG_MAX-19. */
            n = m < 50 ? m : ULONG_MAX - (m - 50);
            inex = mpfr_set_ui (v, n, MPFR_RNDN);
            MPFR_ASSERTN (inex == 0);
            mpfr_log (y, v, MPFR_RNDN);
            if (mpfr_can_round (y, yprec, MPFR_RNDN, MPFR_RNDZ, prec
                                + (rnd == MPFR_RNDN)))
              {
                mpfr_set (t, y, (mpfr_rnd_t) rnd);
                for (i = 0; i <= 1; i++)
                  {
                    if (i)
                      {
                        mpfr_exp_t e;

                        if (MPFR_IS_SINGULAR (t))
                          break;
                        e = mpfr_get_exp (t);
                        set_emin (e);
                        set_emax (e);
                      }
                    TEST_FUNCTION (z, n, (mpfr_rnd_t) rnd);
                    if (i)
                      {
                        set_emin (emin);
                        set_emax (emax);
                      }
                    if (mpfr_cmp (t, z))
                      {
                        printf ("results differ for n = %lu, prec = %u,"
                                " %s%s\n", n, prec,
                                mpfr_print_rnd_mode ((mpfr_rnd_t) rnd),
                                i ? ", reduced exponent range" : "");
                        printf ("  got      ");
                        mpfr_dump (z);
                        printf ("  expected ");
                        mpfr_dump (t);
                        printf ("  approx   ");
                        mpfr_dump (y);
                        exit (1);
                      }
                  }
              }
            else
              {
                /* We are not doing random tests. The precision increase
                   must have be chosen so that this case never occurs. */
                printf ("mpfr_can_round failed for n = %lu, prec = %u, %s\n",
                        n, prec, mpfr_print_rnd_mode ((mpfr_rnd_t) rnd));
                exit (1);
              }
          }
    }

 clear_and_exit:
  mpfr_clears (x, y, z, t, v, (mpfr_ptr) 0);

  tests_end_mpfr ();
  return 0;
}
