/* Test file for mpz_set_fr / mpfr_get_z.

Copyright 2004, 2006-2016 Free Software Foundation, Inc.
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

static void
check_diff (void)
{
  int inex;
  mpfr_t x;
  mpz_t  z;
  mpfr_exp_t emin;

  mpz_init   (z);
  mpfr_init2 (x, 2);

  mpfr_set_ui (x, 2047, MPFR_RNDU);
  mpz_set_fr (z, x, MPFR_RNDN);
  if (mpz_cmp_ui (z, 2048) != 0)
    {
      printf ("get_z RU 2048 failed\n");
      exit (1);
    }

  mpfr_set_prec (x, 6);
  mpfr_set_str (x, "17.5", 10, MPFR_RNDN);
  inex = mpfr_get_z (z, x, MPFR_RNDN);
  if (inex <= 0 || mpz_cmp_ui (z, 18) != 0)
    {
      printf ("get_z RN 17.5 failed\n");
      exit (1);
    }

  /* save default emin */
  emin = mpfr_get_emin ();;

  mpfr_set_emin (17);
  mpfr_set_ui (x, 0, MPFR_RNDN);
  inex = mpfr_get_z (z, x, MPFR_RNDN);
  if (inex != 0 || mpz_cmp_ui (z, 0) != 0)
    {
      printf ("get_z 0 failed\n");
      exit (1);
    }

  /* restore default emin */
  mpfr_set_emin (emin);

  mpfr_clear (x);
  mpz_clear  (z);
}

static void
check_one (mpz_ptr z)
{
  int    inex, ex_inex, same;
  int    sh, neg;
  mpfr_t f;
  mpz_t  got, ex;

  mpfr_init2 (f, MAX (mpz_sizeinbase (z, 2), MPFR_PREC_MIN));
  mpz_init (got);
  mpz_init (ex);

  for (sh = -2*GMP_NUMB_BITS ; sh < 2*GMP_NUMB_BITS ; sh++)
    {
      inex = mpfr_set_z (f, z, MPFR_RNDN);  /* exact */
      MPFR_ASSERTN (inex == 0);

      if (sh < 0)
        {
          mpz_tdiv_q_2exp (ex, z, -sh);
          inex = mpfr_div_2exp (f, f, -sh, MPFR_RNDN);
        }
      else
        {
          mpz_mul_2exp (ex, z, sh);
          inex = mpfr_mul_2exp (f, f, sh, MPFR_RNDN);
        }
      MPFR_ASSERTN (inex == 0);

      for (neg = 0; neg <= 1; neg++)
        {
          /* Test (-1)^neg * z * 2^sh */
          int fi;
          mpfr_flags_t flags[3] = { 0, MPFR_FLAGS_ALL ^ MPFR_FLAGS_ERANGE,
                                    MPFR_FLAGS_ALL }, ex_flags, gt_flags;

          for (fi = 0; fi < numberof (flags); fi++)
            {
              ex_inex = - mpfr_cmp_z (f, ex);
              ex_flags = __gmpfr_flags = flags[fi];
              if (ex_inex != 0)
                ex_flags |= MPFR_FLAGS_INEXACT;
              inex = mpfr_get_z (got, f, MPFR_RNDZ);
              gt_flags = __gmpfr_flags;
              same = SAME_SIGN (inex, ex_inex);

              if (mpz_cmp (got, ex) != 0 || !same || gt_flags != ex_flags)
                {
                  printf ("Error in check_one for sh=%d, fi=%d\n", sh, fi);
                  printf ("     f = "); mpfr_dump (f);
                  printf ("expected "); mpz_dump (ex);
                  printf ("     got "); mpz_dump (got);
                  printf ("Expected inex ~ %d, got %d (%s)\n",
                          inex, ex_inex, same ? "OK" : "wrong");
                  printf ("Flags:\n");
                  printf ("      in"); flags_out (gt_flags);
                  printf ("expected"); flags_out (ex_flags);
                  printf ("     got"); flags_out (gt_flags);
                  exit (1);
                }
            }

          mpz_neg (ex, ex);
          mpfr_neg (f, f, MPFR_RNDN);
        }
    }

  mpfr_clear (f);
  mpz_clear (got);
  mpz_clear (ex);
}

static void
check (void)
{
  mpz_t  z;

  mpz_init (z);

  mpz_set_ui (z, 0L);
  check_one (z);

  mpz_set_si (z, 123L);
  check_one (z);

  mpz_rrandomb (z, RANDS, 2*GMP_NUMB_BITS);
  check_one (z);

  mpz_rrandomb (z, RANDS, 5*GMP_NUMB_BITS);
  check_one (z);

  mpz_clear (z);
}

static void
special (void)
{
  int inex;
  mpfr_t x;
  mpz_t z;
  int i, fi;
  int rnd;
  mpfr_exp_t e;
  mpfr_flags_t flags[3] = { 0, MPFR_FLAGS_ALL ^ MPFR_FLAGS_ERANGE,
                            MPFR_FLAGS_ALL }, ex_flags, gt_flags;

  mpfr_init2 (x, 2);
  mpz_init (z);

  RND_LOOP (rnd)
    for (i = -1; i <= 1; i++)
      for (fi = 0; fi < numberof (flags); fi++)
        {
          ex_flags = flags[fi] | MPFR_FLAGS_ERANGE;
          if (i != 0)
            mpfr_set_nan (x);
          else
            mpfr_set_inf (x, i);
          __gmpfr_flags = flags[fi];
          inex = mpfr_get_z (z, x, (mpfr_rnd_t) rnd);
          gt_flags = __gmpfr_flags;
          if (gt_flags != ex_flags || inex != 0 || mpz_cmp_ui (z, 0) != 0)
            {
              printf ("special() failed on mpfr_get_z"
                      " for %s, i = %d, fi = %d\n",
                      mpfr_print_rnd_mode ((mpfr_rnd_t) rnd), i, fi);
              printf ("Expected z = 0, inex = 0,");
              flags_out (ex_flags);
              printf ("Got      z = ");
              mpz_out_str (stdout, 10, z);
              printf (", inex = %d,", inex);
              flags_out (gt_flags);
              exit (1);
            }
          __gmpfr_flags = flags[fi];
          e = mpfr_get_z_2exp (z, x);
          gt_flags = __gmpfr_flags;
          if (gt_flags != ex_flags || e != __gmpfr_emin ||
              mpz_cmp_ui (z, 0) != 0)
            {
              printf ("special() failed on mpfr_get_z_2exp"
                      " for %s, i = %d, fi = %d\n",
                      mpfr_print_rnd_mode ((mpfr_rnd_t) rnd), i, fi);
              printf ("Expected z = 0, e = %" MPFR_EXP_FSPEC "d,",
                      (mpfr_eexp_t) __gmpfr_emin);
              flags_out (ex_flags);
              printf ("Got      z = ");
              mpz_out_str (stdout, 10, z);
              printf (", e = %" MPFR_EXP_FSPEC "d,", (mpfr_eexp_t) e);
              flags_out (gt_flags);
              exit (1);
            }
        }

  mpfr_clear (x);
  mpz_clear (z);
}

int
main (void)
{
  tests_start_mpfr ();

  check ();
  check_diff ();
  special ();

  tests_end_mpfr ();
  return 0;
}
