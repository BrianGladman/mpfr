/* Test file for mpfr_const_log2.

Copyright 1999, 2001, 2002, 2003, 2004 Free Software Foundation, Inc.

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

#include "mpfr-test.h"

/* tlog2 [prec] [rnd] [0 = no print] */

static void
check (mp_prec_t p0, mp_prec_t p1)
{
  mpfr_t x, y, z;
  mp_rnd_t rnd;

  mpfr_init (x);
  mpfr_init (y);
  mpfr_init2 (z, p1 + 10);
  mpfr_const_log2 (z, GMP_RNDN);
  __gmpfr_const_log2_prec = 1;

  for (; p0<=p1; p0++)
    {
      mpfr_set_prec (x, p0);
      mpfr_set_prec (y, p0);
        {
          rnd = randlimb () % 4;
          mpfr_const_log2 (x, rnd);
          mpfr_set (y, z, rnd);
          if (mpfr_cmp (x, y) && mpfr_can_round (z, mpfr_get_prec(z), GMP_RNDN,
                                                 rnd, p0))
            {
              printf ("mpfr_const_log2 fails for prec=%u, rnd=%s\n",
                      (unsigned int) p0, mpfr_print_rnd_mode (rnd));
              printf ("expected ");
              mpfr_out_str (stdout, 2, 0, y, GMP_RNDN);
              printf ("\ngot      ");
              mpfr_out_str (stdout, 2, 0, x, GMP_RNDN);
              printf ("\n");
              exit (1);
            }
        }
    }

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
}

int
main (int argc, char *argv[])
{
  mpfr_t x;
  int p;
  mp_rnd_t rnd;

  tests_start_mpfr ();

  p = (argc>1) ? atoi(argv[1]) : 53;
  rnd = (argc>2) ? atoi(argv[2]) : GMP_RNDZ;

  mpfr_init (x);

  check (2, 1000);

  /* check precision of 2 bits */
  mpfr_set_prec (x, 2);
  mpfr_const_log2 (x, GMP_RNDN);
  if (mpfr_cmp_ui_2exp(x, 3, -2)) /* 3*2^-2 */
    {
      printf ("mpfr_const_log2 failed for prec=2, rnd=GMP_RNDN\n"
              "expected 0.75, got ");
      mpfr_out_str(stdout, 10, 0, x, GMP_RNDN);
      putchar('\n');
      exit (1);
    }

  if (argc>=2)
    {
      mpfr_set_prec (x, p);
      mpfr_const_log2 (x, rnd);
      printf ("log(2)=");
      mpfr_out_str (stdout, 10, 0, x, rnd);
      puts ("");
    }

  mpfr_set_prec (x, 53);
  mpfr_const_log2 (x, rnd);
  if (mpfr_cmp_str1 (x, "6.9314718055994530941e-1") )
    {
      printf ("mpfr_const_log2 failed for prec=53\n");
      exit (1);
    }

  mpfr_clear(x);

  tests_end_mpfr ();
  return 0;
}
