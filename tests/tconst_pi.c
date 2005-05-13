/* Test file for mpfr_const_pi.

Copyright 1999, 2001, 2002, 2003, 2004, 2005 Free Software Foundation, Inc.

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

/* tconst_pi [prec] [rnd] [0 = no print] */

static void
check_large (void)
{
  mpfr_t x, y;

  mpfr_init2 (x, 20000);
  mpfr_init2 (y, 21000);

  mpfr_const_pi (x, GMP_RNDN); /* First one ! */
  mpfr_const_pi (y, GMP_RNDN); /* Then the other - cache - */
  mpfr_prec_round (y, 20000, GMP_RNDN);
  if (mpfr_cmp (x,y))
    {
      printf ("const_pi: error for large prec\n");
      exit (1);
    }

  /* a worst-case to exercise recomputation */
  if (MPFR_PREC_MAX > 33440) {
    mpfr_set_prec (x, 33440);
    mpfr_const_pi (x, GMP_RNDZ);
  }

  mpfr_clears (x, y, NULL);
}

int
main (int argc, char *argv[])
{
  mpfr_t x;
  int p;
  mp_rnd_t rnd;

  tests_start_mpfr ();

  p = 53;
  if (argc > 1)
    {
      int  a = atoi (argv[1]);
      if (a != 0)
        p = a;
    }

  rnd = (argc > 2) ? (mp_rnd_t) atoi(argv[2]) : GMP_RNDZ;

  mpfr_init2 (x, p);
  mpfr_const_pi (x, rnd);
  if (argc >= 2)
    {
      if (argc < 4 || atoi (argv[3]) != 0)
        {
          printf ("Pi=");
          mpfr_out_str (stdout, 10, 0, x, rnd);
          puts ("");
        }
    }
  else if (mpfr_cmp_str1 (x, "3.141592653589793116") )
    {
      printf ("mpfr_const_pi failed for prec=53\n");
      mpfr_out_str (stdout, 10, 0, x, GMP_RNDN); putchar('\n');
      exit (1);
    }

  mpfr_set_prec (x, 32);
  mpfr_const_pi (x, GMP_RNDN);
  if (mpfr_cmp_str1 (x, "3.141592653468251") )
    {
      printf ("mpfr_const_pi failed for prec=32\n");
      exit (1);
    }

  mpfr_clear (x);

  check_large();

  tests_end_mpfr ();
  return 0;
}
