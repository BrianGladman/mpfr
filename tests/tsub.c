/* Test file for mpfr_sub.

Copyright (C) 2001 Free Software Foundation.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Library General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
License for more details.

You should have received a copy of the GNU Library General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "mpfr.h"
#include "mpfr-impl.h"
#include "mpfr-test.h"

void check_reuse _PROTO((void));

void
check_reuse ()
{
  mpfr_t x, y;

  mpfr_init2 (x, 5);
  mpfr_init2 (y, 5);
  mpfr_set_str_raw (x, "1e-12");
  mpfr_set_ui (y, 1, GMP_RNDN);
  mpfr_sub (x, y, x, GMP_RNDD);
  mpfr_set_str_raw (y, "0.11111");
  if (mpfr_cmp (x, y))
    {
      fprintf (stderr, "Error in mpfr_sub (x, y, x, GMP_RNDD) for x=2^(-12), y=1\n");
      exit (1);
    }

  mpfr_set_prec (x, 24);
  mpfr_set_prec (y, 24);
  mpfr_set_str_raw (x, "-0.100010000000000000000000E19");
  mpfr_set_str_raw (y, "0.100000000000000000000100E15");
  mpfr_add (x, x, y, GMP_RNDD);
  mpfr_set_str_raw (y, "-0.1E19");
  if (mpfr_cmp (x, y))
    {
      fprintf (stderr, "Error in mpfr_add (2)\n");
      exit (1);
    }


  mpfr_clear (x);
  mpfr_clear (y);
}

void
bug_ddefour()
{
    mpfr_t ex, ex1, ex2, ex3, tot, tot1;

    mpfr_init2(ex, 53);
    mpfr_init2(ex1, 53);
    mpfr_init2(ex2, 53);
    mpfr_init2(ex3, 53);
    mpfr_init2(tot, 150);
    mpfr_init2(tot1, 150);

    mpfr_set_ui( ex, 1, GMP_RNDN);
    mpfr_mul_2exp( ex, ex, 906, GMP_RNDN);
    mpfr_log( tot, ex, GMP_RNDN);
    mpfr_set( ex1, tot, GMP_RNDN);
    mpfr_sub( ex2, tot, ex1, GMP_RNDN);
    mpfr_sub( tot1, tot, ex1, GMP_RNDN);
    mpfr_set( ex3, tot1, GMP_RNDN);

    if (!mpfr_cmp(ex2, ex3)) 
      {
	fprintf(stderr, "Error in ddefour test.\n"); exit(-1); 
      }
}


int
main()
{
  check_reuse ();

  return 0;
}
