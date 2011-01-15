/* Test file for mpfr_urandom_gaussian

Copyright 1999, 2000, 2001, 2002, 2003, 2004, 2006, 2007, 2008, 2009, 2010, 2011 Free Software Foundation, Inc.
Contributed by the Arenaire and Caramel projects, INRIA.

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

#include <stdio.h>
#include <stdlib.h>

#include "mpfr-test.h"

static void
test_urandom_gaussian (long nbtests, mpfr_prec_t prec, mpfr_rnd_t rnd,
		       int verbose)
{
  mpfr_t *t;
  mpfr_t av, va, tmp;
  int i, inex, itmp;

  nbtests = (nbtests & 1) ? (nbtests + 1) : nbtests;
  t = malloc (nbtests * sizeof (mpfr_t));
  if (t == NULL)
    {
      fprintf (stderr, "turandom_gaussian: can't allocate memory in test_urandom_gaussian\n");
      exit (1);
    }
  
  for (i = 0; i < nbtests; ++i)
    mpfr_init2 (t[i], prec);
  
  inex = 1;
  for (i = 0; i < nbtests; i += 2)
    {
      itmp = mpfr_urandom_gaussian (t[i], t[i + 1], RANDS, MPFR_RNDN);
      inex = ((itmp & 3) != 0) && ((itmp & 12) != 0) && inex;
    }

  if (inex == 0)
    {
      /* one call in the loop pretended to return an exact number! */
      printf ("Error: mpfr_urandom_gaussian() returns a zero ternary value.\n");
      exit (1);
    }

  if (verbose)
    {
      mpfr_init2 (av, prec);
      mpfr_init2 (va, prec);
      mpfr_init2 (tmp, prec);

      mpfr_set_ui (av, 0, MPFR_RNDN);
      mpfr_set_ui (va, 0, MPFR_RNDN);
      for (i = 0; i < nbtests; ++i)
	{
	  mpfr_add (av, av, t[i], MPFR_RNDN);
	  mpfr_sqr (tmp, t[i], MPFR_RNDN);
	  mpfr_add (va, va, tmp, MPFR_RNDN);
	}
      mpfr_div_ui (av, av, nbtests, MPFR_RNDN);
      mpfr_div_ui (va, va, nbtests, MPFR_RNDN);
      mpfr_sqr (tmp, av, MPFR_RNDN);
      mpfr_sub (va, va, av, MPFR_RNDN);
      
      mpfr_printf ("Average = %.5Rf\nVariance = %.5Rf\n", av, va);
      mpfr_clear (av);
      mpfr_clear (va);
      mpfr_clear (tmp);
    }
  
  for (i = 0; i < nbtests; ++i)
    mpfr_clear (t[i]);
  free (t);
  return;
}


int
main (int argc, char *argv[])
{
  long nbtests;
  int verbose;
  tests_start_mpfr ();

  verbose = 0;
  nbtests = 10;
  if (argc > 1)
    {
      long a = atol (argv[1]);
      verbose = 1;
      if (a != 0)
        nbtests = a;
    }
  
  test_urandom_gaussian (nbtests, 420, MPFR_RNDN, verbose);

  tests_end_mpfr ();
  return 0;
}
