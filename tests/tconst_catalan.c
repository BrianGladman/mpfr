/* Test file for mpfr_const_catalan.

Copyright 2005 Free Software Foundation, Inc.

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

#include <stdlib.h>
#include "mpfr-test.h"

int
main (int argc, char *argv[])
{
  mpfr_t x, y, z, t;
  unsigned int err, prec, yprec, p0 = 2, p1 = 100;
  int rnd;

  tests_start_mpfr ();

  mpfr_init (x);
  mpfr_init (y);
  mpfr_init (z);
  mpfr_init (t);

  mpfr_set_prec (x, 32);
  mpfr_const_catalan (x, GMP_RNDN);
  mpfr_mul_2exp (x, x, 32, GMP_RNDN);
  if (mpfr_cmp_ui (x, 3934042271UL))
    {
      printf ("Error in const_catalan for prec=32\n");
      exit (1);
    }
  
  for (prec = p0; prec <= p1; prec++)
    {
      mpfr_set_prec (z, prec);
      mpfr_set_prec (t, prec);
      yprec = prec + 10;

      for (rnd = 0; rnd < GMP_RND_MAX; rnd++)
	{
	  mpfr_set_prec (y, yprec);
	  mpfr_const_catalan (y, (mp_rnd_t) rnd);
	  err = (rnd == GMP_RNDN) ? yprec + 1 : yprec;
	  if (mpfr_can_round (y, err, (mp_rnd_t) rnd, (mp_rnd_t) rnd, prec))
	    {
	      mpfr_set (t, y, (mp_rnd_t) rnd);
	      mpfr_const_catalan (z, (mp_rnd_t) rnd);
	      if (mpfr_cmp (t, z))
		{
		  printf ("results differ for prec=%u rnd_mode=%s\n", prec,
			  mpfr_print_rnd_mode ((mp_rnd_t) rnd));
		  printf ("   got      ");
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
