/* Test file for mpfr_tan.

Copyright (C) 2001 Free Software Foundation, Inc.

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

#include "gmp.h"
#include "mpfr.h"

void check53 _PROTO ((double, double, mp_rnd_t));

void check53 (double x, double tan_x, mp_rnd_t rnd_mode)
{
  mpfr_t xx, s;

  mpfr_init2 (xx, 53);
  mpfr_init2 (s, 53);
  mpfr_set_d (xx, x, rnd_mode); /* should be exact */
  mpfr_tan (s, xx, rnd_mode);
  if (mpfr_get_d (s) != tan_x && (!isnan(tan_x) || !isnan(mpfr_get_d(s)))) {
    fprintf (stderr, "mpfr_tan failed for x=%1.20e, rnd=%s\n", x,
	     mpfr_print_rnd_mode (rnd_mode));
    fprintf (stderr, "mpfr_tan gives tan(x)=%1.20e, expected %1.20e\n",
	     mpfr_get_d (s), tan_x);
    exit(1);
  }
  mpfr_clear (xx);
  mpfr_clear (s);
}

int
main(int argc, char *argv[])
{
  unsigned int prec, err, yprec, n, p0 = 1, p1 = 100, N = 100;
  mp_rnd_t rnd;
  mpfr_t x, y, z, t;

  mpfr_init (x);
  mpfr_init (y);
  mpfr_init (z);
  mpfr_init (t);

  check53(0.0/0.0, 0.0/0.0, GMP_RNDN); 
  check53(1.0/0.0, 0.0/0.0, GMP_RNDN); 
  check53(-1.0/0.0, 0.0/0.0, GMP_RNDN); 

  mpfr_set_prec (x, 1);
  mpfr_set_d (x, 0.5, GMP_RNDN);
  mpfr_tan (x, x, GMP_RNDD);
  if (mpfr_get_d(x) != 0.50)
    {
      fprintf (stderr, "mpfr_tan(0.5, GMP_RNDD) failed\n");
      exit (1);
    }

  /* generic test */
  for (prec = p0; prec <= p1; prec++)
    {
      mpfr_set_prec (x, prec);
      mpfr_set_prec (z, prec);
      mpfr_set_prec (t, prec);
      yprec = prec + 10;

      for (n=0; n<N; n++)
	{
	  mpfr_random (x);
	  rnd = random () % 4;
	  mpfr_set_prec (y, yprec);
	  mpfr_tan (y, x, rnd);
	  err = (rnd == GMP_RNDN) ? yprec + 1 : yprec;
	  if (mpfr_can_round (y, err, rnd, rnd, prec))
	    {
	      mpfr_set (t, y, rnd);
	      mpfr_tan (z, x, rnd);
	      if (mpfr_cmp (t, z))
		{
		  printf ("results differ for x=");
		  mpfr_out_str (stdout, 2, prec, x, GMP_RNDN);
		  printf (" prec=%u rnd_mode=%s\n", prec,
			  mpfr_print_rnd_mode (rnd));
		  printf ("   got ");
		  mpfr_out_str (stdout, 2, prec, z, GMP_RNDN);
		  putchar ('\n');
		  printf ("   expected ");
		  mpfr_out_str (stdout, 2, prec, t, GMP_RNDN);
		  putchar ('\n');
		  printf ("   approximation was ");
		  mpfr_print_raw (y);
		  putchar ('\n');
		  exit (1);
		}
	    }
	}
    }

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
  mpfr_clear (t);

  return 0;
}
