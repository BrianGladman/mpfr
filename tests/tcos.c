/* Test file for mpfr_cos.

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

void check53 (double x, double cos_x, mp_rnd_t rnd_mode)
{
  mpfr_t xx, c;

  mpfr_init2 (xx, 53);
  mpfr_init2 (c, 53);
  mpfr_set_d (xx, x, rnd_mode); /* should be exact */
  mpfr_cos (c, xx, rnd_mode);
  if (mpfr_get_d (c) != cos_x && (!isnan(cos_x) || !isnan(mpfr_get_d(c))))
    {
      fprintf (stderr, "mpfr_cos failed for x=%1.20e, rnd=%s\n", x,
	       mpfr_print_rnd_mode (rnd_mode));
      fprintf (stderr, "mpfr_cos gives cos(x)=%1.20e, expected %1.20e\n",
	       mpfr_get_d (c), cos_x);
      exit (1);
    }
  mpfr_clear (xx);
  mpfr_clear (c);
}

int
main (int argc, char *argv[])
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

  /* worst case from PhD thesis of Vincent Lefe`vre: x=8980155785351021/2^54 */
  check53 (4.984987858808754279e-1, 8.783012931285841817e-1, GMP_RNDN);
  check53 (4.984987858808754279e-1, 8.783012931285840707e-1, GMP_RNDD);
  check53 (4.984987858808754279e-1, 8.783012931285840707e-1, GMP_RNDZ);
  check53 (4.984987858808754279e-1, 8.783012931285841817e-1, GMP_RNDU);
  check53 (1.00031274099908640274,  0.540039116973283217504, GMP_RNDN);
  check53 (1.00229256850978698523,  0.538371757797526551137, GMP_RNDZ);
  check53 (1.00288304857059840103,  0.537874062022526966409, GMP_RNDZ);
  check53 (1.00591265847407274059,  0.53531755997839769456,  GMP_RNDN);

  check53 (1.00591265847407274059, 0.53531755997839769456,  GMP_RNDN);

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
	  mpfr_cos (y, x, rnd);
	  err = (rnd == GMP_RNDN) ? yprec + 1 : yprec;
	  if (mpfr_can_round (y, err, rnd, rnd, prec))
	    {
	      mpfr_set (t, y, rnd);
	      mpfr_cos (z, x, rnd);
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
