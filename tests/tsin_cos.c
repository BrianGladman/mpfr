/* Test file for mpfr_sin_cos.

Copyright (C) 2000 Free Software Foundation, Inc.

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
#include "cputime.h"

void large_test (int prec, int N)
{
  int i, st;
  mpfr_t x, s, c;

  mpfr_init2 (x, prec);
  mpfr_init2 (s, prec);
  mpfr_init2 (c, prec);
  mpfr_set_d (x, 3.0, GMP_RNDN);
  mpfr_sqrt (x, x, GMP_RNDN);
  st = cputime ();
  for (i=0; i<N; i++) mpfr_sin_cos (s, c, x, GMP_RNDN);
  st = cputime () - st;
  mpfr_out_str (stdout, 10, 0, s, GMP_RNDN); putchar('\n');
  fprintf(stderr, "sin_cos took %dms\n", st);
  mpfr_clear (x);
  mpfr_clear (s);
  mpfr_clear (c);
}

void check53 (double x, double sin_x, double cos_x, mp_rnd_t rnd_mode)
{
  mpfr_t xx, s, c;

  mpfr_init2 (xx, 53);
  mpfr_init2 (s, 53);
  mpfr_init2 (c, 53);
  mpfr_set_d (xx, x, rnd_mode); /* should be exact */
  mpfr_sin_cos (s, c, xx, rnd_mode);
  if (mpfr_get_d (s) != sin_x && (!isnan(sin_x) || !isnan(mpfr_get_d(s)))) {
    fprintf (stderr, "mpfr_sin_cos failed for x=%1.20e, rnd=%s\n", x,
	     mpfr_print_rnd_mode (rnd_mode));
    fprintf (stderr, "mpfr_sin_cos gives sin(x)=%1.20e, expected %1.20e\n",
	     mpfr_get_d (s), sin_x);
    exit(1);
  }
  if (mpfr_get_d (c) != cos_x && (!isnan(cos_x) || !isnan(mpfr_get_d(c)))) {
    fprintf (stderr, "mpfr_sin_cos failed for x=%1.20e, rnd=%s\n", x,
	     mpfr_print_rnd_mode (rnd_mode));
    fprintf (stderr, "mpfr_sin_cos gives cos(x)=%1.20e, expected %1.20e\n",
	     mpfr_get_d (c), cos_x);
    exit(1);
  }
  mpfr_clear (xx);
  mpfr_clear (s);
  mpfr_clear (c);
}

/* tsin_cos prec [N] performs N tests with prec bits */
int main(int argc, char *argv[])
{
  if (argc > 1) {
    large_test (atoi (argv[1]), (argc > 2) ? atoi (argv[2]) : 1);
  }

  check53(0.0/0.0, 0.0/0.0, 0.0/0.0, GMP_RNDN); 
  check53(1.0/0.0, 0.0/0.0, 0.0/0.0, GMP_RNDN); 
  check53(-1.0/0.0, 0.0/0.0, 0.0/0.0, GMP_RNDN); 
  /* worst case from PhD thesis of Vincent Lefe`vre: x=8980155785351021/2^54 */
  check53 (4.984987858808754279e-1, 4.781075595393330379e-1, 
	   8.783012931285841817e-1, GMP_RNDN);
  check53 (4.984987858808754279e-1, 4.781075595393329824e-1,
	   8.783012931285840707e-1, GMP_RNDD);
  check53 (4.984987858808754279e-1, 4.781075595393329824e-1,
	   8.783012931285840707e-1, GMP_RNDZ);
  check53 (4.984987858808754279e-1, 4.781075595393330379e-1,
	   8.783012931285841817e-1, GMP_RNDU);
  return 0;
}
