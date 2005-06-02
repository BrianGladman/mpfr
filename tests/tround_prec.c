/* Test file for mpfr_prec_round.

Copyright 1999, 2000, 2001, 2002, 2003, 2004 Free Software Foundation.

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
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>

#include "mpfr-test.h"

int
main (void)
{
   mpfr_t x;
   mp_exp_t emax;

   tests_start_mpfr ();

   mpfr_init (x);

   mpfr_set_nan (x);
   mpfr_prec_round (x, 2, GMP_RNDN);
   MPFR_ASSERTN(mpfr_nan_p (x));

   mpfr_set_inf (x, 1);
   mpfr_prec_round (x, 2, GMP_RNDN);
   MPFR_ASSERTN(mpfr_inf_p (x) && mpfr_sgn (x) > 0);

   mpfr_set_inf (x, -1);
   mpfr_prec_round (x, 2, GMP_RNDN);
   MPFR_ASSERTN(mpfr_inf_p (x) && mpfr_sgn (x) < 0);

   mpfr_set_ui (x, 0, GMP_RNDN);
   mpfr_prec_round (x, 2, GMP_RNDN);
   MPFR_ASSERTN(mpfr_cmp_ui (x, 0) == 0 && MPFR_IS_POS(x));

   mpfr_set_ui (x, 0, GMP_RNDN);
   mpfr_neg (x, x, GMP_RNDN);
   mpfr_prec_round (x, 2, GMP_RNDN);
   MPFR_ASSERTN(mpfr_cmp_ui (x, 0) == 0 && MPFR_IS_NEG(x));

   emax = mpfr_get_emax ();
   set_emax (0);
   mpfr_set_prec (x, 3);
   mpfr_set_str_binary (x, "0.111");
   mpfr_prec_round (x, 2, GMP_RNDN);
   MPFR_ASSERTN(mpfr_inf_p (x) && mpfr_sgn (x) > 0);
   set_emax (emax);

   mpfr_set_prec (x, mp_bits_per_limb + 2);
   mpfr_set_ui (x, 1, GMP_RNDN);
   mpfr_nextbelow (x);
   mpfr_prec_round (x, mp_bits_per_limb + 1, GMP_RNDN);
   MPFR_ASSERTN(mpfr_cmp_ui (x, 1) == 0);

   mpfr_set_prec (x, 3);
   mpfr_set_ui (x, 5, GMP_RNDN);
   mpfr_prec_round (x, 2, GMP_RNDN);
   if (mpfr_cmp_ui(x, 4))
     {
       printf ("Error in tround: got ");
       mpfr_out_str (stdout, 10, 0, x, GMP_RNDN);
       printf (" instead of 4\n");
       exit (1);
     }

   /* check case when reallocation is needed */
   mpfr_set_prec (x, 3);
   mpfr_set_ui (x, 5, GMP_RNDN); /* exact */
   mpfr_prec_round (x, mp_bits_per_limb + 1, GMP_RNDN);
   if (mpfr_cmp_ui(x, 5))
     {
       printf ("Error in tround: got ");
       mpfr_out_str (stdout, 10, 0, x, GMP_RNDN);
       printf (" instead of 5\n");
       exit (1);
     }

   mpfr_clear(x);
   mpfr_init2 (x, 3);
   mpfr_set_si (x, -5, GMP_RNDN); /* exact */
   mpfr_prec_round (x, mp_bits_per_limb + 1, GMP_RNDN);
   if (mpfr_cmp_si(x, -5))
     {
       printf ("Error in tround: got ");
       mpfr_out_str (stdout, 10, 0, x, GMP_RNDN);
       printf (" instead of -5\n");
       exit (1);
     }

   /* check case when new precision needs less limbs */
   mpfr_set_prec (x, mp_bits_per_limb + 1);
   mpfr_set_ui (x, 5, GMP_RNDN); /* exact */
   mpfr_prec_round (x, 3, GMP_RNDN); /* exact */
   if (mpfr_cmp_ui(x, 5))
     {
       printf ("Error in tround: got ");
       mpfr_out_str (stdout, 10, 0, x, GMP_RNDN);
       printf (" instead of 5\n");
       exit (1);
     }

   mpfr_clear(x);

   tests_end_mpfr ();
   return 0;
}
