/* Test file for mpfr_add1sp.

Copyright 2004 Free Software Foundation.

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

void check_special(void);
void check_random(mpfr_prec_t p);

int main(void)
{
  mpfr_prec_t p;

  tests_start_mpfr ();

  check_special ();
  for(p = 2 ; p < 200 ; p++)
    check_random (p);

  tests_end_mpfr ();
  return 0;
}

#define STD_ERROR \
  {\
    printf("ERROR: for %s and p=%lu and i=%d:\nB=",\
	   mpfr_print_rnd_mode(r), p, i);\
    mpfr_print_binary(b);\
    printf("\nC="); mpfr_print_binary(c);\
    printf("\nadd1  : "); mpfr_print_binary(a1);\
    printf("\nadd1sp: "); mpfr_print_binary(a2);\
    putchar('\n');\
    abort();\
  }

#define STD_ERROR2 \
  {\
    printf("ERROR: Wrong inexact flag for %s and p=%lu and i=%d:\nB=",\
	   mpfr_print_rnd_mode(r), p, i);\
    mpfr_print_binary(b);\
    printf("\nC="); mpfr_print_binary(c);\
    printf("\nA="); mpfr_print_binary(a1);\
    printf("\nAdd1: %d. Add1sp: %d\n", \
	   inexact1, inexact2); \
    abort();\
  }

#define SET_PREC(_p) \
  { p = _p; \
   mpfr_set_prec(a1, _p); mpfr_set_prec(a2, _p); \
   mpfr_set_prec(b,  _p); mpfr_set_prec(c,  _p); \
 }


void check_random(mp_prec_t p)
{
  mpfr_t a1,b,c,a2;
  mp_rnd_t r;
  int i, inexact1, inexact2;

  mpfr_inits2(p, a1,b,c,a2, NULL);

  for(i = 0 ; i < 500 ; i++)
    {
      mpfr_random(b);
      mpfr_random(c);
      if (MPFR_GET_EXP(b) < MPFR_GET_EXP(c))
	mpfr_swap(b, c);
      if (MPFR_IS_PURE_FP(b) && MPFR_IS_PURE_FP(c))
	for(r = 0 ; r < GMP_RND_MAX ; r++)
	  {
	    inexact1 = mpfr_add1(a1, b, c, r);
	    inexact2 = mpfr_add1sp(a2, b, c, r);
	    if (mpfr_cmp(a1, a2))
	      STD_ERROR;
	    if (inexact1 != inexact2)
	      STD_ERROR2;
	  }
    }

  mpfr_clears(a1,a2,b,c,NULL);
}

void check_special(void)
{
  mpfr_t a1,a2,b,c;
  mp_rnd_t r;
  mpfr_prec_t p;
  int i = -1, inexact1, inexact2;

  mpfr_inits(a1,a2,b,c,NULL);

  for(r = 0 ; r < GMP_RND_MAX ; r++)
    {
      SET_PREC(53);
      mpfr_set_str1 (b, "1@100");
      mpfr_set_str1 (c, "1@1");
      inexact1 = mpfr_add1(a1, b, c, r);
      inexact2 = mpfr_add1sp(a2, b, c, r);
      if (mpfr_cmp(a1, a2))
	STD_ERROR;
      if (inexact1 != inexact2)
	STD_ERROR2;
      mpfr_set_str_binary (b, "1E53");
      mpfr_set_str_binary (c, "1E0");
      inexact1 = mpfr_add1(a1, b, c, r);
      inexact2 = mpfr_add1sp(a2, b, c, r);
      if (mpfr_cmp(a1, a2))
        STD_ERROR;
      if (inexact1 != inexact2)
        STD_ERROR2;
    }

  mpfr_clears(a1,a2,b,c,NULL);  
}
