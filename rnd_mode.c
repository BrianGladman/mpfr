/* mpfr_set_machine_rnd_mode -- set the rounding mode for machine floats

Copyright 1999, 2001, 2002 Free Software Foundation, Inc.

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

#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

/* Note: "gmp-impl.h" must be included because MPFR_HAVE_FESETROUND
   may be defined from it (if MPFR is built with GMP). */

#ifdef MPFR_HAVE_FESETROUND
#include <fenv.h>
#endif

/* It is important to test each FE_* macro -- see the ISO C99 standard.
   For instance, with some ARM implementations, only FE_TONEAREST may
   be defined. */

/* sets the machine rounding mode to the value rnd_mode */
int
mpfr_set_machine_rnd_mode (mp_rnd_t rnd_mode)
{
  switch (rnd_mode)
    {
    case GMP_RNDN:
      return
#if defined (MPFR_HAVE_FESETROUND) && defined (FE_TONEAREST)
        fesetround(FE_TONEAREST)
#else
        -1
#endif
        ;
    case GMP_RNDZ:
      return
#if defined (MPFR_HAVE_FESETROUND) && defined (FE_TOWARDZERO)
        fesetround(FE_TOWARDZERO)
#else
        -1
#endif
        ;
    case GMP_RNDU:
      return
#if defined (MPFR_HAVE_FESETROUND) && defined (FE_UPWARD)
        fesetround(FE_UPWARD)
#else
        -1
#endif
        ;
    case GMP_RNDD:
      return
#if defined (MPFR_HAVE_FESETROUND) && defined (FE_DOWNWARD)
        fesetround(FE_DOWNWARD)
#else
        -1
#endif
        ;
    default:
      return -1;
    }
}
