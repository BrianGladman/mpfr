/* mpfr_inits2 -- initialize several floating-point numbers with given
   precision

Copyright 2003 Free Software Foundation, Inc.

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

/* Needed for build with GMP */
#ifndef HAVE_STDARG
# include "config.h"
#endif

#if HAVE_STDARG
# include <stdarg.h>
#else
# include <varargs.h>
#endif

#include "mpfr-impl.h"

/*
 * Contrary to mpfr_init2, mp_prec_t p is the first argument
 */

/* Explicit support for K&R compiler */
void
#if HAVE_STDARG
mpfr_inits2 (mp_prec_t p, mpfr_ptr x, ...)
#else
mpfr_inits2 (va_alist)
 va_dcl
#endif
{
  va_list arg;
#if HAVE_STDARG
  va_start (arg, x);
#else
  mp_prec_t p;
  mpfr_ptr x;
  va_start(arg);
  p =  va_arg (arg, mp_prec_t);
  x =  va_arg (arg, mpfr_ptr);
#endif
  while (x != 0)
    {
      mpfr_init2 (x, p);
      x = (mpfr_ptr) va_arg (arg, mpfr_ptr);
    }
  va_end (arg);
}
