/* mpfr_sgn -- Sign of a floating point number.

Copyright 2003 Free Software Foundation.
Contributed by the Spaces project, INRIA Lorraine.

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


#include "mpfr-impl.h"

int
mpfr_sgn (mpfr_srcptr a)
{
  if (MPFR_UNLIKELY( MPFR_IS_SINGULAR(a) ))
    {
      /* Only infinite is signed */
      if (MPFR_IS_INF(a))
	return MPFR_INT_SIGN(a);
      else
	return 0;
    }
  return MPFR_INT_SIGN(a);
}
