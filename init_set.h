/* Macros for combined initialization and assignement functions.

Copyright (C) 1999 PolKA project, Inria Lorraine and Loria

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

#define mpfr_init_set_si(x, i, p, rnd) {\
 mpfr_init2((x), (p)); \
 mpfr_set_si((x), (i), (rnd)); \
}

#define mpfr_init_set_ui(x, i, p, rnd) {\
 mpfr_init2((x), (p)); \
 mpfr_set_ui((x), (i), (rnd)); \
}

#define mpfr_init_set_d(x, d, p, rnd) {\
 mpfr_init2((x), (p)); \
 mpfr_set_d((x), (d), (rnd)); \
}

#define mpfr_init_set(x, y, p, rnd) {\
 mpfr_init2((x), (p)); \
 mpfr_set((x), (y), (rnd)); \
}

#define mpfr_init_set_f(x, y, p, rnd) {\
 mpfr_init2((x), (p)); \
 mpfr_set_f((x), (y), (rnd)); \
}

#define mpfr_init_set_str(x, y, p, rnd) {\
 mpfr_init2((x), (p)); \
 mpfr_set_str((x), (y), (rnd)); \
}

#define mpfr_init_set_str_raw(x, y, p, rnd) {\
 mpfr_init2((x), (p)); \
 mpfr_set_str_raw((x), (y), (rnd)); \
}

