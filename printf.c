/* mpfr_printf -- printf function and friends.

Copyright 2007 Free Software Foundation, Inc.
Contributed by the Arenaire and Cacao projects, INRIA.

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
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

/* The mpfr_printf-like functions are defined only if stdarg.h exist */
#ifdef HAVE_STDARG

#include <stdarg.h>
#include "mpfr-impl.h"

#ifdef _MPFR_H_HAVE_FILE
int
mpfr_printf (const char *fmt, ...)
{
  va_list ap;
  char *strp;
  int ret;
  va_start (ap, fmt);

  if (mpfr_vasprintf (&strp, fmt, ap) < 0)
    {
      mpfr_free_str (strp);
      va_end (ap);
      return -1;
    }

  ret = printf (strp);

  mpfr_free_str (strp);
  va_end (ap);
  return ret;
}

int
mpfr_vprintf (const char *fmt, va_list ap)
{
  char *strp;
  int ret;

  if (mpfr_vasprintf (&strp, fmt, ap) < 0)
    {
      mpfr_free_str (strp);
      return -1;
    }

  ret = printf (strp);

  mpfr_free_str (strp);
  return ret;
}


int
mpfr_fprintf (FILE *fp, const char *fmt, ...)
{
  va_list ap;
  char *strp;
  int ret;
  va_start (ap, fmt);

  if (mpfr_vasprintf (&strp, fmt, ap) < 0)
    {
      mpfr_free_str (strp);
      va_end (ap);
      return -1;
    }

  ret = fprintf (fp, strp);

  mpfr_free_str (strp);
  va_end (ap);
  return ret;
}

int
mpfr_vfprintf (FILE *fp, const char *fmt, va_list ap)
{
  char *strp;
  int ret;

  if (mpfr_vasprintf (&strp, fmt, ap) < 0)
    {
      mpfr_free_str (strp);
      return -1;
    }

  ret = fprintf (fp, strp);

  mpfr_free_str (strp);
  return ret;
}
#endif /* _MPFR_H_HAVE_FILE */

int
mpfr_sprintf (char *buf, const char *fmt, ...)
{
  va_list ap;
  int ret;
  va_start (ap, fmt);

  ret  = mpfr_vsprintf (buf, fmt, ap);

  va_end (ap);
  return ret;
}

int
mpfr_vsprintf (char *buf, const char *fmt, va_list ap)
{
  char *strp;
  int ret;

  if ((ret = mpfr_vasprintf (&strp, fmt, ap)) < 0)
    {
      mpfr_free_str (strp);
      return -1;
    }

  ret = sprintf (buf, strp);

  mpfr_free_str (strp);
  return ret;
}

int
mpfr_snprintf (char *buf, size_t size, const char *fmt, ...)
{
  va_list ap;
  int ret;
  va_start (ap, fmt);

  ret = mpfr_vsnprintf (buf, size, fmt, ap);

  va_end (ap);
  return ret;
}

int
mpfr_vsnprintf (char *buf, size_t size, const char *fmt, va_list ap)
{
  char *strp;
  int ret;

  if (mpfr_vasprintf (&strp, fmt, ap) < 0)
    {
      mpfr_free_str (strp);
      return -1;
    }

  /* FIXME: snprintf is new in C99. It must be checked with a configure
     test. But since we need our own snprintf function for systems that
     do not have it, it is probably better to write a version that will
     check that the return value is not larger than INT_MAX and use it
     unconditionally. */
  /* FIXME: warning: passing argument 3 of 'snprintf' from
     incompatible pointer type
     Remove the "&"? (The code lacks of comments...) */
  ret = snprintf (buf, size, &strp);
  /* FIXME: if we do not use our own snprintf function, here's a proposed
     solution to detect the case where the ideal return value is larger
     than INT_MAX (for POSIX systems):
       - include <errno.h>;
       - if EOVERFLOW is defined, set errno to 0 before the call;
       - after the call, if EOVERFLOW is defined and if errno is equal
         to EOVERFLOW, then set the erange flag. */

  mpfr_free_str (strp);
  return ret;
}

int
mpfr_asprintf (char **pp, const char *fmt, ...)
{
  va_list ap;
  int ret;
  va_start (ap, fmt);

  if ((ret = mpfr_vasprintf (pp, fmt, ap)) < 0)
    {
      mpfr_free_str (*pp);
      *pp = NULL;
      va_end (ap);
      return -1;
    }

  va_end (ap);
  return ret;
}
#endif /* HAVE_STDARG */
