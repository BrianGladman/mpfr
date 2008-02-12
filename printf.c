/* mpfr_printf -- printf function and friends.

Copyright 2007, 2008 Free Software Foundation, Inc.
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
#include <errno.h>
#include "mpfr-impl.h"

#ifdef _MPFR_H_HAVE_FILE

/* Each printf-like function calls mpfr_vasprintf which
   - returns the number of characters in the returned string excluding the
   terminating null
   - returns -1 and sets the erange flag if the number of produced characters
   exceeds INT_MAX (in that case, also sets errno to EOVERFLOW in POSIX
   systems) */

int
mpfr_printf (const char *fmt, ...)
{
  va_list ap;
  char *strp;
  int ret;
  va_start (ap, fmt);

  if (mpfr_vasprintf (&strp, fmt, ap) < 0)
    {
      if (strp)
        mpfr_free_str (strp);
      va_end (ap);
      return -1;
    }

  ret = printf ("%s", strp);

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
      if (strp)
        mpfr_free_str (strp);
      return -1;
    }

  ret = printf ("%s", strp);

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
      if (strp)
        mpfr_free_str (strp);
      va_end (ap);
      return -1;
    }

  ret = fprintf (fp, "%s", strp);

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
      if (strp)
        mpfr_free_str (strp);
      return -1;
    }

  ret = fprintf (fp, "%s", strp);

  mpfr_free_str (strp);
  return ret;
}
#endif /* _MPFR_H_HAVE_FILE */

int
mpfr_sprintf (char *buf, const char *fmt, ...)
{
  char *strp;
  va_list ap;
  int ret;
  va_start (ap, fmt);

  if ((ret = mpfr_vasprintf (&strp, fmt, ap)) < 0)
    {
      if (strp)
        mpfr_free_str (strp);
      return -1;
    }

  ret = sprintf (buf, "%s", strp);
  mpfr_free_str (strp);

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
      if (strp)
        mpfr_free_str (strp);
      return -1;
    }

  ret = sprintf (buf, "%s", strp);

  mpfr_free_str (strp);
  return ret;
}

int
mpfr_snprintf (char *buf, size_t size, const char *fmt, ...)
{
  char *strp;
  va_list ap;
  int ret;
  int min_size;

  /* C99 allows SIZE to be null */
  if (size == 0)
    return 0;

  MPFR_ASSERTD (buf != NULL);

  va_start (ap, fmt);
  if ((ret = mpfr_vasprintf (&strp, fmt, ap)) < 0)
    {
      if (strp)
        mpfr_free_str (strp);
      return -1;
    }
  va_end (ap);

  min_size = ret < size ? ret : size - 1;
  strncpy (buf, strp, min_size);
  buf[min_size + 1] = '\0';
  mpfr_free_str (strp);
  return ret;
}

int
mpfr_vsnprintf (char *buf, size_t size, const char *fmt, va_list ap)
{
  char *strp;
  int ret;
  int min_size;

  /* C99 allows SIZE to be null */
  if (size == 0)
    return 0;

  MPFR_ASSERTD (buf != NULL);

  if ((ret = mpfr_vasprintf (&strp, fmt, ap)) < 0)
    {
      if (strp)
        mpfr_free_str (strp);
      return -1;
    }

  min_size = ret < size ? ret : size - 1;
  strncpy (buf, strp, min_size);
  buf[min_size + 1] = '\0';
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
      if (*pp)
        mpfr_free_str (*pp);
      *pp = NULL;
      va_end (ap);
      return -1;
    }

  va_end (ap);
  return ret;
}
#endif /* HAVE_STDARG */
