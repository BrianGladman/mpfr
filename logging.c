/* MPFR Logging functions.

Copyright 2005 Free Software Foundation, Inc.

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
#include "mpfr-impl.h"

/* Logging MPFR needs GCC >= 3.0 and GLIBC >= 2.0. */

#ifdef MPFR_USE_LOGGING

/* Can't include them before (in particular, printf.h) */
#include <stdlib.h>
#include <printf.h>
#include <stdarg.h>
#include <time.h>

# if !__MPFR_GLIBC(2,0)
#  error "Logging not supported. (Glibc >= 2.0)"
# endif

/* Define LOGGING variables */

FILE *mpfr_log_file;
int   mpfr_log_type;
int   mpfr_log_level;
int   mpfr_log_base;
int   mpfr_log_current;
int   mpfr_log_worstcase_limit;
mp_prec_t mpfr_log_prec;

static int 
mpfr_printf_mpfr_print (FILE *stream, const struct printf_info *info, 
			const void * const *arg)
{
  int length;
  int org_type_logging;

  /* TODO: Use much more flag from info  */
  mpfr_srcptr w = *((mpfr_srcptr *) (arg[0]));
  mp_prec_t prec = mpfr_log_prec != 0 ? mpfr_log_prec
    : info->width == -1 ? 0 : (mp_prec_t) info->width;

  org_type_logging = mpfr_log_type;
  mpfr_log_type = 0; /* We disable the logging during this print! */
  if (info->alt)
    length = fprintf (stream, "%lu", MPFR_PREC (w));
  else
    length = mpfr_out_str (stream, mpfr_log_base, prec, w, GMP_RNDN);
  mpfr_log_type = org_type_logging;

  return length;
}

static int
mpfr_printf_mpfr_arginfo (const struct printf_info *info, size_t n,
			  int *argtypes)
{
  if (n > 0)
    argtypes[0] = PA_POINTER;
  return 1;
}

static void mpfr_log_begin (void) __attribute__((constructor));

/* We let the system close the LOG itself 
   (Otherwise functions called by destructor can't use LOG File */
static void 
mpfr_log_begin (void) 
{
  const char *filename;
  time_t tt;

  filename = getenv ("MPFR_LOG_FILE");
  filename = filename != NULL && *filename != 0 ? filename : "./mpfr.log";
  mpfr_log_file = fopen (filename, "w");
  if (mpfr_log_file == NULL) {
    fprintf (stderr, "MPFR LOG: Can't open '%s' with w.\n", filename);
    abort ();
  }

  filename = getenv ("MPFR_LOG_BASE");
  mpfr_log_base = filename == NULL || *filename == 0 ? 10 : atoi (filename);

  filename = getenv ("MPFR_LOG_LEVEL");
  mpfr_log_level = filename == NULL || *filename == 0 ? 7 : atoi (filename);
  mpfr_log_current = 0;

  filename = getenv ("MPFR_LOG_PREC");
  mpfr_log_prec = filename == NULL || *filename == 0 ? 0 : atol (filename);

  mpfr_log_type = 0;
  if (getenv ("MPFR_LOG_INPUT") != NULL)
    mpfr_log_type |= MPFR_LOG_INPUT_F;
  if (getenv ("MPFR_LOG_OUTPUT") != NULL)
    mpfr_log_type |= MPFR_LOG_OUTPUT_F;
  if (getenv ("MPFR_LOG_TIME") != NULL)
    mpfr_log_type |= MPFR_LOG_TIME_F;
  if (getenv ("MPFR_LOG_INTERNAL") != NULL)
    mpfr_log_type |= MPFR_LOG_INTERNAL_F;
  if (getenv ("MPFR_LOG_MSG") != NULL)
    mpfr_log_type |= MPFR_LOG_MSG_F;
  if (getenv ("MPFR_LOG_ZIV") != NULL)
    mpfr_log_type |= MPFR_LOG_BADCASE_F;

  time (&tt);
  fprintf (mpfr_log_file, "MPFR LOG FILE %s\n", ctime (&tt));  

  /* Register printf functions */
  register_printf_function ('R', mpfr_printf_mpfr_print, 
			    mpfr_printf_mpfr_arginfo);
}

/* Return user CPU time measured in milliseconds. Thanks to Torbjorn. */

#if defined (ANSIONLY) || defined (USG) || defined (__SVR4) \
 || defined (_UNICOS) || defined(__hpux)

int mpfr_get_cputime (void) {
  return (int) ((unsigned long long) clock () * 1000 / CLOCKS_PER_SEC);
}

#else /* Use getrusage for cputime */

#include <sys/types.h>
#include <sys/resource.h>

int mpfr_get_cputime (void) {
  struct rusage rus;
  getrusage (0, &rus);
  return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}

#endif /* cputime */

#endif /* MPFR_USE_LOGGING */
