/* mpfr_vasprintf -- main function for the printf functions family 
   plus helper macros & functions.

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

#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif

#include <stdlib.h>  /* for abs */
#include <ctype.h>   /* for toupper */
#include <stddef.h>
#include <math.h>

#include "mpfr-impl.h"
/* #include "mpfr-gmp.h" */

#if   _MPFR_PREC_FORMAT == 1
#define _MPFR_PREC_FORMAT_SPEC "hu"
#elif _MPFR_PREC_FORMAT == 2
#define _MPFR_PREC_FORMAT_SPEC "u"
#elif _MPFR_PREC_FORMAT == 3
#define _MPFR_PREC_FORMAT_SPEC "lu"
#else
#error "mpfr_prec_t size not supported"
#endif

#if (__GMP_MP_SIZE_T_INT == 1)
#define _MP_EXP_FORMAT_SPEC "i"
#elif (__GMP_MP_SIZE_T_INT == 0)
#define _MP_EXP_FORMAT_SPEC "li"
#else
#error "mp_exp_t size not supported"
#endif


/* We assume that a single conversion specifier produces at most 4095 chars
   (Rationale for International Standard —Programming Languages— C 
   Revision 5.10 April-2003, 7.19.6.1 p152). */
#define MAX_CHAR_PRODUCED_BY_SPEC 4096

/* some macro and functions for parsing format string */
#define READ_INT(format, specinfo, field, label_out)    \
  do {							\
    while (*(format))					\
      {							\
	switch (*(format))				\
	  {						\
	  case '0':					\
	  case '1':					\
	  case '2':					\
	  case '3':					\
	  case '4':					\
	  case '5':					\
	  case '6':					\
	  case '7':					\
	  case '8':					\
	  case '9':					\
	    (specinfo).field *= 10;			\
	    (specinfo).field += *(format) - '0';	\
	    ++(format);					\
	    break;					\
	  case '*':					\
	    (specinfo).field = va_arg (ap, int);	\
	    ++(format);					\
	  default:					\
	    goto label_out;				\
	  }						\
      }							\
  } while (0)

/* __arg_type contains all the types described by the 'type' field of the 
   format string */
enum __arg_type
{
  NONE,
  CHAR_ARG,
  SHORT_ARG,
  LONG_ARG,
  LONG_LONG_ARG,
#ifdef HAVE_STDINT_H
  INTMAX_ARG,
#endif
  SIZE_ARG,
  PTRDIFF_ARG,
  LONG_DOUBLE_ARG,
  MPF_ARG,
  MPQ_ARG,
  MP_LIMB_ARG,
  MP_LIMB_ARRAY_ARG,
  MPZ_ARG,
  MPFR_PREC_ARG,
  MPFR_ARG
};

/* Each conversion specification of the format string will be translated in a
   printf_spec structure by the parser.
   This structure is adapted from the gnu libc one. */
struct printf_spec
{
  unsigned int alt:1;           /* # flag */
  unsigned int space:1;         /* Space flag */
  unsigned int left:1;          /* - flag */
  unsigned int showsign:1;      /* + flag */
  unsigned int group:1;         /* ' flag (not for gmp nor mpfr types) */

  int width;                    /* Width */
  int prec;                     /* Precision */

  enum __arg_type arg_type;     /* Type of argument */
  mp_rnd_t rnd_mode;            /* Rounding mode */
  char spec;                    /* Conversion specifier */

  char pad;                     /* Padding character */
};

static void
specinfo_init (struct printf_spec *specinfo)
{
  specinfo->alt = 0;
  specinfo->space = 0;
  specinfo->left = 0;
  specinfo->showsign = 0;
  specinfo->group = 0;
  specinfo->width = 0;
  specinfo->prec = 0;
  specinfo->arg_type = NONE;
  specinfo->rnd_mode = GMP_RNDN;
  specinfo->spec = 'i';
  specinfo->pad = ' ';
}

static __gmp_const char *
parse_flags (__gmp_const char *format, struct printf_spec *specinfo)
{
  while (*format)
    {
      switch (*format)
        {
        case '0':
          specinfo->pad = '0';
          ++format;
          break;
        case '#':
          specinfo->alt = 1;
          ++format;
          break;
        case '+':
          specinfo->showsign = 1;
          ++format;
          break;
        case ' ':
          specinfo->space = 1;
          ++format;
          break;
        case '-':
          specinfo->left = 1;
          ++format;
          break;
        case '\'':
          ++format;
          break;
        default:
          return format;
        }
    }
  return format;
}

static __gmp_const char *
parse_arg_type (__gmp_const char *format, struct printf_spec *specinfo)
{
  switch (*format)
    {
    case '\0':
      break;
    case 'h':
      if (*++format == 'h')
        {
          ++format;
          specinfo->arg_type = CHAR_ARG;
        }
      else
        specinfo->arg_type = SHORT_ARG;
      break;
    case 'l':
      if (*++format == 'l')
        {
          ++format;
          specinfo->arg_type = LONG_LONG_ARG;
          break;
        }
      else
        {
          specinfo->arg_type = LONG_ARG;
          break;
        }
    case 'q':
      ++format;
      specinfo->arg_type = LONG_LONG_ARG;
      break;
#ifdef HAVE_STDINT_H
    case 'j':
      ++format;
      specinfo->arg_type = INTMAX_ARG;
      break;
#endif
    case 'z':
      ++format;
      specinfo->arg_type = SIZE_ARG;
      break;
    case 't':
      ++format;
      specinfo->arg_type = PTRDIFF_ARG;
      break;
    case 'L':
      ++format;
      specinfo->arg_type = LONG_DOUBLE_ARG;
      break;
    case 'F':
      ++format;
      specinfo->arg_type = MPF_ARG;
      break;
    case 'Q':
      ++format;
      specinfo->arg_type = MPQ_ARG;
      break;
    case 'M':
      ++format;
      specinfo->arg_type = MP_LIMB_ARG;
      break;
    case 'N':
      ++format;
      specinfo->arg_type = MP_LIMB_ARRAY_ARG;
      break;
    case 'Z':
      ++format;
      specinfo->arg_type = MPZ_ARG;
      break;
      /* mpfr specific specifiers */
    case 'P':
      ++format;
      specinfo->arg_type = MPFR_PREC_ARG;
      break;
    case 'R':
      ++format;
      specinfo->arg_type = MPFR_ARG;
    }
  return format;
}

/* some macros and functions for filling the buffer */
/* CONSUME_VA_ARG removes from va_list AP the type expected by SPECINFO */

#ifdef HAVE_STDINT_H
#define CASE_INTMAX_ARG							\
      case INTMAX_ARG:							\
	(void) va_arg ((ap), intmax_t);					\
	break;								
#endif

#define CONSUME_VA_ARG(specinfo, ap)					\
  do {									\
    switch ((specinfo).arg_type)					\
      {									\
      case CHAR_ARG:							\
      case SHORT_ARG:							\
	(void) va_arg ((ap), int);					\
	break;								\
      case LONG_ARG:							\
	if (((specinfo).spec == 'd') || ((specinfo).spec == 'i')	\
	    || ((specinfo).spec == 'o') || ((specinfo).spec == 'u')	\
	    || ((specinfo).spec == 'x') || ((specinfo).spec == 'X'))	\
	  (void) va_arg ((ap), long);					\
	else if ((specinfo).spec == 'c')				\
	  (void) va_arg ((ap), wint_t);					\
	else if ((specinfo).spec == 's')				\
	  (void) va_arg ((ap), wchar_t);				\
	break;								\
      case LONG_LONG_ARG:						\
	(void) va_arg ((ap), long long);				\
	break;								\
      CASE_INTMAX_ARG							\
      case SIZE_ARG:							\
	(void) va_arg ((ap), size_t);					\
	break;								\
      case PTRDIFF_ARG:							\
	(void) va_arg ((ap), ptrdiff_t);				\
	break;								\
      case LONG_DOUBLE_ARG:						\
	(void) va_arg ((ap), long double);				\
	break;								\
      case MPF_ARG:							\
	(void) va_arg ((ap), mpf_srcptr);				\
	break;								\
      case MPQ_ARG:							\
	(void) va_arg ((ap), mpq_srcptr);				\
	break;								\
      case MP_LIMB_ARG:							\
      case MP_LIMB_ARRAY_ARG:						\
	(void) va_arg ((ap), mp_ptr);					\
	(void) va_arg ((ap), mp_size_t);				\
	break;								\
      case MPZ_ARG:							\
	(void) va_arg ((ap), mpz_srcptr);				\
	break;								\
      default:								\
	switch ((specinfo).spec)					\
	  {								\
	  case 'd':							\
	  case 'i':							\
	  case 'o':							\
	  case 'u':							\
	  case 'x':							\
	  case 'X':							\
	  case 'c':							\
	    (void) va_arg ((ap), int);					\
	    break;							\
	  case 'f':							\
	  case 'F':							\
	  case 'e':							\
	  case 'E':							\
	  case 'g':							\
	  case 'G':							\
	  case 'a':							\
	  case 'A':							\
	    (void) va_arg ((ap), double);				\
	    break;							\
	  case 's':							\
	    (void) va_arg ((ap), char *);				\
	    break;							\
	  case 'p':							\
	    (void) va_arg ((ap), void *);				\
	  }								\
      }									\
  } while (0)

/* process the format part which does not deal with mpfr types */
#define FLUSH(flag, start, end, ap, buf_ptr)			\
  do {								\
    __gmp_const size_t n = (end) - (start);			\
    if ((flag))							\
      /* previous specifiers are understood by gmp_printf */	\
      {								\
	char fmt_copy[n + 1];					\
	strncpy (fmt_copy, (start), n);				\
	fmt_copy[n] = '\0';					\
	sprntf_gmp ((buf_ptr), (fmt_copy), (ap));		\
	(flag) = 0;						\
      }								\
    else if ((start) != (end))					\
      /* no conversion specification, just simple characters */	\
      buffer_cat ((buf_ptr), (start), n);			\
  } while (0)

struct string_buffer
{
  char *start;                  /* beginning of the buffer */
  char *curr;                   /* last character (!= '\0') written */
  size_t size;                  /* buffer capacity */
};

static void
buffer_init (struct string_buffer *b, size_t s)
{
  b->start = (char *) mpfr_default_allocate (s);
  b->start[0] = '\0';
  b->curr = b->start;
  b->size = s;
}

/* Concatenate the string to the buffer and expand it if needed. */
static void
buffer_cat (struct string_buffer *b, __gmp_const char *s, size_t l)
{
  if (MPFR_UNLIKELY ((b->curr + l + 1) > (b->start + b->size)))
    {
      __gmp_const ptrdiff_t pos = b->curr - b->start;
      __gmp_const size_t n =
        ((l + 1 >
          MAX_CHAR_PRODUCED_BY_SPEC * sizeof (char)) ? l +
         1 : MAX_CHAR_PRODUCED_BY_SPEC) * sizeof (char);
      b->start =
        (char *) mpfr_default_reallocate (b->start, b->size, b->size + n);
      b->size += n;
      b->curr = b->start + pos;
    }
  strncat (b->curr, s, l);
  b->curr += l;
}

/* let gmp_xprintf process the part it can understand */
static void
sprntf_gmp (struct string_buffer *b, __gmp_const char *fmt, va_list ap)
{
  int l;
  char *s;
  l = gmp_vasprintf (&s, fmt, ap);
  buffer_cat (b, s, l);
  mpfr_default_free (s, l + 1);
}

/* internal function printing the mpfr_t P into the buffer BUF according to
   specification described by SPEC.
   Note: SPEC must be a integer conversion specification. */
static void
sprnt_int (struct string_buffer *buf, mpfr_srcptr p, struct printf_spec spec)
{
  char format[20];
  char *s;
  int f;
  mpz_t z;

  if (MPFR_UNLIKELY (MPFR_IS_NAN (p)))
    {
      buffer_cat (buf, spec.spec == 'X' ? "NAN" : "nan", 3);
      return ;
    }

  if (MPFR_UNLIKELY (MPFR_IS_INF (p)))
    {
      char inf_str[5];
      int neg = MPFR_SIGN (p) < 0;      /* 0 if positive, 1 if negative */
      if (spec.spec == 'X')
        strcpy (inf_str, neg ? "-INF" : "INF");
      else
        strcpy (inf_str, neg ? "-inf" : "inf");
      buffer_cat (buf, inf_str, neg + 3);
      return ;
    }

  mpz_init (z);
  mpfr_get_z (z, p, spec.rnd_mode);

  memset (format, '\0', 20);
  f = 0;
  format[f++] = '%';
  if (spec.left)
    format[f++] = '-';
  else if (spec.pad == '0')
    format[f++] = '0';
  if (spec.showsign)
    format[f++] = '+';
  else if (spec.space)
    format[f++] = ' ';
  if (spec.alt)
    format[f++] = '#';
  format[f] = '\0';
  strncat (format, "*.*Z", 4);
  f += 4;
  format[f++] = spec.spec;
  format[f] = '\0';

  f = gmp_asprintf (&s, format, spec.width, spec.prec, z);

  buffer_cat (buf, s, f);

  mpfr_default_free (s, f + 1);
  mpz_clear (z);
}

/* sprnt_fp is an internal function printing the mpfr_t P into the buffer BUF
   according to specification described by SPEC.
   note: SPEC must be a floating point conversion specification. */
static void
sprnt_fp (struct string_buffer *buf, mpfr_srcptr p, struct printf_spec spec)
{
  __gmp_const char d_point[] = { MPFR_DECIMAL_POINT, '\0' };
  mp_exp_t exp;
  int base;
  int remove_trailing_zeros = 0;
  int nsd;                      /* number of significant digits */
  int str_len;
  char *str;
  char *str_curr = NULL;

  /* char_fp details how much characters are needed in each part of a float 
     print. */
  struct char_fp
  {
    int sgn;                    /* 1 if sign char is present */
    int base_prefix;            /* '0', 2 if '0x' or '0X' */
    int int_part;               /* Size of integral part given by get_str */
    int point;                  /* Decimal point char */
    int frac_part;              /* Size of fractional part given by get_str */
    int exp_part;               /* Size of exponent (always in base ten) */

    unsigned int total;         /* Total size */
  };

  struct char_fp nbc;

  if (MPFR_UNLIKELY (MPFR_IS_NAN (p)))
    {
      switch (spec.spec)
        {
        case 'a':
        case 'b':
        case 'e':
        case 'f':
        case 'g':
          buffer_cat (buf, "nan", 3);   /* [FIXME: padding !] */
          break;
        case 'A':
        case 'E':
        case 'F':
        case 'G':
          buffer_cat (buf, "NAN", 3);
        }
      return ;
    }

  if (MPFR_UNLIKELY (MPFR_IS_INF (p)))
    {
      int neg = MPFR_SIGN (p) < 0;      /* 0 if positive, 1 if negative */
      switch (spec.spec)
        {
        case 'a':
        case 'b':
        case 'e':
        case 'f':
        case 'g':
          buffer_cat (buf, neg ? "-inf" : "inf", neg + 3);
          break;
        case 'A':
        case 'E':
        case 'F':
        case 'G':
          buffer_cat (buf, neg ? "-INF" : "INF", neg + 3);
        }
      return ;
    }

  nbc.sgn = (MPFR_SIGN (p) < 0) || spec.showsign || spec.space ? 1 : 0;
  nbc.total = nbc.sgn;

  if ((spec.spec == 'g') || (spec.spec == 'G'))
    {
      mp_exp_t threshold =
        (spec.prec < 0) ? 6 : (spec.prec == 0) ? 1 : spec.prec;
      if (MPFR_UNLIKELY (MPFR_IS_ZERO (p)))
	{
          spec.spec = (spec.spec == 'g') ? 'f' : 'F';
          spec.prec = 6;
	}
      else if ((threshold > MPFR_GET_EXP (p) - 1) && (MPFR_GET_EXP (p) > -4))
        {
          spec.spec = (spec.spec == 'g') ? 'f' : 'F';
          spec.prec = threshold - MPFR_GET_EXP (p);
        }
      else
        {
          spec.spec = (spec.spec == 'g') ? 'e' : 'E';
          spec.prec = threshold - 1;
        }
      remove_trailing_zeros = spec.alt ? 0 : 1;
    }

  switch (spec.spec)
    {
    case 'a':
    case 'A':
      nbc.base_prefix = 2;
      nbc.int_part = 1;
      nbc.point = (spec.prec == 0) && (spec.alt == 0) ? 0 : 1;
      nbc.frac_part = spec.prec;
      if (MPFR_UNLIKELY (MPFR_IS_ZERO (p)))
	{
	  nbc.exp_part = 3;
	}
      else
	{
	  exp = (mp_exp_t) abs (MPFR_GET_EXP (p));
	  nbc.exp_part = ceil (log10 (exp)) + 2;
	  if (nbc.exp_part < 3)
	    /* the exponent always contains at least on digit in hexadecimal */
	    nbc.exp_part = 3;
	}
      nbc.total += nbc.base_prefix + nbc.int_part + nbc.point + nbc.exp_part
        + ((spec.prec < 0) ? 0 : nbc.frac_part);

      base = 16;
      /* If precision is missing, we ask for (use)full precision */
      nsd = (spec.prec < 0) ? 0 : nbc.int_part + nbc.frac_part;
      if (nsd == 1)
        /* mpfr_get_str() do not allow asking for 1 digit in 2^n bases */
        nsd = 2;
      remove_trailing_zeros = (spec.prec < 0);
      break;

    case 'b':
      nbc.base_prefix = 0;

      nbc.int_part = 1;
      nbc.point = (spec.prec == 0) && (spec.alt == 0) ? 0 : 1;
      nbc.frac_part = spec.prec;
      if (MPFR_UNLIKELY (MPFR_IS_ZERO (p)))
	{
	  nbc.exp_part = 3;
	}
      else
	{
	  exp = (mp_exp_t) abs (MPFR_GET_EXP (p));
	  nbc.exp_part = ceil (log10 (exp)) + 2;
	  if (nbc.exp_part < 3)
	    /* the exponent always contains at least one digit in base 2 */
	    nbc.exp_part = 3;
	}
      nbc.total += nbc.base_prefix + nbc.int_part + nbc.point + nbc.exp_part
        + ((spec.prec < 0) ? 0 : nbc.frac_part);

      base = 2;
      /* If precision is missing, we ask for (use)full precision */
      nsd = (spec.prec < 0) ? 0 : nbc.int_part + nbc.frac_part;
      if (nsd == 1)
        /* mpfr_get_str() do not allow asking for 1 digit in 2^n bases */
        nsd = 2;
      remove_trailing_zeros = (spec.prec < 0);
      break;

    case 'f':
    case 'F':
      nbc.base_prefix = 0;
      nbc.exp_part = 0;
      if (MPFR_UNLIKELY (MPFR_IS_ZERO (p)))
	{
	  nbc.int_part = 1;
	}
      else
	{
	  mpfr_t l;
	  mpfr_init (l);
	  mpfr_abs (l, p, GMP_RNDN);
	  mpfr_log10 (l, l, GMP_RNDN);
	  nbc.int_part = mpfr_get_si (l, GMP_RNDU);
	  if (nbc.int_part < 0)
	    nbc.int_part = 1;
	  mpfr_clear (l);
	}
      nbc.point = (spec.prec == 0) && (spec.alt == 0) ? 0 : 1;
      nbc.frac_part = (spec.prec < 0) ? 6 : spec.prec;  /* [FIXME] */
      nbc.total += nbc.int_part + nbc.point + nbc.frac_part;

      base = 10;
      nsd = nbc.int_part + nbc.frac_part;
      break;

    case 'e':
    case 'E':
    default:
      nbc.base_prefix = 0;

      nbc.int_part = 1;
      nbc.point = (spec.prec == 0) && (spec.alt == 0) ? 0 : 1;
      nbc.frac_part = spec.prec;
      if (MPFR_UNLIKELY (MPFR_IS_ZERO (p)))
	{
	  nbc.exp_part = 4;
	}
      else
	{
	  exp = (mp_exp_t) abs (MPFR_GET_EXP (p));
	  nbc.exp_part = ceil (log10 (exp)) + 2;
	  if (nbc.exp_part < 4)
	    /* the exponent always contains at least two digits in base ten */
	    nbc.exp_part = 4;
	}
      nbc.total += nbc.base_prefix + nbc.int_part + nbc.point + nbc.exp_part
        + ((spec.prec < 0) ? 0 : nbc.frac_part);

      base = 10;
      nsd = (spec.prec < 0) ? 0 : nbc.int_part + nbc.frac_part;
      if (spec.prec < 0)
        remove_trailing_zeros = 1;
    }

  /* get NSD significant digits from mpfr */
  str = mpfr_get_str (0, &exp, base, nsd, p, spec.rnd_mode);
  MPFR_ASSERTN (str != NULL);

  if (spec.spec == 'A')
    {
      char *s1 = str;
      while (*s1)
        {
          *s1 = toupper (*s1);
          s1++;
        }
    }

  str_len = strlen (str);
  str_curr = str;
  if (spec.prec < 0)
    spec.prec = str_len - nbc.int_part - nbc.sgn;

  if (nbc.frac_part < 0)
    {
      nbc.frac_part = str_len - nbc.int_part - nbc.sgn;
      if (remove_trailing_zeros)
        {
          int str_frac_len = str_len - 1;
          while (str_frac_len > nbc.int_part + nbc.sgn)
            {
              if (str_curr[str_frac_len] != '0')
                break;
              --str_frac_len;
              --nbc.frac_part;
            }
        }
      nbc.total += nbc.frac_part;
    }
  MPFR_ASSERTD (nbc.total < MAX_CHAR_PRODUCED_BY_SPEC);

  /* right justification with spaces */
  if ((spec.left == 0) && (spec.pad == ' ') && (nbc.total < spec.width))
    {
      __gmp_const int n = spec.width - nbc.total;
      char *padding = (char *) mpfr_default_allocate (n + 1);
      memset (padding, ' ', n);
      padding[n] = '\0';
      buffer_cat (buf, padding, n);
      mpfr_default_free (padding, n + 1);
    }

  /* build the string */
  if (nbc.sgn)
    /* sign character */
    {
      buffer_cat (buf, (MPFR_SIGN (p) < 0) ? "-" : "+", 1);
      if (MPFR_SIGN (p) < 0)
        str_curr++;
    }
  if ((spec.left == 0) && (spec.pad == '0') && (nbc.total < spec.width))
    /* leading zeros in integral part */
    {
      __gmp_const int n = spec.width - nbc.total;
      char *padding = (char *) mpfr_default_allocate (n + 1);
      memset (padding, '0', n);
      padding[n] = '\0';
      buffer_cat (buf, padding, n);
      mpfr_default_free (padding, n + 1);
    }
  /* integer part */
  if (exp < 0 && (spec.spec == 'f' || spec.spec == 'F'))
    /* there is always a digit before the decimal point */
    buffer_cat (buf, "0", 1);
  else
    {
      buffer_cat (buf, str_curr, nbc.int_part);
      str_curr += nbc.int_part;
    }

  if (nbc.point)
    /* decimal point */
    buffer_cat (buf, d_point, 1);

  if (nbc.frac_part)
    /* fractionnal part */
    {
      if ((spec.spec == 'f' || spec.spec == 'F') && (exp < 0))
        /* leading zeros in fractional part when p < 1 */
        {
          char *zeros = (char *) mpfr_default_allocate (1 - exp);
          memset (zeros, '0', -exp);
          zeros[-exp] = '\0';
          buffer_cat (buf, zeros, -exp);
          mpfr_default_free (zeros, -exp + 1);
        }
      buffer_cat (buf, str_curr, nbc.frac_part);

      if ((remove_trailing_zeros == 0) && (nbc.frac_part < spec.prec))
        /* add trailing zeros */
        {
          __gmp_const int n = spec.spec - nbc.frac_part;
          char *zeros = (char *) mpfr_default_allocate (n + 1);
          memset (zeros, '0', n);
          zeros[n] = '\0';
          buffer_cat (buf, zeros, n);
          mpfr_default_free (zeros, n + 1);
        }
    }

  if (nbc.exp_part)
    /* exponent part */
    {
      char exp_fmt[6];
      char *exp_str;
      switch (spec.spec)
        {
        case 'a':
        case 'A':
          buffer_cat (buf, spec.spec == 'A' ? "P" : "p", 1);
          exp = MPFR_IS_ZERO (p)? 0: (exp - 1) * 4;
	  strcpy (exp_fmt, "%+.1"_MP_EXP_FORMAT_SPEC);
          break;
        case 'b':
          buffer_cat (buf, "p", 1);
          exp = MPFR_IS_ZERO (p)? 0: exp - 1;
	  strcpy (exp_fmt, "%+.1"_MP_EXP_FORMAT_SPEC);
          break;
        case 'e':
        case 'E':
        case 'f':
        case 'F':
          buffer_cat (buf,
                      ((spec.spec == 'E') || (spec.spec == 'F')) ? "E" : "e",
                      1);
          exp = MPFR_IS_ZERO (p)? 0: exp - nbc.int_part;
	  strcpy (exp_fmt, "%+.2"_MP_EXP_FORMAT_SPEC);
        }

      MPFR_ASSERTN (exp - 1 >= LONG_MIN);
      MPFR_ASSERTN (exp - 1 <= LONG_MAX);
      exp_str = (char *) mpfr_default_allocate (nbc.exp_part + 1);
      snprintf (exp_str, nbc.exp_part, exp_fmt, exp);

      MPFR_ASSERTD (strlen (exp_str) == nbc.exp_part - 1);
      buffer_cat (buf, exp_str, nbc.exp_part - 1);
      mpfr_default_free (exp_str, nbc.exp_part + 1);
    }
  /* left justification with spaces */
  if (spec.left && (spec.pad == ' ') && (nbc.total < spec.width))
    {
      __gmp_const int n = spec.width - nbc.total;
      char *padding = (char *) mpfr_default_allocate (n + 1);
      memset (padding, ' ', n);
      padding[n] = '\0';
      buffer_cat (buf, padding, n);
      mpfr_default_free (padding, n + 1);
    }

  mpfr_free_str (str);
}


int
mpfr_vasprintf (char **ptr, __gmp_const char *fmt, va_list ap)
{
  struct string_buffer buf;
  ptrdiff_t nbchar;

  /* informations on the conversion specification filled by the parser */
  struct printf_spec spec;
  /* flag raised when previous part of fmt need to be processed by 
     gmp_vsnprintf */
  int gmp_fmt_flag;
  /* beginning and end of the previous unprocessed part of fmt */
  __gmp_const char *start, *end;
  /* pointer to arguments for gmp_vasprintf */
  va_list ap2;

  nbchar = 0;
  buffer_init (&buf, MAX_CHAR_PRODUCED_BY_SPEC * sizeof (char));
  gmp_fmt_flag = 0;
  va_copy (ap2, ap);
  start = fmt;
  while (*fmt)
    {
      /* Look for the next format specification */
      while ((*fmt) && (*fmt != '%'))
        ++fmt;

      if (*fmt == '\0')
        break;

      if (*++fmt == '%')
        continue;

      end = fmt - 1;

      /* format string analysis */
      specinfo_init (&spec);
      fmt = parse_flags (fmt, &spec);

      READ_INT (fmt, spec, width, width_analysis);
    width_analysis:
      if (spec.width < 0)
        {
          spec.left = 1;
          spec.width = -spec.width;
        }
      if (*fmt == '.')
        {
          __gmp_const char *f = ++fmt;
          READ_INT (fmt, spec, prec, prec_analysis);
        prec_analysis:
          if (f == fmt)
            spec.prec = -1;
        }
      else
        spec.prec = -1;

      fmt = parse_arg_type (fmt, &spec);
      if (spec.arg_type == MPFR_ARG)
        {
          switch (*fmt)
            {
            case '\0':
              break;
            case '*':
              ++fmt;
              spec.rnd_mode = va_arg (ap, mp_rnd_t);
              break;
            case 'D':
              ++fmt;
              spec.rnd_mode = GMP_RNDD;
              break;
            case 'U':
              ++fmt;
              spec.rnd_mode = GMP_RNDU;
              break;
            case 'Z':
              ++fmt;
              spec.rnd_mode = GMP_RNDZ;
              break;
            case 'N':
              ++fmt;
            default:
              spec.rnd_mode = GMP_RNDN;
            }
        }

      spec.spec = *fmt;
      if (*fmt)
        fmt++;

      /* Format processing */
      if (spec.spec == '\0')
        break;
      else if (spec.spec == 'n')
        {
          int *int_ptr;         /* [FIXME] pas forcement un entier (eventuellement un mp{zqf}) */
          int_ptr = va_arg (ap, int *);
          FLUSH (gmp_fmt_flag, start, end, ap2, &buf);
          va_end (ap2);
          va_copy (ap2, ap);
          start = fmt;
          *int_ptr = (int) (buf.curr - buf.start);
        }
      else if (spec.arg_type == MPFR_PREC_ARG)
        {
          char *s;
          int l;
          mpfr_prec_t prec;
          prec = va_arg (ap, mpfr_prec_t);

          FLUSH (gmp_fmt_flag, start, end, ap2, &buf);
          va_end (ap2);
          va_copy (ap2, ap);
          start = fmt;

          l = gmp_asprintf (&s, "%"_MPFR_PREC_FORMAT_SPEC, prec);
          buffer_cat (&buf, s, l);
          mpfr_default_free (s, l + 1);
        }
      else if (spec.arg_type == MPFR_ARG)
        {
          mpfr_srcptr p;
          p = va_arg (ap, mpfr_srcptr);

          FLUSH (gmp_fmt_flag, start, end, ap2, &buf);
          va_end (ap2);
          va_copy (ap2, ap);
          start = fmt;

          switch (spec.spec)
            {
            case 'd':
            case 'i':
            case 'o':
            case 'u':
            case 'x':
            case 'X':
              sprnt_int (&buf, p, spec);
              break;
            case 'f':
            case 'F':
            case 'e':
            case 'E':
            case 'g':
            case 'G':
            case 'a':
            case 'A':
	    case 'b':  
              sprnt_fp (&buf, p, spec);
            }
        }
      else
        {
          CONSUME_VA_ARG (spec, ap);
          gmp_fmt_flag = 1;
        }
    }

  if (start != fmt)
    FLUSH (gmp_fmt_flag, start, fmt, ap2, &buf);

  va_end (ap2);
  nbchar = buf.curr - buf.start;
  buf.start =
    (char *) mpfr_default_reallocate (buf.start, buf.size, nbchar + 1);
  *ptr = buf.start;
  return nbchar;
}

#endif /* HAVE_STDARG */
