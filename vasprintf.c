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

/* The mpfr_printf-like functions are defined only if stdarg.h exists */
#ifdef HAVE_STDARG

#include <stdarg.h>

#ifdef HAVE_WCHAR_H
#include <wchar.h>
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif

#include <stddef.h>             /* for ptrdiff_t */
#include <string.h>             /* for strlen */

#include "mpfr-impl.h"

#if   _MPFR_PREC_FORMAT == 1
#define MPFR_PREC_FORMAT_SPEC "hu"
#elif _MPFR_PREC_FORMAT == 2
#define MPFR_PREC_FORMAT_SPEC "u"
#elif _MPFR_PREC_FORMAT == 3
#define MPFR_PREC_FORMAT_SPEC "lu"
#else
#error "mpfr_prec_t size not supported"
#endif

#if (__GMP_MP_SIZE_T_INT == 1)
#define MPFR_EXP_FORMAT_SPEC "i"
#elif (__GMP_MP_SIZE_T_INT == 0)
#define MPFR_EXP_FORMAT_SPEC "li"
#else
#error "mp_exp_t size not supported"
#endif

/* Output for special values defined in the C99 standard */
#define MPFR_NAN_STRING_LC "nan"
#define MPFR_NAN_STRING_UC "NAN"
#define MPFR_INF_STRING_LC "inf"
#define MPFR_INF_STRING_UC "INF"

/* We assume that a single conversion specifier produces at most 4095 chars
   (Rationale for International Standard -Programming Languages- C
   Revision 5.10 April-2003, 7.19.6.1 p.152). */
#define MAX_CHAR_PRODUCED_BY_SPEC 4096

static const char num_to_text[16] = "0123456789abcdef";

/* some macro and functions for parsing format string */
#define READ_INT(format, specinfo, field, label_out)            \
  do {                                                          \
    while (*(format))                                           \
      {                                                         \
        switch (*(format))                                      \
          {                                                     \
          case '0':                                             \
          case '1':                                             \
          case '2':                                             \
          case '3':                                             \
          case '4':                                             \
          case '5':                                             \
          case '6':                                             \
          case '7':                                             \
          case '8':                                             \
          case '9':                                             \
            MPFR_ASSERTN ((specinfo).field < INT_MAX / 10);     \
            (specinfo).field *= 10;                             \
            (specinfo).field += *(format) - '0';                \
            ++(format);                                         \
            break;                                              \
          case '*':                                             \
            (specinfo).field = va_arg ((ap), int);              \
            ++(format);                                         \
          default:                                              \
            goto label_out;                                     \
          }                                                     \
      }                                                         \
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
  unsigned int group:1;         /* ' gnu flag (not for gmp/mpfr types) */

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

static const char *
parse_flags (const char *format, struct printf_spec *specinfo)
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
          /* gnu extension for thousand separator */
          ++format;
          break;
        default:
          return format;
        }
    }
  return format;
}

static const char *
parse_arg_type (const char *format, struct printf_spec *specinfo)
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


/* some macros and functions filling the buffer */

/* CONSUME_VA_ARG removes from va_list AP the type expected by SPECINFO */
#ifdef HAVE_STDINT_H
#define CASE_INTMAX_ARG(specinfo, ap)           \
  case INTMAX_ARG:                              \
  (void) va_arg ((ap), intmax_t);               \
  break;
#endif

#ifdef HAVE_WCHAR_H
#define CASE_LONG_ARG(specinfo, ap)                             \
  case LONG_ARG:                                                \
  if (((specinfo).spec == 'd') || ((specinfo).spec == 'i')      \
      || ((specinfo).spec == 'o') || ((specinfo).spec == 'u')   \
      || ((specinfo).spec == 'x') || ((specinfo).spec == 'X'))  \
    (void) va_arg ((ap), long);                                 \
  else if ((specinfo).spec == 'c')                              \
    (void) va_arg ((ap), wint_t);                               \
  else if ((specinfo).spec == 's')                              \
    (void) va_arg ((ap), wchar_t);                              \
  break;
#else
#define CASE_LONG_ARG(specinfo, ap)             \
  case LONG_ARG:                                \
  (void) va_arg ((ap), long);                   \
  break;
#endif

#define CONSUME_VA_ARG(specinfo, ap)            \
  do {                                          \
    switch ((specinfo).arg_type)                \
      {                                         \
      case CHAR_ARG:                            \
      case SHORT_ARG:                           \
        (void) va_arg ((ap), int);              \
        break;                                  \
        CASE_LONG_ARG (specinfo, ap)            \
      case LONG_LONG_ARG:                       \
        (void) va_arg ((ap), long long);        \
        break;                                  \
        CASE_INTMAX_ARG (specinfo, ap)          \
      case SIZE_ARG:                            \
        (void) va_arg ((ap), size_t);           \
        break;                                  \
      case PTRDIFF_ARG:                         \
        (void) va_arg ((ap), ptrdiff_t);        \
        break;                                  \
      case LONG_DOUBLE_ARG:                     \
        (void) va_arg ((ap), long double);      \
        break;                                  \
      case MPF_ARG:                             \
        (void) va_arg ((ap), mpf_srcptr);       \
        break;                                  \
      case MPQ_ARG:                             \
        (void) va_arg ((ap), mpq_srcptr);       \
        break;                                  \
      case MP_LIMB_ARG:                         \
      case MP_LIMB_ARRAY_ARG:                   \
        (void) va_arg ((ap), mp_ptr);           \
        (void) va_arg ((ap), mp_size_t);        \
        break;                                  \
      case MPZ_ARG:                             \
        (void) va_arg ((ap), mpz_srcptr);       \
        break;                                  \
      default:                                  \
        switch ((specinfo).spec)                \
          {                                     \
          case 'd':                             \
          case 'i':                             \
          case 'o':                             \
          case 'u':                             \
          case 'x':                             \
          case 'X':                             \
          case 'c':                             \
            (void) va_arg ((ap), int);          \
            break;                              \
          case 'f':                             \
          case 'F':                             \
          case 'e':                             \
          case 'E':                             \
          case 'g':                             \
          case 'G':                             \
          case 'a':                             \
          case 'A':                             \
            (void) va_arg ((ap), double);       \
            break;                              \
          case 's':                             \
            (void) va_arg ((ap), char *);       \
            break;                              \
          case 'p':                             \
            (void) va_arg ((ap), void *);       \
          }                                     \
      }                                         \
  } while (0)

/* process the format part which does not deal with mpfr types */
#define FLUSH(flag, start, end, ap, buf_ptr)                    \
  do {                                                          \
    const size_t n = (end) - (start);                           \
    if ((flag))                                                 \
      /* previous specifiers are understood by gmp_printf */    \
      {                                                         \
        char fmt_copy[n + 1];                                   \
        strncpy (fmt_copy, (start), n);                         \
        fmt_copy[n] = '\0';                                     \
        sprntf_gmp ((buf_ptr), (fmt_copy), (ap));               \
        (flag) = 0;                                             \
      }                                                         \
    else if ((start) != (end))                                  \
      /* no conversion specification, just simple characters */ \
      buffer_cat ((buf_ptr), (start), n);                       \
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
  b->start = (char *) (*__gmp_allocate_func) (s);
  b->start[0] = '\0';
  b->curr = b->start;
  b->size = s;
}

/* Concatenate the string to the buffer and expand it if needed. */
static void
buffer_cat (struct string_buffer *b, const char *s, size_t len)
{
  MPFR_ASSERTD (len <= strlen (s));
  if (MPFR_UNLIKELY ((b->curr + len + 1) > (b->start + b->size)))
    {
      const size_t pos = b->curr - b->start;
      const size_t n = sizeof (char)
        * ((len + 1 > MAX_CHAR_PRODUCED_BY_SPEC) ?
           len + 1 : MAX_CHAR_PRODUCED_BY_SPEC);
      b->start =
        (char *) (*__gmp_reallocate_func) (b->start, b->size, b->size + n);
      b->size += n;
      b->curr = b->start + pos;
    }
  strncat (b->curr, s, len);
  b->curr += len;
}

/* Add N characters C to the end of buffer */
static void
buffer_pad (struct string_buffer *b, const char c, const size_t n)
{
  char *padding;

  padding = (char *) (*__gmp_allocate_func) (n + 1);
  memset (padding, c, n);
  padding[n] = '\0';
  buffer_cat (b, padding, n);
  (*__gmp_free_func) (padding, n + 1);
}

/* Print NaN with padding */
static void
sprnt_nan (struct string_buffer *buf, const struct printf_spec spec)
{
  /* right justification padding */
  if ((spec.left == 0) && (spec.width > 3))
    buffer_pad (buf, ' ', spec.width - 3);

  switch (spec.spec)
    {
    case 'A':
    case 'E':
    case 'F':
    case 'G':
    case 'X':
      buffer_cat (buf, MPFR_NAN_STRING_UC, strlen (MPFR_NAN_STRING_UC));
      break;

    default:
      buffer_cat (buf, MPFR_NAN_STRING_LC, strlen (MPFR_NAN_STRING_LC));
    }

  /* left justification padding */
  if ((spec.left == 1) && (spec.width > 3))
    buffer_pad (buf, ' ', spec.width - 3);
  return;
}

/* Print Infinities with padding
   NEG = -1 for '-inf', 0 for 'inf' */
static void
sprnt_inf (struct string_buffer *buf, const struct printf_spec spec, int neg)
{
  /* right justification padding */
  if ((spec.left == 0) && (spec.width > 3 - neg))
    buffer_pad (buf, ' ', spec.width - 3 + neg);

  switch (spec.spec)
    {
    case 'A':
    case 'E':
    case 'F':
    case 'G':
    case 'X':
      if (neg < 0)
        buffer_cat (buf, "-", 1);
      buffer_cat (buf, MPFR_INF_STRING_UC, strlen (MPFR_INF_STRING_UC));
      break;
    default:
      if (neg < 0)
        buffer_cat (buf, "-", 1);
      buffer_cat (buf, MPFR_INF_STRING_LC, strlen (MPFR_INF_STRING_LC));
    }

  /* left justification padding */
  if ((spec.left == 1) && (spec.width > 3 - neg))
    buffer_pad (buf, ' ', spec.width - 3 + neg);
}

/* let gmp_xprintf process the part it can understand */
static void
sprntf_gmp (struct string_buffer *b, const char *fmt, va_list ap)
{
  int len;
  char *s;

  len = gmp_vasprintf (&s, fmt, ap);
  buffer_cat (b, s, len);
  mpfr_free_str (s);
}

/* internal function printing the mpfr_t P into the buffer BUF according to
   specification described by SPEC.
   Note: SPEC must be a integer conversion specification. */
static void
sprnt_int (struct string_buffer *buf, mpfr_srcptr p, struct printf_spec spec)
{
  char format[10];  /* buffer for the format string corresponding to spec
                       (see below) */
  char *s;
  int f;
  mpz_t z;

  if (MPFR_UNLIKELY (MPFR_IS_NAN (p)))
    {
      sprnt_nan (buf, spec);
      return;
    }

  if (MPFR_UNLIKELY (MPFR_IS_INF (p)))
    {
      sprnt_inf (buf, spec, MPFR_SIGN (p) < 0 ? -1 : 0);
      return;
    }

  mpz_init (z);
  mpfr_get_z (z, p, spec.rnd_mode);

  /* format contains at most 9 characters plus the terminating '\0'
     like in "%-+#*.*Zo" */
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
  MPFR_ASSERTD (f < 10);

  f = gmp_asprintf (&s, format, spec.width, spec.prec, z);

  buffer_cat (buf, s, f);

  mpfr_free_str (s);
  mpz_clear (z);
}

/* return ceil(log(|MPFR_EXP(p)|)/log(10)) + 2 assuming the exponent >= 1
   for exp_part. FIXME: The exponent that will actually be written is not
   the exponent of the MPFR number! */
static unsigned long
uceil_log10_exp_p2 (mpfr_srcptr p)
{
  unsigned long c = 2;  /* 1 for exponent char + 1 for sign char */
  mp_exp_t e = MPFR_GET_EXP (p);
  mp_exp_unsigned_t x;

  x = SAFE_ABS (mp_exp_unsigned_t, e);
  while (x > 1)
    {
      c++;
      x = (x + 9) / 10;
    }
  return c;
}

/* sprnt_fp_a prints a mpfr_t in "%a" and "%A" cases */
static void
sprnt_fp_a (struct string_buffer *buf, mpfr_srcptr p, struct printf_spec spec)
{
  char *str;
  mp_exp_t exp;
  const char d_point[] = { MPFR_DECIMAL_POINT, '\0' };

  /* char_fp details how much characters are needed in each part of a float
     print. */
  struct char_fp
  {
    int sgn;                    /* 1 if sign char is present */
    int point;                  /* 1 if decimal point char   */
    long int_part;              /* Size of integral part     */
    long frac_part;             /* Size of fractional part   */
    unsigned long exp_part;     /* Size of exponent (always in base ten) */

    unsigned int total;         /* Total size: sum of the precedent fields
                                   plus the base prefix plus the first digit
                                   before the decimal point */
  };

  struct char_fp nbc;

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (p)))
    {
      if (MPFR_IS_NAN (p))
        {
          if (spec.spec == 'A')
            {
              nbc.total = strlen (MPFR_NAN_STRING_UC);
              str = (char *) (*__gmp_allocate_func) (1 + nbc.total);
              str[0] = '\0';
              strcat (str, MPFR_NAN_STRING_UC);
            }
          else
            {
              nbc.total = strlen (MPFR_NAN_STRING_LC);
              str = (char *) (*__gmp_allocate_func) (1 + nbc.total);
              str[0] = '\0';
              strcat (str, MPFR_NAN_STRING_LC);
            }
        }
      else if (MPFR_IS_INF (p))
        {
          if (spec.spec == 'A')
            {
              nbc.total = strlen (MPFR_INF_STRING_UC);
              nbc.total += (MPFR_SIGN (p) < 0) ? 1 : 0;
              str = (char *) (*__gmp_allocate_func) (1 + nbc.total);
              str[0] = '\0';
              if (MPFR_SIGN (p) < 0)
                strcat (str, "-");
              strcat (str, MPFR_INF_STRING_UC);
            }
          else
            {
              nbc.total = strlen (MPFR_INF_STRING_LC);
              nbc.total += (MPFR_SIGN (p) < 0) ? 1 : 0;
              str = (char *) (*__gmp_allocate_func) (1 + nbc.total);
              str[0] = '\0';
              if (MPFR_SIGN (p) < 0)
                strcat (str, "-");
              strcat (str, MPFR_INF_STRING_LC);
            }
        }
      else
        /* p == 0 */
        {
          /* First, we compute how much characters will be displayed */
          nbc.sgn =
            (MPFR_SIGN (p) < 0) || spec.showsign || spec.space ? 1 : 0;
          nbc.frac_part = (spec.prec < 0) ? 0 : spec.prec;
          nbc.point = (nbc.frac_part == 0) && (spec.alt == 0) ? 0 : 1;
          nbc.exp_part = 3;
          nbc.total = nbc.sgn + nbc.point + nbc.frac_part + nbc.exp_part;
          nbc.total += 2; /* 2 characters "0x" or "0X" */
          if ((spec.left == 0) && (spec.pad == '0')
              && (nbc.total < spec.width))
            /* pad with leading zeros */
            nbc.int_part =
              (spec.width > nbc.total) ? 1 : spec.width - nbc.total;
          else
            /* we only need the character '0' */
            nbc.int_part = 1;
          nbc.total += nbc.int_part;

          str = (char *)(*__gmp_allocate_func) (nbc.total + 1);

          /* Now, we build the string */
          if (nbc.sgn != 0)
            /* signed zero */
            {
              str[0] =
                (MPFR_SIGN (p) < 0) ? '-' : (spec.showsign) ? '+' : ' ';
              str[1] = '\0';
            }
          else
            /* positive zero, no sign displayed */
            str[0] = '\0';
          /* prefix */
          strcat (str, (spec.spec == 'A') ? "0X" : "0x");
          /* integral part, with optional padding */
          {
            int i;
            for (i = 1; i < nbc.int_part; ++i)
              strcat (str, "0");
          }
          /* digit */
          strcat (str, "0");
          if (nbc.point == 1)
            strcat (str, d_point);
          if (nbc.frac_part > 0)
            {
              /* fill frac_part with '0' */
              int i;
              for (i = 0; i < nbc.frac_part; i++)
                strcat (str, "0");
            }
          /* exponent */
          strcat (str, (spec.spec == 'A') ? "P+0" : "p+0");
          MPFR_ASSERTD (nbc.total == strlen (str));
        }
    }
  else
    /* p is a regular number, p != 0*/
    {
      char *raw_str;
      char *raw_str_cur; /* cursor in raw_str */

      /* get sign - and significant digits in raw_str */
      if (spec.prec < 0)
        /* Let mpfr_get_str determine precision. */
        {
          raw_str = mpfr_get_str (0, &exp, 16, 0, p, spec.rnd_mode);
          /* exp is the exponent for radix sixteen with decimal point BEFORE
             the first digit, we want the exponent for radix two and the decimal
             point AFTER the first digit */
          exp = (exp - 1) * 4;

        }
      else if (spec.prec > 0)
        /* One digit before decimal point plus SPEC.PREC after it. */
        {
          raw_str = mpfr_get_str (0, &exp, 16, 1 + spec.prec, p,
                                  spec.rnd_mode);
          /* exp is the exponent for radix sixteen with decimal point BEFORE
             the first digit, we want the exponent for radix two and the decimal
             point AFTER the first digit */
          exp = (exp - 1) * 4;
        }
      else
        /* One digit is sufficient but mpfr_get_str returns at least two
           digits. So, in order to avoid double rounding, we will build
           ourselves the raw string. */
        {
          mp_limb_t *pm = MPFR_MANT (p);
          mp_size_t ps;
          size_t n;
          int digit, shift;
          char digit_str[2] = {'\0', '\0'};
          int rnd_away;

          /* rnd_away:
             1 if round away from zero,
             0 if round to zero,
             -1 if not decided yet. */
          rnd_away =
            spec.rnd_mode == GMP_RNDD ? MPFR_IS_NEG (p) :
            spec.rnd_mode == GMP_RNDU ? MPFR_IS_POS (p) :
            spec.rnd_mode == GMP_RNDZ ? 0 : -1;

          /* Find out the radix-16 digit by grouping the 4 first digits. Even
             if MPFR_PREC (p) < 4, we can read 4 bits in its first limb */
          exp = MPFR_GET_EXP (p);
          shift = BITS_PER_MP_LIMB - 4;
          ps = (MPFR_PREC (p) - 1) / BITS_PER_MP_LIMB;
          digit = pm[ps]>>shift;

          /* exponent for radix-2 with the decimal point after the first
             hexadecimal digit */
          exp -= 4;

          if (MPFR_PREC (p) > 4)
            /* round taking into account bits outside the first f ones */
            {
              if (rnd_away == -1)
                /* Round to nearest mode: we have to decide in that particular
                   case if we have to round to away from zero or not */
                {
                  mp_limb_t limb, rb, mask;
                  /* compute rounding bit */
                  mask = MPFR_LIMB_ONE << (shift - 1);
                  rb = pm[ps] & mask;
                  if (rb == 0)
                    rnd_away = 0;
                  else
                    {
                      mask = MPFR_LIMB_MASK (shift);
                      limb = pm[ps] & mask;
                      while ((ps > 0) && (limb == 0))
                        limb = pm[--ps];
                      if (limb == 0)
                        /* tie case, round to even */
                        rnd_away = (digit & 1) ? 1 : 0;
                      else
                        rnd_away = 1;
                    }
                }

              if (rnd_away == 1)
                {
                  digit++;
                  if (digit > 15)
                    /* As we want only the first significant digit, we have
                       to shift one position to the left */
                    {
                      digit /= 16;
                      exp += 4;  /* no possible overflow because
                                    exp < MPFR_EXP_MAX/4 */
                    }
                }
            }

          n = 2 * sizeof (char);
          n += (MPFR_SIGN (p) < 0) ? 1 : 0;
          raw_str = (char *)(*__gmp_allocate_func) (n);
          raw_str[0] = '\0';

          if (MPFR_SIGN (p) < 0)
            strcat (raw_str, "-");

          MPFR_ASSERTD ((0 <= digit) && (digit <= 15));
          digit_str[0] = num_to_text [digit];

          strcat (raw_str, digit_str);
        }

      raw_str_cur = raw_str;
      if (spec.spec == 'A')
        /* All digits in upper case */
        {
          char *s1 = raw_str;
          while (*s1)
            {
              switch (*s1)
                {
                case 'a':
                  *s1 = 'A';
                  break;
                case 'b':
                  *s1 = 'B';
                  break;
                case 'c':
                  *s1 = 'C';
                  break;
                case 'd':
                  *s1 = 'D';
                  break;
                case 'e':
                  *s1 = 'E';
                  break;
                case 'f':
                  *s1 = 'F';
                  break;
                }
              s1++;
            }
        }

      /* compute the number of characters in each part of str */
      nbc.sgn = (MPFR_SIGN (p) < 0) || (spec.showsign) || (spec.space) ? 1 : 0;
      if (spec.prec < 0)
        {
          char *end;
          /* raw_str contains the sign if negative, the first digit that will
             be printed before the decimal point and all other digits are in
             the fractional part */
          nbc.frac_part = strlen (raw_str) - 1;
          nbc.frac_part -= (MPFR_SIGN (p) < 0) ? 1 : 0;
          /* we remove the trailing zeros if any */
          for (end = raw_str + strlen (raw_str) - 1;
               (*end == '0') && (nbc.frac_part > 0); --end)
            --nbc.frac_part;
        }
      else if (spec.prec == 0)
        nbc.frac_part = 0;
      else
        nbc.frac_part = spec.prec;
      nbc.point = (nbc.frac_part == 0) && (spec.alt == 0) ? 0 : 1;
      /* the exp_part contains the base character plus the sign character
         plus (exp == 0) ? 1 : floor(log10(abs(exp))) + 1 digits */
      nbc.exp_part = 3;
      {
        mp_exp_unsigned_t x;

        x = SAFE_ABS (mp_exp_unsigned_t, exp);
        while (x > 9)
          {
            nbc.exp_part++;
            x /= 10;
          }
      }
      /* Number of characters to be printed (2 for "0x") */
      nbc.total = 2 + nbc.sgn + nbc.point + nbc.frac_part
        + nbc.exp_part;
      /* optional heading zeros are counted in int_part */
      if ((spec.left == 0) && (spec.pad == '0') && (nbc.total < spec.width))
        nbc.int_part = spec.width - nbc.total;
      else
        nbc.int_part = 1;
      nbc.total += nbc.int_part;

      str = (char *) (*__gmp_allocate_func) (nbc.total + 1);
      str[0] = '\0';

      /* build str from raw_str */
      /* sign */
      if (nbc.sgn != 0)
        {
          if (MPFR_SIGN (p) < 0)
            {
              strcat (str, "-");
              ++raw_str_cur;
            }
          else if (spec.showsign)
            strcat (str, "+");
          else if (spec.space)
            strcat (str, " ");
        }
      /* base prefix */
      strcat (str, (spec.spec == 'A') ? "0X" : "0x");
      /* optional padding with heading zero */
      {
        long int i;
        for (i = nbc.int_part - 1; i > 0; i--)
          strcat (str, "0");
      }
      /* first digit */
      {
        char d[2];
        d[0] = *raw_str_cur++;
        d[1] = '\0';
        strcat (str, d);
      }

      if (nbc.frac_part == 0)
        /* no fractional part, optional decimal point */
        {
          if (spec.alt)
            strcat (str, d_point);
        }
      else
        /* a decimal point followed by fractional part */
        {
          char c;

          strcat (str, d_point);
          /* Some trailing zeros may be not taken into account */
          c = raw_str_cur[nbc.frac_part];
          raw_str_cur[nbc.frac_part] = '\0';
          strcat (str, raw_str_cur);
          /* we need to replace raw_str in its initial state in order to be
             correctly freed by mpfr_free_str */
          raw_str_cur[nbc.frac_part] = c;
        }

      /* exponent part */
      {
        char exp_fmt[8];
        char *exp_str;

        exp_fmt[0] = (spec.spec == 'A') ? 'P' : 'p';
        exp_fmt[1] = '\0';
        strcat (exp_fmt, "%+.1" MPFR_EXP_FORMAT_SPEC);
        exp_str = (char *) (*__gmp_allocate_func) (nbc.exp_part + 1);
        MPFR_ASSERTN (exp >= LONG_MIN);
        MPFR_ASSERTN (exp <= LONG_MAX);
        sprintf (exp_str, exp_fmt, exp);
        strcat (str, exp_str);
        (*__gmp_free_func) (exp_str, nbc.exp_part + 1);
      }

      mpfr_free_str (raw_str);
      MPFR_ASSERTD (nbc.total == strlen (str));
    }

  /* right justification padding */
  if ((spec.left == 0) && (spec.width > nbc.total))
    buffer_pad (buf, ' ', spec.width - nbc.total);

  buffer_cat (buf, str, nbc.total);

  /* left justification padding */
  if ((spec.left == 1) && (spec.width > nbc.total))
    buffer_pad (buf, ' ', spec.width - nbc.total);

  mpfr_free_str (str);
  return;
}

/* sprnt_fp is an internal function printing the mpfr_t P into the buffer BUF
   according to specification described by SPEC.
   note: SPEC must be a floating point conversion specification. */
static void
sprnt_fp (struct string_buffer *buf, mpfr_srcptr p, struct printf_spec spec)
{
  const char d_point[] = { MPFR_DECIMAL_POINT, '\0' };
  mp_exp_t exp;
  int base;
  int remove_trailing_zeros = 0;
  unsigned long nsd;            /* number of significant digits */
  size_t str_len;
  char *str;
  char *str_curr = NULL;

  /* char_fp details how much characters are needed in each part of a float
     print. */
  struct char_fp
  {
    int sgn;                    /* 1 if sign char is present */
    int base_prefix;            /* 0, 1 (for '0'), or 2 (for '0x' and '0X') */
    unsigned long int_part;     /* Size of integral part given by get_str */
    int point;                  /* 1 if decimal point char */
    long frac_part;             /* Size of fractional part given by get_str */
    unsigned long exp_part;     /* Size of exponent (always in base ten) */

    unsigned int total;         /* Total size */
  };

  struct char_fp nbc;

  /* special values */
  if (MPFR_UNLIKELY (MPFR_IS_NAN (p)))
    {
      sprnt_nan (buf, spec);
      return;
    }

  if (MPFR_UNLIKELY (MPFR_IS_INF (p)))
    {
      sprnt_inf (buf, spec, MPFR_SIGN (p) < 0 ? -1 : 0);
      return;
    }

  /* ordinary numbers */
  nbc.sgn = (MPFR_SIGN (p) < 0) || spec.showsign || spec.space ? 1 : 0;
  nbc.total = nbc.sgn;

  /* Replace 'g'/'G' by 'e'/'E' or 'f'/'F' following the C99 rules:
     if P > X >= -4 then the conversion is with style 'f'/'F'
     and precision P-(X+1).
     otherwise, the conversion is with style 'e'/'E'
     and precision P-1.
     where P is the threshold computed below and X is the exponent
     that would be displayed with style 'e'. */
  if ((spec.spec == 'g') || (spec.spec == 'G'))
    {
      if (MPFR_UNLIKELY (MPFR_IS_ZERO (p)))
        {
          spec.spec = (spec.spec == 'g') ? 'f' : 'F';
          spec.prec = 6;
        }
      else
        {
          mp_exp_t threshold;
          unsigned long exp;

          threshold = (spec.prec < 0) ? 6 : (spec.prec == 0) ? 1 : spec.prec;
          exp = uceil_log10_exp_p2 (p);

          if ((threshold > MPFR_GET_EXP (p) - 1) && (MPFR_GET_EXP (p) > -4))
            {
              spec.spec = (spec.spec == 'g') ? 'f' : 'F';
              spec.prec = threshold - MPFR_GET_EXP (p);
            }
          else
            {
              spec.spec = (spec.spec == 'g') ? 'e' : 'E';
              spec.prec = threshold - 1;
            }
        }
      remove_trailing_zeros = spec.alt ? 0 : 1;
    }

  /* Note: case 'g'/'G' has been changed into 'e'/'E' or 'f'/'F' above. */
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
          nbc.exp_part = uceil_log10_exp_p2 (p);
          if (nbc.exp_part < 3)
            /* the exponent always contains at least 1 digit in hexadecimal */
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
          nbc.exp_part = uceil_log10_exp_p2 (p);
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
          nbc.int_part = (MPFR_SIGN (l) < 0) || MPFR_IS_ZERO (l) ?
            1 : mpfr_get_si (l, GMP_RNDU);
          mpfr_clear (l);
        }
      nbc.point = (spec.prec == 0) && (spec.alt == 0) ? 0 : 1;
      nbc.frac_part = (spec.prec < 0) ? 6 : spec.prec;
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
          nbc.exp_part = uceil_log10_exp_p2 (p);
          if (nbc.exp_part < 4)
            /* the exponent always contains at least two digits in base 10,
               as required by the C standard */
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
          switch (*s1)
            {
            case 'a':
              *s1 = 'A';
              break;
            case 'b':
              *s1 = 'B';
              break;
            case 'c':
              *s1 = 'C';
              break;
            case 'd':
              *s1 = 'D';
              break;
            case 'e':
              *s1 = 'E';
              break;
            case 'f':
              *s1 = 'F';
              break;
            }
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

  /* build the string */
  if ((spec.left == 0) && (spec.pad == ' ') && (nbc.total < spec.width))
    /* right justification with spaces */
    buffer_pad (buf, ' ', spec.width - nbc.total);

  if (nbc.sgn)
    /* sign character */
    {
      char s[2];
      s[0] = (MPFR_SIGN (p) < 0) ? '-' : (spec.showsign) ? '+' : ' ';
      s[1] = '\0';
      buffer_cat (buf, s, 1);
      if (MPFR_SIGN (p) < 0)
        str_curr++;
    }

  if ((spec.left == 0) && (spec.pad == '0') && (nbc.total < spec.width))
    /* leading zeros in integral part */
    buffer_pad (buf, '0', spec.width - nbc.total);

  if (exp < 1 && (spec.spec == 'f' || spec.spec == 'F'))
    /* there is always a digit before the decimal point */
    buffer_cat (buf, "0", 1);
  else
    /* integral part */
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
        buffer_pad (buf, '0', -exp);

      buffer_cat (buf, str_curr, nbc.frac_part);

      if ((remove_trailing_zeros == 0) && (nbc.frac_part < spec.prec))
        /* add trailing zeros */
        buffer_pad (buf, '0', spec.prec - nbc.frac_part);
    }

  /* Note: case 'g'/'G' has been changed into 'e'/'E' or 'f'/'F' above. */
  if (nbc.exp_part)
    /* exponent part */
    {
      char exp_fmt[7];
      char *exp_str;
      switch (spec.spec)
        {
        case 'a':
        case 'A':
          buffer_cat (buf, spec.spec == 'A' ? "P" : "p", 1);
          exp = MPFR_IS_ZERO (p) ? 0 : (exp - 1) * 4;
          strcpy (exp_fmt, "%+.1" MPFR_EXP_FORMAT_SPEC);
          break;
        case 'b':
          buffer_cat (buf, "p", 1);
          exp = MPFR_IS_ZERO (p) ? 0 : exp - 1;
          strcpy (exp_fmt, "%+.1" MPFR_EXP_FORMAT_SPEC);
          break;
        case 'e':
        case 'E':
        case 'f':
        case 'F':
          buffer_cat (buf,
                      ((spec.spec == 'E') || (spec.spec == 'F')) ? "E" : "e",
                      1);
          exp = MPFR_IS_ZERO (p) ? 0 : exp - nbc.int_part;
          strcpy (exp_fmt, "%+.2" MPFR_EXP_FORMAT_SPEC);
        }

      MPFR_ASSERTN (exp - 1 >= LONG_MIN);
      MPFR_ASSERTN (exp - 1 <= LONG_MAX);
      exp_str = (char *) (*__gmp_allocate_func) (nbc.exp_part + 1);
      snprintf (exp_str, nbc.exp_part, exp_fmt, exp);

      MPFR_ASSERTD (nbc.exp_part == 1 + (unsigned long) strlen (exp_str));
      buffer_cat (buf, exp_str, nbc.exp_part - 1);
      (*__gmp_free_func) (exp_str, nbc.exp_part + 1);
    }

  if (spec.left && (spec.pad == ' ') && (nbc.total < spec.width))
    /* left justification with spaces */
    buffer_pad (buf, ' ', spec.width - nbc.total);

  mpfr_free_str (str);
}


int
mpfr_vasprintf (char **ptr, const char *fmt, va_list ap)
{
  struct string_buffer buf;
  size_t nbchar;

  /* informations on the conversion specification filled by the parser */
  struct printf_spec spec;
  /* flag raised when previous part of fmt need to be processed by
     gmp_vsnprintf */
  int gmp_fmt_flag;
  /* beginning and end of the previous unprocessed part of fmt */
  const char *start, *end;
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
          const char *f = ++fmt;
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
          int *int_ptr;         /* [FIXME] do mp{zqf} cases */
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

          l = gmp_asprintf (&s, "%" MPFR_PREC_FORMAT_SPEC, prec);
          buffer_cat (&buf, s, l);
          mpfr_free_str (s);
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
            case 'a':
            case 'A':
              sprnt_fp_a (&buf, p, spec);
              break;
            case 'f':
            case 'F':
            case 'e':
            case 'E':
            case 'g':
            case 'G':
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
  MPFR_ASSERTD (nbchar == strlen (buf.start));
  buf.start =
    (char *) (*__gmp_reallocate_func) (buf.start, buf.size, nbchar + 1);
  *ptr = buf.start;
  /* FIXME: If nbchar is larger than INT_MAX, the ISO C99 standard
     is silent, but POSIX says concerning the snprintf() function:
     "[EOVERFLOW] The value of n is greater than {INT_MAX} or the
     number of bytes needed to hold the output excluding the
     terminating null is greater than {INT_MAX}." See:
     http://www.opengroup.org/onlinepubs/009695399/functions/fprintf.html
     But it doesn't say anything concerning the other printf-like
     functions. A defect report has been submitted to austin-review-l
     (item 2532). So, for the time being, I (VL) propose that we should
     return a negative value and set the erange flag. */
  MPFR_ASSERTN (nbchar <= INT_MAX);
  return (int) nbchar;
}

#endif /* HAVE_STDARG */
