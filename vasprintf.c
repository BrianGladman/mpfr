/* mpfr_vasprintf -- main function for the printf functions family
   plus helper macros & functions.

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

#define TRUE -1
#define FALSE 0

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
#define MPFR_NAN_STRING_LENGTH 3
#define MPFR_INF_STRING_LC "inf"
#define MPFR_INF_STRING_UC "INF"
#define MPFR_INF_STRING_LENGTH 3

/* We assume that a single conversion specifier produces at most 4095 chars
   (Rationale for International Standard -Programming Languages- C
   Revision 5.10 April-2003, 7.19.6.1 p.152).
   MAX_CHAR_BY_SPEC must be less than INT_MAX to be compatible with
   mpfr_vasprintf() return type. */
#define MAX_CHAR_BY_SPEC 4096

static const char num_to_text[16] = "0123456789abcdef";

/* some macro and functions for parsing format string */
#define READ_INT(ap, format, specinfo, field, label_out)                \
  do {                                                                  \
    while (*(format))                                                   \
      {                                                                 \
        switch (*(format))                                              \
          {                                                             \
          case '0':                                                     \
          case '1':                                                     \
          case '2':                                                     \
          case '3':                                                     \
          case '4':                                                     \
          case '5':                                                     \
          case '6':                                                     \
          case '7':                                                     \
          case '8':                                                     \
          case '9':                                                     \
            MPFR_ASSERTN (specinfo.field < MAX_CHAR_BY_SPEC / 10);      \
            specinfo.field *= 10;                                       \
            MPFR_ASSERTN (specinfo.field < MAX_CHAR_BY_SPEC-*(format)+'0'); \
            specinfo.field += *(format) - '0';                          \
            ++(format);                                                 \
            break;                                                      \
          case '*':                                                     \
            specinfo.field = va_arg ((ap), int);                        \
            MPFR_ASSERTN (specinfo.field < MAX_CHAR_BY_SPEC);           \
            ++(format);                                                 \
          default:                                                      \
            goto label_out;                                             \
          }                                                             \
      }                                                                 \
  } while (0)

/* arg_t contains all the types described by the 'type' field of the
   format string */
enum arg_t
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
   This structure is adapted from the GNU libc one. */
struct printf_spec
{
  unsigned int alt:1;           /* # flag */
  unsigned int space:1;         /* Space flag */
  unsigned int left:1;          /* - flag */
  unsigned int showsign:1;      /* + flag */
  unsigned int group:1;         /* ' gnu flag (not for gmp/mpfr types) */

  int width;                    /* Width */
  int prec;                     /* Precision */

  enum arg_t arg_type;          /* Type of argument */
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
#else
#define CASE_INTMAX_ARG(specinfo, ap)
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
        * ((len + 1 > MAX_CHAR_BY_SPEC) ?
           len + 1 : MAX_CHAR_BY_SPEC);
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

  if (n == 0)
    return;
  padding = (char *) (*__gmp_allocate_func) (n + 1);
  memset (padding, c, n);
  padding[n] = '\0';
  buffer_cat (b, padding, n);
  (*__gmp_free_func) (padding, n + 1);
}

/* Print NaN with padding. */
static int
sprnt_nan (struct string_buffer *buf, const struct printf_spec spec)
{
  /* right justification padding */
  if ((spec.left == 0) && (spec.width > MPFR_NAN_STRING_LENGTH))
    buffer_pad (buf, ' ', spec.width - MPFR_NAN_STRING_LENGTH);

  switch (spec.spec)
    {
    case 'A':
    case 'E':
    case 'F':
    case 'G':
    case 'X':
      buffer_cat (buf, MPFR_NAN_STRING_UC, MPFR_NAN_STRING_LENGTH);
      break;

    default:
      buffer_cat (buf, MPFR_NAN_STRING_LC, MPFR_NAN_STRING_LENGTH);
    }

  /* left justification padding */
  if ((spec.left == 1) && (spec.width > MPFR_NAN_STRING_LENGTH))
    buffer_pad (buf, ' ', spec.width - MPFR_NAN_STRING_LENGTH);

  return (spec.width > MPFR_NAN_STRING_LENGTH) ? spec.width
    : MPFR_NAN_STRING_LENGTH;
}

/* Print Infinities with padding
   NEG = -1 for '-inf', 0 for 'inf' */
static int
sprnt_inf (struct string_buffer *buf, const struct printf_spec spec,
           const int neg)
{
  const int length = MPFR_INF_STRING_LENGTH - neg;

  /* right justification padding */
  if ((spec.left == 0) && (spec.width > length))
    buffer_pad (buf, ' ', spec.width - length);

  switch (spec.spec)
    {
    case 'A':
    case 'E':
    case 'F':
    case 'G':
    case 'X':
      if (neg < 0)
        buffer_cat (buf, "-", 1);
      buffer_cat (buf, MPFR_INF_STRING_UC, MPFR_INF_STRING_LENGTH);
      break;
    default:
      if (neg < 0)
        buffer_cat (buf, "-", 1);
      buffer_cat (buf, MPFR_INF_STRING_LC, MPFR_INF_STRING_LENGTH);
    }

  /* left justification padding */
  if ((spec.left == 1) && (spec.width > length))
    buffer_pad (buf, ' ', spec.width - length);

  return (spec.width > length) ? spec.width : length;
}

/* let gmp_xprintf process the part it can understand */
static int
sprntf_gmp (struct string_buffer *b, const char *fmt, va_list ap)
{
  int length;
  char *s;

  length = gmp_vasprintf (&s, fmt, ap);
  if (length > 0)
    buffer_cat (b, s, length);

  mpfr_free_str (s);
  return length;
}

/* internal function printing the mpfr_t P into the buffer BUF according to
   specification described by SPEC.
   Note: SPEC must be a integer conversion specification.
   Return the number of written characters.
   Return -1 if the built string has more than MAX_CHAR_BY_SPEC
   characters. */
static int
sprnt_int (struct string_buffer *buf, mpfr_srcptr p, struct printf_spec spec)
{
  char format[10];  /* buffer for the format string corresponding to spec
                       (see below) */
  int f;
  char *s;
  int length;
  mpz_t z;

  if (MPFR_UNLIKELY (MPFR_IS_NAN (p)))
    return sprnt_nan (buf, spec);

  if (MPFR_UNLIKELY (MPFR_IS_INF (p)))
    return sprnt_inf (buf, spec, MPFR_SIGN (p) < 0 ? -1 : 0);

  mpz_init (z);
  mpfr_get_z (z, p, spec.rnd_mode);

  /* FORMAT contains at most 9 characters plus the terminating '\0'
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

  length = gmp_asprintf (&s, format, spec.width, spec.prec, z);

  if ((0 < length) && (length < MAX_CHAR_BY_SPEC))
    buffer_cat (buf, s, length);
  else
    length = -1;

  mpfr_free_str (s);
  mpz_clear (z);
  return length;
}

/* struct for easy string clearing */
struct string_list
{
  char *string;
  struct string_list *next; /* NULL in last node */
};

/* initialisation */
static void
init_string_list (struct string_list *sl)
{
  sl->string = NULL;
  sl->next = NULL;
}

/* clear all strings in the list */
static void
clear_string_list (struct string_list *sl)
{
  struct string_list *n;

  while (sl)
    {
      if (sl->string)
        mpfr_free_str (sl->string);
      n = sl->next;
      (*__gmp_free_func) (sl, sizeof(struct string_list));
      sl = n;
    }
}

/* add a string in the list */
static char *
register_string (struct string_list *sl, char *new_string)
{
  /* look for the last node */
  while (sl->next)
    sl = sl->next;

  sl->next = (struct string_list*)
    (*__gmp_allocate_func) (sizeof (struct string_list));

  sl = sl->next;
  sl->next = NULL;
  return sl->string = new_string;
}

/* padding type: where are the padding characters */
enum pad_t
  {
    LEFT,          /* spaces in left hand side for right justification */
    LEADING_ZEROS, /* padding with '0' characters in integral part */
    RIGHT          /* spaces in right hand side for left justification */
  };

/* number_parts details how much characters are needed in each part of a float
   print.  */
struct number_parts
{
  enum pad_t pad_type;    /* Padding type */
  int pad_size;           /* Number of padding characters */

  char sign;              /* Sign character */

  char *prefix_ptr;       /* Pointer to prefix part */
  int prefix_size;        /* number of characters in *prefix_ptr */

  char *ip_ptr;           /* Pointer to integral part characters*/
  int ip_size;            /* Number of digits in *ip_ptr */
  int ip_trailing_zeros;  /* Number of additional null digits in integral
                             part */

  char point;             /* Decimal point character */

  int fp_leading_zeros;   /* Number of additional leading zeros in fractional
                             part */
  char *fp_ptr;           /* Pointer to fractional part characters */
  int fp_size;            /* Number of digits in *fp_ptr */
  int fp_trailing_zeros;  /* Number of additional trailing zeros in fractional
                             part */

  char *exp_ptr;          /* Pointer to exponent part */
  int exp_size;           /* Number of characters in *exp_ptr */

  struct string_list *sl; /* List of string buffers in use: we need such a
                             mechanism because fp_ptr may point into the same
                             string as ip_ptr */
};

/* Determine the different parts of the string representation of the regular
   number p when spec.spec is 'a', 'A', or 'b'.
   return -1 if some field > MAX_CHAR_BY_SPEC */
static int
regular_ab (struct number_parts *np, mpfr_srcptr p,
            const struct printf_spec spec)
{
  int uppercase;
  int base;
  char *str;
  mp_exp_t exp;

  uppercase = spec.spec == 'A';

  /* sign */
  if (MPFR_IS_NEG (p))
    np->sign = '-';
  else if (spec.showsign || spec.space)
    np->sign = spec.showsign ? '+' : ' ';

  if (spec.spec == 'a' || spec.spec == 'A')
    /* prefix part */
    {
      np->prefix_size = 2;
      str = (char *) (*__gmp_allocate_func) (1 + np->prefix_size);
      str[0] = '0';
      str[1] = uppercase ? 'X' : 'x';
      str[2] = '\0';
      np->prefix_ptr = register_string (np->sl, str);
    }

  /* integral part */
  np->ip_size = 1;
  base = (spec.spec == 'b') ? 2 : 16;

  if (spec.spec == 'b' || spec.prec != 0)
    /* In order to avoid ambiguity in rounding to even case, we will ever
       output at least one fractional digit in binary mode */
    {
      int nsd;

      /* Number of significant digits:
         - if no given precision, let mpfr_get_str determine it;
         - if a zero precision is specified and if we are in binary mode, then
         ask for two binary digits, one before decimal point, and one after;
         - if a non-zero precision is specified, then one digit before decimal
         point plus SPEC.PREC after it. */
      nsd = spec.prec < 0 ? 0
        : (spec.prec == 0 && spec.spec == 'b') ? 2 : spec.prec + np->ip_size;
      str = mpfr_get_str (0, &exp, base, nsd, p, spec.rnd_mode);
      register_string (np->sl, str);
      np->ip_ptr = MPFR_IS_NEG (p) ? ++str : str;  /* skip sign if any */

      if (base == 16)
        /* EXP is the exponent for radix sixteen with decimal point BEFORE the
           first digit, we want the exponent for radix two and the decimal
           point AFTER the first digit. */
        exp = (exp - 1) * 4;
      else
        /* EXP is the exponent for decimal point BEFORE the first digit, we
           want the exponent for decimal point AFTER the first digit. */
        --exp;
    }
  else
    /* One hexadecimal digit is sufficient but mpfr_get_str returns at least
       two digits when the base is a power of two.
       So, in order to avoid double rounding, we will build our own string. */
    {
      mp_limb_t *pm = MPFR_MANT (p);
      mp_size_t ps;
      int digit, shift;
      int rnd_away;

      /* rnd_away:
         1 if round away from zero,
         0 if round to zero,
         -1 if not decided yet. */
      rnd_away =
        spec.rnd_mode == GMP_RNDD ? MPFR_IS_NEG (p) :
        spec.rnd_mode == GMP_RNDU ? MPFR_IS_POS (p) :
        spec.rnd_mode == GMP_RNDZ ? 0 : -1;

      /* exponent for radix-2 with the decimal point after the first
         hexadecimal digit */
      exp = MPFR_GET_EXP (p) - 4; /* FIXME: possible overflow */

      /* Determine the radix-16 digit by grouping the 4 first digits. Even
         if MPFR_PREC (p) < 4, we can read 4 bits in its first limb */
      shift = BITS_PER_MP_LIMB - 4;
      ps = (MPFR_PREC (p) - 1) / BITS_PER_MP_LIMB;
      digit = pm[ps]>>shift;

      if (MPFR_PREC (p) > 4)
        /* round taking into account bits outside the first 4 ones */
        {
          if (rnd_away == -1)
            /* Round to nearest mode: we have to decide in that particular
               case if we have to round away from zero or not */
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
                  digit >>= 1;
                  ++exp;  /* no possible overflow because
                             exp == EXP(p)-3 */
                }
            }
        }

      MPFR_ASSERTD ((0 <= digit) && (digit <= 15));
      np->ip_size = 1;
      str = (char *)(*__gmp_allocate_func) (1 + np->ip_size);
      str[0] = num_to_text [digit];
      str[1] = '\0';

      np->ip_ptr = register_string (np->sl, str);
    }

  if (uppercase)
    /* All digits in upper case */
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

  if (spec.spec == 'b' || spec.prec != 0)
    /* compute the number of digits in fractional part */
    {
      char *ptr;
      size_t str_len;

      /* the sign has been skipped, skip also the first digit */
      ++str;
      str_len = strlen (str);
      ptr = str + str_len - 1; /* points to the end of str */

      if (spec.prec < 0)
        /* remove trailing zeros, if any */
        {
          while ((*ptr == '0') && str_len)
            {
              --ptr;
              --str_len;
            }
        }

      if (str_len > INT_MAX || str_len > MAX_CHAR_BY_SPEC)
        /* too much digits in fractional part */
        return -1;

      if (str_len)
        /* there are some non-zero digits in fractional part */
        {
          np->fp_ptr = str;
          np->fp_size = (int) str_len;
          if ((int) str_len < spec.prec)
            np->fp_trailing_zeros = spec.prec - str_len;
        }
    }

  /* decimal point */
  if (np->fp_size || spec.alt)
    np->point = MPFR_DECIMAL_POINT;

  /* the exponent part contains the character 'p', or 'P' plus the sign
     character plus at least one digit and only as many more digits as
     necessary to represent the exponent.
     We assume that |EXP| < 10^INT_MAX. */
  np->exp_size = 3;
  {
    mp_exp_unsigned_t x;

    x = SAFE_ABS (mp_exp_unsigned_t, exp);
    while (x > 9)
      {
        np->exp_size++;
        x /= 10;
      }
  }
  str = (char *) (*__gmp_allocate_func) (1 + np->exp_size);
  np->exp_ptr = register_string (np->sl, str);
  {
    char exp_fmt[8];  /* contains at most 7 characters like in "p%+.1i",
                         or "P%+.2li" */

    exp_fmt[0] = uppercase ? 'P' : 'p';
    exp_fmt[1] = '\0';
    strcat (exp_fmt, "%+.1" MPFR_EXP_FORMAT_SPEC);

    if (sprintf (str, exp_fmt, exp) < 0)
      return -1;
  }

  return 0;
}

/* Determine the different parts of the string representation of the regular
   number p when spec.spec is 'e', 'E', 'g', or 'G'.

   return -1 if some field > MAX_CHAR_BY_SPEC */
static int
regular_eg (struct number_parts *np, mpfr_srcptr p,
            const struct printf_spec spec)
{
  int uppercase;
  char *str;
  mp_exp_t exp;

  uppercase = spec.spec == 'E' || spec.spec == 'G';

  /* sign */
  if (MPFR_IS_NEG (p))
    np->sign = '-';
  else if (spec.showsign || spec.space)
    np->sign = spec.showsign ? '+' : ' ';

  /* integral part */
  np->ip_size = 1;
  {
    int nsd;

    /* Number of significant digits:
       - if no given precision, then let mpfr_get_str determine it,
       - if a precision is specified, then one digit before decimal point
       plus SPEC.PREC after it.
       We use the fact here that mpfr_get_exp allows us to ask for only one
       significant digit when the base is not a power of 2. */
    nsd = (spec.prec < 0) ? 0 : spec.prec + np->ip_size;
    str = mpfr_get_str (0, &exp, 10, nsd, p, spec.rnd_mode);
  }
  register_string (np->sl, str);
  np->ip_ptr = MPFR_IS_NEG (p) ? ++str : str;  /* skip sign if any */

  if (spec.prec != 0)
    /* compute the number of digits in fractional part */
    {
      char *ptr;
      size_t str_len;

      /* the sign has been skipped, skip also the first digit */
      ++str;
      str_len = strlen (str);
      ptr = str + str_len - 1; /* points to the end of str */

      if ((spec.prec < 0)
          && ((spec.spec != 'g' && spec.spec != 'G') || !spec.alt))
        /* remove trailing zeros, if any */
        {
          while ((*ptr == '0') && str_len)
            {
              --ptr;
              --str_len;
            }
        }

      if (str_len > INT_MAX || str_len > MAX_CHAR_BY_SPEC)
        /* too much digits in fractional part */
        return -1;

      if (str_len)
        /* there are some non-zero digits in fractional part */
        {
          np->fp_ptr = str;
          np->fp_size = (int) str_len;
          if ((int) str_len < spec.prec)
            np->fp_trailing_zeros = spec.prec - str_len;
        }
    }

  /* decimal point */
  if (np->fp_size || spec.alt)
    np->point = MPFR_DECIMAL_POINT;

  /* EXP is the exponent for decimal point BEFORE the first digit, we want
     the exponent for decimal point AFTER the first digit. */
  exp--;

  /* the exponent part contains the character 'e', or 'E' plus the sign
     character plus at least two digits and only as many more digits as
     necessary to represent the exponent.
     We assume that |EXP| < 10^INT_MAX. */
  np->exp_size = 3;
  {
    mp_exp_unsigned_t x;

    x = SAFE_ABS (mp_exp_unsigned_t, exp);
    while (x > 9)
      {
        np->exp_size++;
        x /= 10;
      }
  }
  if (np->exp_size < 4)
    np->exp_size = 4;

  str = (char *) (*__gmp_allocate_func) (1 + np->exp_size);
  np->exp_ptr = register_string (np->sl, str);

  {
    char exp_fmt[8];  /* e.g. "e%+.2i", or "E%+.2li" */

    exp_fmt[0] = uppercase ? 'E' : 'e';
    exp_fmt[1] = '\0';
    strcat (exp_fmt, "%+.2" MPFR_EXP_FORMAT_SPEC);

    if (sprintf (str, exp_fmt, exp) < 0)
      return -1;
  }

  return 0;
}

/* Determine the different parts of the string representation of the regular
   number p when spec.spec is 'f', 'F', 'g', or 'G'.

   return -1 if some field > MAX_CHAR_BY_SPEC */
static int
regular_fg (struct number_parts *np, mpfr_srcptr p,
            const struct printf_spec spec)
{
  int keep_trailing_zeros;
  char * str;
  mpfr_t x;

  keep_trailing_zeros = (spec.spec == 'g' || spec.spec == 'G') && spec.alt;

  /* sign */
  if (MPFR_IS_NEG (p))
    np->sign = '-';
  else if (spec.showsign || spec.space)
    np->sign = spec.showsign ? '+' : ' ';

  /* Determine the position of the most significant decimal digit. */
  {
    /* Let p = m*10^e with 1 <= m < 10 and p = n*2^d with 0.5 <= d < 1.
       We need at most 1+log2(floor(d/3)+1) bits of precision in order to
       represent the exact value of e+1 if p >= 1, or |e| if p < 1. */
    mp_prec_t m;
    mp_prec_t n;

    m = (mp_prec_t) SAFE_ABS (mp_exp_unsigned_t, MPFR_GET_EXP (p));
    m /= 3;
    m++;
    n = 1;
    while (m)
      {
        m >>= 1;
        n++;
      }

    if (n < MPFR_PREC (p))
      mpfr_init2 (x, MPFR_PREC (p));
    else
      mpfr_init2 (x, n);
  }
  mpfr_abs (x, p, GMP_RNDD); /* With our choice of precision,
                                x == |p| exactly. */

  if (MPFR_GET_EXP (x) < 0)
    /* 0 < p < 1 */
    {
      if (spec.prec != 0)
        {
          /* integral part, one digit: "0" */
          np->ip_size = 1;
          str = (char *) (*__gmp_allocate_func) (1 + np->ip_size);
          str[0] = '0';
          str[1] = '\0';
          np->ip_ptr = register_string (np->sl, str);

          /* decimal point */
          np->point = MPFR_DECIMAL_POINT;

          mpfr_log10 (x, x, GMP_RNDD);
          mpfr_floor (x, x);
          mpfr_abs (x, x, GMP_RNDD);
          /* We have rounded away from zero so that x == |e| (with p = m*10^e,
             see above). Now, x is the number of digits in the fractional
             part. */

          if ((spec.prec > 0) && (mpfr_cmp_si (x, spec.prec + 1) > 0))
            /* p is too small for the given precision,
               output "0.0_00" or "0.0_01" depending on rnd_mode */
            {
              int rnd_away;
              /* rnd_away:
                 1 if round away from zero,
                 0 if round to zero,
                 -1 if not decided yet. */
              rnd_away =
                spec.rnd_mode == GMP_RNDD ? MPFR_IS_NEG (p) :
                spec.rnd_mode == GMP_RNDU ? MPFR_IS_POS (p) :
                spec.rnd_mode == GMP_RNDZ ? 0 : -1;

              if (rnd_away < 0)
                {
                  int c;

                  /* let compare p with x = sign(p)*5*10^(-spec.prec-1) */
                  mpfr_set_si (x, -spec.prec - 1, GMP_RNDN);
                  mpfr_exp10 (x, x, GMP_RNDN);
                  mpfr_mul_ui (x, x, 5, GMP_RNDN);
                  MPFR_SET_SAME_SIGN (x, p);
                  c = mpfr_cmp (p, x);

                  if (MPFR_LIKELY (c != 0))
                    /* round to zero if p and c have same sign */
                    rnd_away = MPFR_IS_POS (p) ^ (c < 0);
                  else
                    /* tie case, round to even (i.e. zero) */
                    rnd_away = 0;
                }

              if (rnd_away)
                /* the last output digit is '1' */
                {
                  np->fp_leading_zeros = spec.prec - 1;

                  np->fp_size = 1;
                  str = (char *) (*__gmp_allocate_func) (1 + np->fp_size);
                  str[0] = '1';
                  str[1] = '\0';
                  np->fp_ptr = register_string (np->sl, str);
                }
              else
                /* only spec.prec zeros in fractional part */
                np->fp_leading_zeros = spec.prec;
            }
          else
            /* some significant digits can be output in the fractional part */
            {
              int n;
              int nsd;
              mp_exp_t exp;
              size_t str_len;

              n = mpfr_get_ui (x, GMP_RNDZ);

              np->fp_leading_zeros = n - 1;
              /* WARNING: nsd may equal 1, we use here the fact that
                 mpfr_get_str can return one digit with base ten
                 (undocumented feature, see comments in get_str.c) */
              nsd = spec.prec < 0 ? 0 : spec.prec - n + 1;
              str = mpfr_get_str (NULL, &exp, 10, nsd, p, spec.rnd_mode);
              register_string (np->sl, str);
              np->fp_ptr = MPFR_IS_NEG (p) ? ++str : str; /* skip sign */
              np->point = MPFR_DECIMAL_POINT;

              str_len = strlen (str); /* the sign has been skipped */
              str += str_len - 1; /* points to the end of str */

              if (!keep_trailing_zeros)
                /* remove trailing zeros, if any */
                {
                  while ((*str == '0') && str_len)
                    {
                      --str;
                      --str_len;
                    }
                }

              if ((str_len > INT_MAX)
                  || (str_len > MAX_CHAR_BY_SPEC))
                /* too much digits in fractional part */
                {
                  mpfr_clear (x);
                  return -1;
                }
              np->fp_size = str_len;

              MPFR_ASSERTD (n == 1 - exp);
            }
        }
      else
        /* spec.prec == 0 */
        {
          int rnd_away;
          /* rnd_away:
             1 if round away from zero,
             0 if round to zero,
             -1 if not decided yet. */
          rnd_away =
            spec.rnd_mode == GMP_RNDD ? MPFR_IS_NEG (p) :
            spec.rnd_mode == GMP_RNDU ? MPFR_IS_POS (p) :
            spec.rnd_mode == GMP_RNDZ ? 0 : -1;

          if (rnd_away < 0)
            {
              int c;

              /* let compare p with x = sign(p)*5*10^(-spec.prec-1) */
              mpfr_set_si (x, -spec.prec - 1, GMP_RNDN);
              mpfr_exp10 (x, x, GMP_RNDN);
              mpfr_mul_ui (x, x, 5, GMP_RNDN);
              MPFR_SET_SAME_SIGN (x, p);
              c = mpfr_cmp (p, x);

              if (MPFR_LIKELY (c != 0))
                /* round to zero if p and c have same sign */
                rnd_away = MPFR_IS_POS (p) ^ (c < 0);
              else
                /* tie case, round to even (i.e. zero) */
                rnd_away = 0;
            }

          /* integral part only */
          np->ip_size = 1;
          str = (char *) (*__gmp_allocate_func) (1 + np->ip_size);
          str[0] = rnd_away ? '1' : '0';
          str[1] = '\0';
          np->ip_ptr = register_string (np->sl, str);

          if (spec.alt)
            np->point = MPFR_DECIMAL_POINT;
        }
    }
  else
    /* 1 <= p */
    {
      mp_exp_t exp;
      int nsd;  /* Number of significant digits */

      mpfr_log10 (x, x, GMP_RNDZ);
      mpfr_floor (x, x);
      mpfr_add_ui (x, x, 1, GMP_RNDZ);
      /* We have rounded towards zero so that x == e + 1 (with p = m*10^e,
         see above). x is now the number of digits in the integral part. */

      if (mpfr_cmp_si (x, MAX_CHAR_BY_SPEC) > 0)
        /* P is too large to print all its integral part digits */
        {
          mpfr_clear (x);
          return -1;
        }

      np->ip_size = mpfr_get_si (x, GMP_RNDN);

      nsd = spec.prec < 0 ? 0 : spec.prec + np->ip_size;
      str = mpfr_get_str (NULL, &exp, 10, nsd, p, spec.rnd_mode);
      register_string (np->sl, str);
      np->ip_ptr = MPFR_IS_NEG (p) ? ++str : str; /* skip sign */

      if (nsd == 0)
        /* compute how much non-zero digits in integral and fractional
           parts */
        {
          size_t str_len;
          str_len = strlen (str); /* note: the sign has been skipped */

          if (np->ip_size > str_len)
            /* mpfr_get_str doesn't give the trailing zeros when p is a
               multiple of 10 (p integer, so no fractional part) */
            {
              np->ip_trailing_zeros = np->ip_size - str_len;
              np->ip_size = str_len;
            }
          else
            /* str may contain some digits which are in fractional part */
            {
              char *ptr;

              ptr = str + str_len - 1; /* points to the end of str */
              str_len -= np->ip_size;  /* number of digits in fractional
                                          part */

              if (!keep_trailing_zeros)
                /* remove trailing zeros, if any */
                {
                  while ((*ptr == '0') && str_len)
                    {
                      --ptr;
                      --str_len;
                    }
                }

              if ((str_len > INT_MAX)
                  || (str_len > MAX_CHAR_BY_SPEC))
                /* too much digits in fractional part */
                {
                  mpfr_clear (x);
                  return -1;
                }

              if (str_len)
                /* some digits in fractional part */
                {
                  np->point = MPFR_DECIMAL_POINT;
                  np->fp_ptr = str + np->ip_size;
                  np->fp_size = str_len;
                }
              else if (spec.alt)
                np->point = MPFR_DECIMAL_POINT;
            }
        }
      else
        /* spec.prec digits in fractional part */
        {
          MPFR_ASSERTD (np->ip_size == exp);

          if (spec.prec)
            {
              np->point = MPFR_DECIMAL_POINT;
              np->fp_ptr = str + np->ip_size;
              np->fp_size = spec.prec;
            }
          else if (spec.alt)
            np->point = MPFR_DECIMAL_POINT;
        }
    }

  mpfr_clear (x);
  return 0;
}

/* partition_number determines the different parts of the string
   representation of the number p according to the given specification.
   partition_number initializes the given structure np, so all previous
   information in that variable is lost.
   return the total number of characters to be written.
   return -1 if an error occured, in that case np's fields are in an undefined
   state but all string buffers have been freed. */
static int
partition_number (struct number_parts *np, mpfr_srcptr p,
                  struct printf_spec spec)
{
  char *str;
  int total;
  int uppercase;

  /* WARNING: left justification means right space padding */
  np->pad_type = spec.left ? RIGHT : spec.pad == '0' ? LEADING_ZEROS : LEFT;
  np->pad_size = 0;
  np->sign = '\0';
  np->prefix_ptr =NULL;
  np->prefix_size = 0;
  np->ip_ptr = NULL;
  np->ip_size = 0;
  np->ip_trailing_zeros = 0;
  np->point = '\0';
  np->fp_leading_zeros = 0;
  np->fp_ptr = NULL;
  np->fp_size = 0;
  np->fp_trailing_zeros = 0;
  np->exp_ptr = NULL;
  np->exp_size = 0;
  np->sl = (struct string_list *)
    (*__gmp_allocate_func) (sizeof (struct string_list));
  init_string_list (np->sl);

  uppercase = spec.spec == 'A' || spec.spec == 'E' || spec.spec == 'F'
    || spec.spec == 'G';

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (p)))
    {
      if (MPFR_IS_NAN (p))
        {
          if (np->pad_type == LEADING_ZEROS)
            /* don't want "0000nan", change to right justification padding
               with left spaces instead */
            np->pad_type = LEFT;

          if (uppercase)
            {
              np->ip_size = MPFR_NAN_STRING_LENGTH;
              str = (char *) (*__gmp_allocate_func) (1 + np->ip_size);
              strcpy (str, MPFR_NAN_STRING_UC);
              np->ip_ptr = register_string (np->sl, str);
            }
          else
            {
              np->ip_size = MPFR_NAN_STRING_LENGTH;
              str = (char *) (*__gmp_allocate_func) (1 + np->ip_size);
              strcpy (str, MPFR_NAN_STRING_LC);
              np->ip_ptr = register_string (np->sl, str);
            }
        }
      else if (MPFR_IS_INF (p))
        {
          if (np->pad_type == LEADING_ZEROS)
            /* don't want "0000inf", change to right justification padding
               with left spaces instead */
            np->pad_type = LEFT;

          if (uppercase)
            {
              if (MPFR_IS_NEG (p))
                np->sign = '-';
              np->ip_size = MPFR_INF_STRING_LENGTH;
              str = (char *) (*__gmp_allocate_func) (1 + np->ip_size);
              strcpy (str, MPFR_INF_STRING_UC);
              np->ip_ptr = register_string (np->sl, str);
            }
          else
            {
              if (MPFR_IS_NEG (p))
                np->sign = '-';
              np->ip_size = MPFR_INF_STRING_LENGTH;
              str = (char *) (*__gmp_allocate_func) (1 + np->ip_size);
              strcpy (str, MPFR_INF_STRING_LC);
              np->ip_ptr = register_string (np->sl, str);
            }
        }
      else
        /* p == 0 */
        {
          if (MPFR_IS_NEG (p))
            /* signed zero */
            np->sign = '-';
          else if (spec.showsign || spec.space)
            np->sign = spec.showsign ? '+' : ' ';

          if (spec.spec == 'a' || spec.spec == 'A')
            /* prefix part */
            {
              np->prefix_size = 2;
              str = (char *) (*__gmp_allocate_func) (1 + np->prefix_size);
              str[0] = '0';
              str[1] = uppercase ? 'X' : 'x';
              str[2] = '\0';
              np->prefix_ptr = register_string (np->sl, str);
            }

          /* integral part */
          np->ip_size = 1;
          str = (char *) (*__gmp_allocate_func) (1 + np->ip_size);
          str[0] = '0';
          str[1] = '\0';
          np->ip_ptr = register_string (np->sl, str);

          if (spec.prec > 0)
            /* fractional part */
            {
              np->point = MPFR_DECIMAL_POINT;
              np->fp_trailing_zeros = spec.prec;
            }
          else if (spec.alt)
            np->point = MPFR_DECIMAL_POINT;

          if (spec.spec == 'a' || spec.spec == 'A' || spec.spec == 'b'
              || spec.spec == 'e' || spec.spec == 'E')
            /* exponent part */
            {
              np->exp_size = (spec.spec == 'e' || spec.spec == 'E') ? 4 : 3;
              str = (char *) (*__gmp_allocate_func) (1 + np->exp_size);
              if (spec.spec == 'e' || spec.spec == 'E')
                strcpy (str, uppercase ? "E+00" : "e+00");
              else
                strcpy (str, uppercase ? "P+0" : "p+0");
              np->exp_ptr = register_string (np->sl, str);
            }
        }
    }
  else
    /* regular p, p != 0 */
    {
      if (spec.spec == 'a' || spec.spec == 'A' || spec.spec == 'b')
        {
          if (regular_ab (np, p, spec))
            goto error;
        }
      else if (spec.spec == 'f' || spec.spec == 'F')
        {
          if (regular_fg (np, p, spec))
            goto error;
        }
      else if (spec.spec == 'e' || spec.spec == 'E')
        {
          if (regular_eg (np, p, spec))
            goto error;
        }
      else
        /* %g case */
        {
          /* Replace 'g'/'G' by 'e'/'E' or 'f'/'F' following the C99 rules:
             if P > X >= -4 then the conversion is with style 'f'/'F' and
             precision P-(X+1).
             otherwise, the conversion is with style 'e'/'E' and
             precision P-1.
             where P is the threshold computed below and X is the exponent
             that would be displayed with style 'e'. */
          int threshold;
          mpfr_t x;

          threshold = (spec.prec < 0) ? 6 : (spec.prec == 0) ? 1 : spec.prec;

          mpfr_init2 (x, 53);
          mpfr_log10 (x, p, GMP_RNDD);
          mpfr_rint (x, x, GMP_RNDD);
          if (mpfr_cmp_si (x, threshold) < 0 && mpfr_cmp_si (x, -4) >= 0)
            {
              spec.prec = threshold - 1 - mpfr_get_si (x, GMP_RNDD);
              mpfr_clear (x);

              if (regular_fg (np, p, spec))
                goto error;
            }
          else
            {
              spec.prec = threshold - 1;
              mpfr_clear (x);

              if (regular_eg (np, p, spec))
                goto error;
            }
        }
    }

  /* compute the number of characters to be written verifying it is not too
     much */
  total = np->sign ? 1 : 0;
  total += np->prefix_size;
  total += np->ip_size;
  if (total < 0 || total > MAX_CHAR_BY_SPEC)
    goto error;
  total += np->ip_trailing_zeros;
  if (total < 0 || total > MAX_CHAR_BY_SPEC)
    goto error;
  total += np->point ? 1 : 0;
  total += np->fp_leading_zeros;
  if (total < 0 || total > MAX_CHAR_BY_SPEC)
    goto error;
  total += np->fp_size;
  if (total < 0 || total > MAX_CHAR_BY_SPEC)
    goto error;
  total += np->fp_trailing_zeros;
  if (total < 0 || total > MAX_CHAR_BY_SPEC)
    goto error;
  total += np->exp_size;
  if (total < 0 || total > MAX_CHAR_BY_SPEC)
    goto error;

  if (spec.width > total)
    /* pad with spaces or zeros depending on np->pad_type */
    {
      np->pad_size = spec.width - total;
      total += np->pad_size; /* here total == spec.width,
                                so 0 < total < MAX_CHAR_BY_SPEC */
    }

  return total;

 error:
  clear_string_list (np->sl);
  np->ip_ptr = NULL;
  np->fp_ptr = NULL;
  np->exp_ptr = NULL;
  return -1;
}

/* sprnt_fp prints a mpfr_t according to spec.spec specification.
   return the size of the string (not counting the terminating '\0')
   return -1 if the built string is too long (i.e. has more than
   MAX_CHAR_BY_SPEC characters). */
static int
sprnt_fp (struct string_buffer *buf, mpfr_srcptr p,
          const struct printf_spec spec)
{
  int length;
  struct number_parts np;

  length = partition_number (&np, p, spec);
  if (length < 0)
    return -1;

  /* right justification padding with left spaces */
  if (np.pad_type == LEFT && np.pad_size)
    buffer_pad (buf, ' ', np.pad_size);

  /* sign character (may be '-', '+', or ' ') */
  if (np.sign)
    buffer_pad (buf, np.sign, 1);

  /* prefix part */
  if (np.prefix_ptr)
    buffer_cat (buf, np.prefix_ptr, np.prefix_size);

  /* right justification  padding with leading zeros */
  if (np.pad_type == LEADING_ZEROS && np.pad_size)
    buffer_pad (buf, '0', np.pad_size);

  /* integral part (may also be "nan" or "inf") */
  if (np.ip_ptr)
    buffer_cat (buf, np.ip_ptr, np.ip_size);

  /* trailing zeros in integral part */
  if (np.ip_trailing_zeros)
    buffer_pad (buf, '0', np.ip_trailing_zeros);

  /* decimal point */
  if (np.point)
    buffer_pad (buf, np.point, 1);

  /* leading zeros in fractional part */
  if (np.fp_leading_zeros)
    buffer_pad (buf, '0', np.fp_leading_zeros);

  /* significant digits in fractional part */
  if (np.fp_ptr)
    buffer_cat (buf, np.fp_ptr, np.fp_size);

  /* trailing zeros in fractional part */
  if (np.fp_trailing_zeros)
    buffer_pad (buf, '0', np.fp_trailing_zeros);

  /* exponent part */
  if (np.exp_ptr)
    buffer_cat (buf, np.exp_ptr, np.exp_size);

  /* left justication padding with right spaces */
  if (np.pad_type == RIGHT && np.pad_size)
    buffer_pad (buf, ' ', np.pad_size);

  clear_string_list (np.sl);
  return length;
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

  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_SAVE_EXPO_MARK (expo);

  nbchar = 0;
  buffer_init (&buf, (MAX_CHAR_BY_SPEC + 1) * sizeof (char));
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

      READ_INT (ap, fmt, spec, width, width_analysis);
    width_analysis:
      if (spec.width < 0)
        {
          spec.left = 1;
          spec.width = -spec.width;
          MPFR_ASSERTN (spec.width < MAX_CHAR_BY_SPEC);
        }
      if (*fmt == '.')
        {
          const char *f = ++fmt;
          READ_INT (ap, fmt, spec, prec, prec_analysis);
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
        /* end of the format string */
        break;
      else if (spec.spec == 'n')
        /* put the number of characters written so far in the location pointed
           by the next va_list argument */
        {
          int *int_ptr;         /* [TODO] do mp{zqf} cases */
          int_ptr = va_arg (ap, int *);
          FLUSH (gmp_fmt_flag, start, end, ap2, &buf);
          va_end (ap2);
          va_copy (ap2, ap);
          start = fmt;
          *int_ptr = (int) (buf.curr - buf.start);
        }
      else if (spec.arg_type == MPFR_PREC_ARG)
        /* output mp_prec_t variable */
        {
          char *s;
          int length;
          mpfr_prec_t prec;
          prec = va_arg (ap, mpfr_prec_t);

          FLUSH (gmp_fmt_flag, start, end, ap2, &buf);
          va_end (ap2);
          va_copy (ap2, ap);
          start = fmt;

          length = gmp_asprintf (&s, "%" MPFR_PREC_FORMAT_SPEC, prec);
          if ((long)buf.size + length <= INT_MAX)
            {
              buffer_cat (&buf, s, length);
              mpfr_free_str (s);
            }
          else
            {
              mpfr_free_str (s);
              goto error;
            }
        }
      else if (spec.arg_type == MPFR_ARG)
        /* output a mpfr_t variable */
        {
          int length = 0;
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
              length = sprnt_int (&buf, p, spec);
              break;
            case 'a':
            case 'A':
            case 'b':
            case 'e':
            case 'E':
            case 'f':
            case 'F':
            case 'g':
            case 'G':
              length = sprnt_fp (&buf, p, spec);
            }

          if ((length < 0) || ((long)buf.size + length > INT_MAX))
            goto error;
        }
      else
        /* gmp_printf specification, step forward in the va_list */
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

  /* If nbchar is larger than INT_MAX, the ISO C99 standard is silent, but
     POSIX says concerning the snprintf() function:
     "[EOVERFLOW] The value of n is greater than {INT_MAX} or the
     number of bytes needed to hold the output excluding the
     terminating null is greater than {INT_MAX}." See:
     http://www.opengroup.org/onlinepubs/009695399/functions/fprintf.html
     But it doesn't say anything concerning the other printf-like functions.
     A defect report has been submitted to austin-review-l (item 2532).
     So, for the time being, we return a negative value and set the erange
     flag. */
  if (nbchar <= INT_MAX)
    {
      MPFR_SAVE_EXPO_FREE (expo);
      return (int) nbchar;
    }

 error:
  MPFR_SAVE_EXPO_FREE (expo);
  mpfr_set_erangeflag ();
  *ptr = NULL;
  (*__gmp_free_func) (buf.start, buf.size);
  return -1;
}

#endif /* HAVE_STDARG */
