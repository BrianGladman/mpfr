/* mpfr_fpif -- Binary export & import of MPFR numbers
   (floating-point interchange format)

Copyright 2012-2014 Free Software Foundation, Inc.
Contributed by Olivier Demengeon.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include "mpfr-impl.h"

#if !defined (HAVE_BIG_ENDIAN) && !defined (HAVE_LITTLE_ENDIAN)
#error "Endianness is unknown. Not supported yet."
#endif

#define MPFR_KIND_ZERO 119
#define MPFR_KIND_INF 120
#define MPFR_KIND_NAN 121
#define MPFR_MAX_EMBEDDED_PRECISION 248
#define MPFR_MIN_EMBEDDED_EXPONENT -47
#define MPFR_MAX_EMBEDDED_EXPONENT 47
#define MPFR_EXTERNAL_EXPONENT 94

/* Begin: Low level helper functions */
#define COUNT_NB_BYTE(storage, size)            \
  do                                            \
    {                                           \
      (storage) >>= 8;                          \
      (size)++;                                 \
    }                                           \
  while ((storage) != 0)

#define ALLOC_RESULT(buffer, buffer_size, wanted_size)                  \
  do                                                                    \
    {                                                                   \
      if ((buffer) == NULL || *(buffer_size) < (wanted_size))           \
        {                                                               \
          (buffer) = (unsigned char *) (*__gmp_reallocate_func)         \
            ((buffer), *(buffer_size), (wanted_size));                  \
          if ((buffer) == NULL)                                         \
            {                                                           \
              *(buffer_size) = 0;                                       \
              return NULL;                                              \
            }                                                           \
        }                                                               \
      *(buffer_size) = (wanted_size);                                   \
    }                                                                   \
  while (0)

/*
 * size in byte of a MPFR number in a binary object of a variable size
 */
#define MAX_VARIABLE_STORAGE(exponent_size, precision) \
  ((size_t)(((precision) >> 3) + (exponent_size) +     \
            ((precision) > 248 ? sizeof(mpfr_prec_t) : 0) + 3))

/* copy in result[] the values in data[] with a different endianness,
   where data_size might be smaller than data_max_size, so that we only
   copy data_size bytes from the end of data[]. */
static void
#if defined (HAVE_BIG_ENDIAN)
putLittleEndianData (unsigned char * result, unsigned char * data,
                     size_t data_max_size, size_t data_size)
#elif defined (HAVE_LITTLE_ENDIAN)
putBigEndianData (unsigned char * result, unsigned char * data,
                  size_t data_max_size, size_t data_size)
#endif
{
  unsigned int j;

  MPFR_ASSERTD (data_size <= data_max_size);
  for (j = 0; j < data_size; j++)
    result[j] = data[data_max_size - j - 1];
}

/* copy in result[] the values in data[] with the same endianness */
static void
#if defined (HAVE_BIG_ENDIAN)
putBigEndianData (unsigned char * result, unsigned char * data,
                  size_t data_max_size, size_t data_size)
#elif defined (HAVE_LITTLE_ENDIAN)
putLittleEndianData (unsigned char * result, unsigned char * data,
                     size_t data_max_size, size_t data_size)
#endif
{
  MPFR_ASSERTD (data_size <= data_max_size);
  memcpy (result, data, data_size);
}

/* copy in result[] the values in data[] with a different endianness;
   the data are written at the end of the result[] buffer (if
   data_size < data_max_size, the first bytes of result[] are
   left untouched). */
static void
#if defined (HAVE_BIG_ENDIAN)
getLittleEndianData (unsigned char * result, unsigned char * data,
                     size_t data_max_size, size_t data_size)
#elif defined (HAVE_LITTLE_ENDIAN)
getBigEndianData (unsigned char * result, unsigned char * data,
                  size_t data_max_size, size_t data_size)
#endif
{
  unsigned int j;

  MPFR_ASSERTD (data_size <= data_max_size);
  for (j = 0; j < data_size; j++)
    result[data_max_size - j - 1] = data[j];
}

/* copy in result[] the values in data[] with the same endianness */
static void
#if defined (HAVE_BIG_ENDIAN)
getBigEndianData (unsigned char * result, unsigned char * data,
                  size_t data_max_size, size_t data_size)
#elif defined (HAVE_LITTLE_ENDIAN)
getLittleEndianData (unsigned char * result, unsigned char * data,
                     size_t data_max_size, size_t data_size)
#endif
{
  MPFR_ASSERTD (data_size <= data_max_size);
  memcpy (result, data, data_size);
}

/* End: Low level helper functions */

/* Internal Function */
/*
 * buffer : OUT : store the precision in binary format, can be NULL
 *               (maybe reallocated if too small)
 * buffer_size : IN/OUT : size of the buffer => size used in the buffer
 * precision : IN : precision to store
 * return pointer to a buffer storing the precision in binary format
 */
static unsigned char*
mpfr_fpif_store_precision (unsigned char * buffer, size_t * buffer_size,
                           mpfr_prec_t precision)
{
  unsigned char * result;
  size_t size_precision;

  size_precision = 0;

  if (precision > MPFR_MAX_EMBEDDED_PRECISION)
    {
      mpfr_prec_t copy_precision;

      copy_precision = precision - (MPFR_MAX_EMBEDDED_PRECISION + 1);
      COUNT_NB_BYTE(copy_precision, size_precision);
    }

  result = buffer;
  ALLOC_RESULT(result, buffer_size, size_precision + 1);

  if (precision > MPFR_MAX_EMBEDDED_PRECISION)
    {
      result[0] = size_precision - 1;
      precision -= (MPFR_MAX_EMBEDDED_PRECISION + 1);
      putLittleEndianData (result + 1, (unsigned char *) &precision,
                           sizeof(mpfr_prec_t), size_precision);
    }
  else
    result[0] = precision + 7;

  return result;
}

/*
 * fh : IN : file handler
 * return the precision stored in the binary buffer, 0 in case of error
 */
static mpfr_prec_t
mpfr_fpif_read_precision_from_file (FILE *fh)
{
  mpfr_prec_t precision;
  size_t precision_size;
  unsigned char buffer[8];
  size_t status;

  if (fh == NULL)
    return 0;

  status = fread (buffer, 1, 1, fh);
  if (status != 1)
    return 0;

  precision_size = buffer[0];
  if (precision_size >= 8)
    return precision_size - 7;

  precision_size++;

  /* Read the precision in little-endian format. */
  status = fread (buffer, precision_size, 1, fh);
  if (status != 1)
    return 0;

  while (precision_size > sizeof(mpfr_prec_t))
    {
      if (buffer[precision_size-1] != 0)
        return 0;  /* the read precision doesn't fit in a mpfr_prec_t */
      precision_size--;
    }

  if (precision_size == sizeof(mpfr_prec_t) &&
      buffer[precision_size-1] >= 0x80)
    return 0;  /* the read precision doesn't fit in a mpfr_prec_t */

  precision = 0;  /* to pad with 0's if data_size < data_max_size */

  /* On big-endian machines, the data must be copied at the end of the
     precision object in the memory; thus data_max_size (3rd argument)
     must be sizeof(mpfr_prec_t). */
  getLittleEndianData ((unsigned char *) &precision, buffer,
                       sizeof(mpfr_prec_t), precision_size);

  return precision + (MPFR_MAX_EMBEDDED_PRECISION + 1);
}

/*
 * buffer : OUT : store the kind of the MPFR number x, his sign, the size of
 *                his exponent and his exponent value in a binary format,
 *                can be NULL (maybe reallocated if too small)
 * buffer_size : IN/OUT : size of the buffer => size used in the buffer
 * x : IN : MPFR number
 * return pointer to a buffer storing the kind of the MPFR number x, his sign,
 *        the size of his exponent and his exponent value in a binary format,
 */
/* TODO
 *   exponents that use more than 16 bytes are not managed
*/
static unsigned char*
mpfr_fpif_store_exponent (unsigned char *buffer, size_t *buffer_size, mpfr_t x)
{
  unsigned char *result;
  mpfr_exp_t exponent;
  size_t exponent_size;
  char exponent_sign, sign;

  exponent = mpfr_get_exp (x);
  exponent_sign = 0;
  exponent_size = 0;

  if (mpfr_regular_p (x))
    {
      if ((exponent > MPFR_MAX_EMBEDDED_EXPONENT) ||
          (exponent < MPFR_MIN_EMBEDDED_EXPONENT))
        {
          mpfr_exp_t copy_exponent;

          if (exponent < 0)
            {
              exponent = -exponent;
              exponent_sign = 1;
            }

          exponent -= MPFR_MAX_EMBEDDED_EXPONENT;

          copy_exponent = exponent << 1;
          COUNT_NB_BYTE(copy_exponent, exponent_size);

          exponent |= exponent_sign << ((exponent_size << 3) - 1);
        }
      else
        exponent += MPFR_MAX_EMBEDDED_EXPONENT;
    }

  result = buffer;
  ALLOC_RESULT(result, buffer_size, exponent_size + 1);

  if (mpfr_regular_p (x))
    {
      if (exponent_size == 0)
        result[0] = exponent;
      else
        {
          result[0] = MPFR_EXTERNAL_EXPONENT + exponent_size;

          putLittleEndianData (result + 1, (unsigned char *) &exponent,
                               sizeof(mpfr_exp_t), exponent_size);
        }
    }
  else if (mpfr_zero_p (x))
    result[0] = MPFR_KIND_ZERO;
  else if (mpfr_inf_p (x))
    result[0] = MPFR_KIND_INF;
  else if (mpfr_nan_p (x))
    result[0] = MPFR_KIND_NAN;
  else
    {
      result = NULL;
      *buffer_size = 0;
    }

  sign = mpfr_sgn (x) > 0 ? 0 : 1;
  result[0] |= sign << 7;

  return result;
}

/*
 * x : OUT : MPFR number extracted from the binary buffer
 * fh : IN : file handler
 * return 0 if successful
 */
/* TODO
 *   exponents that use more than 16 bytes are not managed
*/
static int
mpfr_fpif_read_exponent_from_file (mpfr_t x, FILE * fh)
{
  mpfr_exp_t exponent;
  size_t exponent_size;
  int sign;
  unsigned char buffer[sizeof(mpfr_exp_t)];
  int status;
  size_t statusFile;

  status = 0;

  if (fh == NULL)
    return 1;

  statusFile = fread (buffer, 1, 1, fh);
  if (statusFile != 1)
    return 1;

  sign = (buffer[0] & 0x80) ? -1 : 0;
  exponent = buffer[0] & 0x7F;
  exponent_size = 1;

  if ((exponent > MPFR_EXTERNAL_EXPONENT) && (exponent < MPFR_KIND_ZERO))
    {
      mpfr_exp_t exponent_sign;

      exponent_size = exponent - MPFR_EXTERNAL_EXPONENT;

      if (exponent_size > sizeof(mpfr_exp_t))
        return 1;

      statusFile = fread (buffer, exponent_size, 1, fh);
      if (statusFile != 1)
        return 1;

      exponent = 0;
      getLittleEndianData ((unsigned char *) &exponent, buffer,
                           sizeof(mpfr_exp_t), exponent_size);

      exponent_sign = exponent & (1 << ((exponent_size << 3) - 1));

      exponent &= ~exponent_sign;
      exponent += MPFR_MAX_EMBEDDED_EXPONENT;

      if (exponent_sign != 0)
        exponent = -exponent;

      mpfr_setsign (x, x, sign, MPFR_RNDN);
      status = mpfr_set_exp (x, exponent);

      exponent_size++;
    }
  else if (exponent == MPFR_KIND_ZERO)
    mpfr_set_zero (x, sign);
  else if (exponent == MPFR_KIND_INF)
    mpfr_set_inf (x, sign);
  else if (exponent == MPFR_KIND_NAN)
    mpfr_set_nan (x);
  else if (exponent < 95)
    status = mpfr_set_exp (x, exponent - MPFR_MAX_EMBEDDED_EXPONENT);
  else
    return 1;

  return status;
}

/*
 * buffer : OUT : store the limb of the MPFR number x in a binary format,
 *                can be NULL (maybe reallocated if too small)
 * buffer_size : IN/OUT : size of the buffer => size used in the buffer
 * x : IN : MPFR number
 * return pointer to a buffer storing the limb of the MPFR number x in a binary
 *        format
 */
static unsigned char*
mpfr_fpif_store_limbs (unsigned char *buffer, size_t *buffer_size, mpfr_t x)
{
  unsigned char *result;
  mpfr_prec_t precision;
  size_t nb_byte;
  size_t nb_limb, mp_bytes_per_limb;
  size_t nb_partial_byte;
  size_t i, j;

  precision = mpfr_get_prec (x);
  nb_byte = (precision + 7) >> 3;
  mp_bytes_per_limb = mp_bits_per_limb >> 3;
  nb_partial_byte = nb_byte % mp_bytes_per_limb;
  nb_limb = (nb_byte + mp_bytes_per_limb - 1) / mp_bytes_per_limb;

  result = buffer;
  ALLOC_RESULT(result, buffer_size, nb_byte);

  putBigEndianData (result, (unsigned char*) MPFR_MANT(x),
                    sizeof(mp_limb_t), nb_partial_byte);
  for (i = nb_partial_byte, j = (nb_partial_byte == 0) ? 0 : 1; j < nb_limb;
       i += mp_bytes_per_limb, j++)
    putLittleEndianData (result + i, (unsigned char*) (MPFR_MANT(x) + j),
                         sizeof(mp_limb_t), sizeof(mp_limb_t));

  return result;
}

/*
 * x : OUT : MPFR number extracted from the binary buffer, should have the same
 *           precision than the number in the binary format
 * buffer : IN : limb of the MPFR number x in a binary format,
 * buffer_size : IN/OUT : size of the buffer => size used in the buffer
 * return 0 if successful
 */
static int
mpfr_fpif_read_limbs (mpfr_t x, unsigned char *buffer, size_t *buffer_size)
{
  mpfr_prec_t precision;
  size_t nb_byte;
  size_t mp_bytes_per_limb;
  size_t nb_partial_byte;
  size_t i, j;

  precision = mpfr_get_prec (x);
  nb_byte = (precision + 7) >> 3;
  mp_bytes_per_limb = mp_bits_per_limb >> 3;
  nb_partial_byte = nb_byte % mp_bytes_per_limb;

  if ((buffer == NULL) || (*buffer_size < nb_byte))
    {
      *buffer_size = 0;
      return 1;
    }
  *buffer_size = nb_byte;

  if (nb_partial_byte > 0)
    {
      memset (MPFR_MANT(x), 0, sizeof(mp_limb_t));
      getBigEndianData ((unsigned char*) MPFR_MANT(x), buffer,
                        sizeof(mp_limb_t), nb_partial_byte);
    }
  for (i = nb_partial_byte, j = (nb_partial_byte == 0) ? 0 : 1; i < nb_byte;
       i += mp_bytes_per_limb, j++)
    getLittleEndianData ((unsigned char*) (MPFR_MANT(x) + j), buffer + i,
                         sizeof(mp_limb_t), sizeof(mp_limb_t));

  return 0;
}

/* External Function */
/*
 * fh : IN : file hander
 * x : IN : MPFR number to put in the file
 * return 0 if successful
 */
int
mpfr_fpif_export (FILE *fh, mpfr_t x)
{
  int status;
  unsigned char *buf;
  unsigned char *bufResult;
  size_t used_size, buf_size;

  if (fh == NULL)
    return -1;

  buf_size = MAX_VARIABLE_STORAGE(sizeof(mpfr_exp_t), mpfr_get_prec (x));
  buf = (unsigned char*) (*__gmp_allocate_func) (buf_size);
  if (buf == NULL)
    return -1;

  used_size = buf_size;
  buf = mpfr_fpif_store_precision (buf, &used_size, mpfr_get_prec (x));
  used_size > buf_size ? buf_size = used_size : 0;
  status = fwrite (buf, used_size, 1, fh);
  if (status != 1)
    {
      (*__gmp_free_func) (buf, buf_size);
      return -1;
    }
  used_size = buf_size;
  bufResult = mpfr_fpif_store_exponent (buf, &used_size, x);
  if (bufResult == NULL)
    {
      (*__gmp_free_func) (buf, buf_size);
      return -1;
    }
  buf = bufResult;
  used_size > buf_size ? buf_size = used_size : 0;
  status = fwrite (buf, used_size, 1, fh);
  if (status != 1)
    {
      (*__gmp_free_func) (buf, buf_size);
      return -1;
    }

  if (mpfr_regular_p (x))
    {
      used_size = buf_size;
      buf = mpfr_fpif_store_limbs (buf, &used_size, x);
      used_size > buf_size ? buf_size = used_size : 0;
      status = fwrite (buf, used_size, 1, fh);
      if (status != 1)
        {
          (*__gmp_free_func) (buf, buf_size);
          return -1;
        }
    }

  (*__gmp_free_func) (buf, buf_size);
  return 0;
}

/*
 * x : IN/OUT : MPFR number extracted from the file, his precision is reset to
 *              be able to hold the number
 * fh : IN : file hander
 * Return 0 if the import was successful.
 */
int
mpfr_fpif_import (mpfr_t x, FILE *fh)
{
  int status;
  mpfr_prec_t precision;
  unsigned char *buffer;
  size_t used_size;

  precision = mpfr_fpif_read_precision_from_file (fh);
  if (precision == 0) /* precision = 0 means an error */
    return -1;
  if (precision > MPFR_PREC_MAX)
    return -1;
  MPFR_ASSERTN (MPFR_PREC_MIN <= 8);
  if (precision < MPFR_PREC_MIN)
    precision = MPFR_PREC_MIN;
  mpfr_set_prec (x, precision);

  status = mpfr_fpif_read_exponent_from_file (x, fh);
  if (status != 0)
    return -1;

  if (mpfr_regular_p (x))
    {
      used_size = (precision + 7) >> 3; /* ceil(precision/8) */
      buffer = (unsigned char*) (*__gmp_allocate_func) (used_size);
      if (buffer == NULL)
        {
          return -1;
        }
      status = fread (buffer, used_size, 1, fh);
      if (status != 1)
        {
          (*__gmp_free_func) (buffer, used_size);
          return -1;
        }
      status = mpfr_fpif_read_limbs (x, buffer, &used_size);
      (*__gmp_free_func) (buffer, used_size);
      if (status != 0)
        return -1;
    }

  return 0;
}
