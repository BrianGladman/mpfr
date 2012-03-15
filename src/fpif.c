/* mpfr_fpif -- Binary export / import of mpfr number.

Copyright 2002, 2003, 2004, 2006, 2007, 2008, 2009, 2010, 2011, 2012 Free Software Foundation, Inc.
Contributed by the AriC and Caramel projects, INRIA.

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

#include <mpfr.h>

/* FIXME:
    - <endian.h> is not standard.
    - Use the specific alloc/free functions.
    - Use () around arguments in macros.
    - Use the GNU coding style.
*/

#include <endian.h>
#ifndef __BYTE_ORDER
#error "__BYTE_ORDER is not defined"
#endif
#if __BYTE_ORDER == __BIG_ENDIAN
#define IS_BIG_ENDIAN
#elif __BYTE_ORDER == __LITTLE_ENDIAN
#else
#error "__BYTE_ORDER not managed"
#endif

#define MPFR_KIND_ZERO 119
#define MPFR_KIND_INF 120
#define MPFR_KIND_NAN 121
#define MPFR_MAX_EMBEDDED_PRECISION 248
#define MPFR_MIN_EMBEDDED_EXPONENT -47
#define MPFR_MAX_EMBEDDED_EXPONENT 47
#define MPFR_EXTERNAL_EXPONENT 94

/* Begin: Low level helper functions */
#define COUNT_NB_BYTE(storage, size)\
  do\
  {\
    storage >>= 8;\
    size++;\
  }\
  while(storage != 0);
  
#define ALLOC_RESULT(buffer, buffer_size, wanted_size)\
  if ((buffer == NULL) || (*buffer_size < (wanted_size)))\
  {\
    buffer = (unsigned char *) realloc (buffer, wanted_size);\
    if (buffer == NULL)\
    {\
      *buffer_size = 0;\
      return NULL;\
    }\
  }\
  *buffer_size = wanted_size;

/*
 * size in byte of a mpfr number in a binary object of a variable size
 */
#define MAX_VARIABLE_STORAGE(exponent_size, precision) \
  ((size_t)((precision >> 3) + exponent_size + \
    (precision > 248 ? sizeof(mpfr_prec_t) : 0) + 3))

#ifdef IS_BIG_ENDIAN
static void putLittleEndianData (unsigned char * result, unsigned char * data, 
  size_t data_max_size, size_t data_size, unsigned int nb_data)
#else
/* Intel */
static void putBigEndianData (unsigned char * result, unsigned char * data, 
  size_t data_max_size, size_t data_size, unsigned int nb_data)
#endif
{
  signed int i,j;

  for(i = 0; i < nb_data; i++)
  {
    for(j = 0; j < data_size; j++)
    {
      result[(i * data_size) + j] =
        data[(i * data_max_size) + data_max_size - j - 1];
    }
  }
}

#ifdef IS_BIG_ENDIAN
static void putBigEndianData (unsigned char * result, unsigned char * data, 
  size_t data_max_size, size_t data_size, unsigned int nb_data)
#else
// Intel
static void putLittleEndianData (unsigned char * result, unsigned char * data, 
  size_t data_max_size, size_t data_size, unsigned int nb_data)
#endif
{
  int i;
  for(i = 0; i < nb_data; i++)
  {
    memcpy (&result[i * data_size], &data[i * data_max_size], data_size);
  }
}

#ifdef IS_BIG_ENDIAN
static void getLittleEndianData (unsigned char * result, unsigned char * data, 
  size_t data_max_size, size_t data_size, unsigned int nb_data)
#else
// Intel
static void getBigEndianData (unsigned char * result, unsigned char * data, 
  size_t data_max_size, size_t data_size, unsigned int nb_data)
#endif
{
  signed int i,j;

  for(i = 0; i < nb_data; i++)
  {
    for(j = 0; j < data_size; j++)
    {
      result[(i * data_max_size) + data_max_size - j - 1] =
        data[(i * data_size) + j];
    }
  }
}

#ifdef IS_BIG_ENDIAN
static void getBigEndianData (unsigned char * result, unsigned char * data, 
  size_t data_max_size, size_t data_size, unsigned int nb_data)
#else
/* Intel */
static void getLittleEndianData (unsigned char * result, unsigned char * data, 
  size_t data_max_size, size_t data_size, unsigned int nb_data)
#endif
{
  signed int i;

  for(i = 0; i < nb_data; i++)
  {
    memcpy (&result[i * data_size], &data[i * data_max_size], data_size);
  }
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
static unsigned char * mpfr_fpif_store_precision
  (unsigned char * buffer, size_t * buffer_size, mpfr_prec_t precision)
{
  unsigned char * result;
  size_t size_precision;
  
  size_precision = 0;
  
  if (precision > MPFR_MAX_EMBEDDED_PRECISION)
  {
    mpfr_prec_t copy_precision;

    copy_precision = precision - (MPFR_MAX_EMBEDDED_PRECISION + 1);
    COUNT_NB_BYTE(copy_precision, size_precision)
  }
  
  result = buffer;
  ALLOC_RESULT(result, buffer_size, size_precision + 1)
  
  if (precision > MPFR_MAX_EMBEDDED_PRECISION)
  {
    result[0] = size_precision - 1;
    precision -= (MPFR_MAX_EMBEDDED_PRECISION + 1);
    putLittleEndianData (&result[1], (unsigned char *) &precision, 
      sizeof(precision), size_precision, 1);
  }
  else
  {
    result[0] = precision + 7;
  }
  
  return result;
}

/*
 * fh : IN : file handler
 * return the precision stored in the binary buffer
 */
static mpfr_prec_t mpfr_fpif_read_precision_from_file(FILE *fh)
{
  mpfr_prec_t precision;
  size_t precision_size;
  unsigned char buffer[sizeof(mpfr_prec_t)];
  size_t status;
  
  if (fh == NULL)
  {
    return 0;
  }

  status = fread(&buffer[0], 1, 1, fh);  
  if (status != 1)
  {
    return 0;
  }
  precision_size = buffer[0];
    
  if (precision_size < 8)
  {
    status = fread(&buffer[0], precision_size + 1, 1, fh);
    if (status != 1)
    {
      return 0;
    }

    precision = 0;
    getLittleEndianData ((unsigned char *)&precision, &buffer[0], 
      sizeof(precision), precision_size + 1, 1);
    precision += (MPFR_MAX_EMBEDDED_PRECISION + 1);
  }
  else
  {
    precision = precision_size - 7;
  }
    
  return precision;
}

/*
 * buffer : OUT : store the kind of the mpfr number x, his sign, the size of 
 *                his exponent and his exponent value in a binary format, 
 *                can be NULL (maybe reallocated if too small)
 * buffer_size : IN/OUT : size of the buffer => size used in the buffer
 * x : IN : mpfr number
 * return pointer to a buffer storing the kind of the mpfr number x, his sign, 
 *        the size of his exponent and his exponent value in a binary format, 
 */
/* Todo
 *   exposant that use more than 16 bytes are not managed
*/
static unsigned char * mpfr_fpif_store_exponent
  (unsigned char * buffer, size_t * buffer_size, mpfr_t x)
{
  unsigned char * result;
  mpfr_exp_t exponent;
  size_t exponent_size;
  signed char exponent_sign, sign;

  exponent = mpfr_get_exp (x);
  exponent_sign = 0;
  exponent_size = 0;
  
  if (mpfr_regular_p (x) != 0)
  {
    if ((exponent > MPFR_MAX_EMBEDDED_EXPONENT) || (exponent < MPFR_MIN_EMBEDDED_EXPONENT))
    {
      mpfr_exp_t copy_exponent;
      
      if (exponent < 0)
      {
        exponent = -exponent;
        exponent_sign = 1;
      }

      exponent -= MPFR_MAX_EMBEDDED_EXPONENT;
            
      copy_exponent = exponent << 1;
      COUNT_NB_BYTE(copy_exponent, exponent_size)

      exponent |= exponent_sign << ((exponent_size << 3) - 1);
    }
    else
    {
      exponent += MPFR_MAX_EMBEDDED_EXPONENT;
    }
  }
  
  result = buffer;
  ALLOC_RESULT(result, buffer_size, exponent_size + 1)
  
  if (mpfr_regular_p (x) != 0)
  {
    if (exponent_size == 0)
    {
      result[0] = exponent;
    }
    else
    {
      result[0] = MPFR_EXTERNAL_EXPONENT + exponent_size;
      
      putLittleEndianData (&result[1], (unsigned char *) &exponent, 
        sizeof (exponent), exponent_size, 1);
    }
  }
  else if (mpfr_zero_p (x) != 0)
  {
    result[0] = MPFR_KIND_ZERO;
  }
  else if (mpfr_inf_p (x) != 0)
  {
    result[0] = MPFR_KIND_INF;
  }
  else if (mpfr_nan_p (x) != 0)
  {
    result[0] = MPFR_KIND_NAN;
  }
  else
  {
    result = NULL;
    *buffer_size = 0;
  }
  
  sign = mpfr_sgn(x) > 0 ? 0 : 1;
  result[0] |= sign << 7;
  
  return result;
}

/*
 * x : OUT : mpfr number extracted from the binary buffer
 * fh : IN : file handler
 * return 0 if successfull
 */
/* Todo
 *   exposant that use more than 16 bytes are not managed
*/
static int mpfr_fpif_read_exponent_from_file(mpfr_t x, FILE * fh)
{
  mpfr_exp_t exponent;
  size_t exponent_size;
  signed char sign;
  unsigned char buffer[sizeof(mpfr_exp_t)];
  signed int status;
  size_t statusFile;
    
  status = 0;

  if (fh == NULL)
  {
    return 1;
  }

  statusFile = fread(&buffer[0], 1, 1, fh);  
  if (statusFile != 1)
  {
    return 1;
  }  
  
  sign = -(buffer[0] & 0x80);
  exponent = buffer[0] & 0x7F;
  exponent_size = 1;
  
  if ((exponent > MPFR_EXTERNAL_EXPONENT) && (exponent < MPFR_KIND_ZERO))
  {
    mpfr_exp_t exponent_sign;
    
    exponent_size = exponent - MPFR_EXTERNAL_EXPONENT;
    
    if (exponent_size > sizeof(exponent))
    {
      return 1;
    }


    statusFile = fread(&buffer[0], exponent_size, 1, fh);  
    if (statusFile != 1)
    {
      return 1;
    } 

    exponent = 0;
    getLittleEndianData ((unsigned char *) &exponent, &buffer[0],
      sizeof (exponent), exponent_size, 1);
    
    exponent_sign = exponent & (1 << ((exponent_size << 3) - 1));
    
    exponent &= ~exponent_sign;
    exponent += MPFR_MAX_EMBEDDED_EXPONENT;
    
    if (exponent_sign != 0)
    {
      exponent = -exponent;
    }
    
    mpfr_setsign(x, x, sign, MPFR_RNDN);
    status = mpfr_set_exp (x, exponent);

    exponent_size++;
  }
  else if (exponent == MPFR_KIND_ZERO)
  {
    mpfr_set_zero(x, sign);
  }
  else if (exponent == MPFR_KIND_INF)
  {
    mpfr_set_inf(x, sign);
  }
  else if (exponent == MPFR_KIND_NAN)
  {
    mpfr_set_nan(x);
  }
  else if (exponent < 95)
  {
    status = mpfr_set_exp(x, exponent - MPFR_MAX_EMBEDDED_EXPONENT);
  }
  else
  {
    return 1;
  }
  
  return status;
}

/*
 * buffer : OUT : store the limb of the mpfr number x in a binary format, 
 *                can be NULL (maybe reallocated if too small)
 * buffer_size : IN/OUT : size of the buffer => size used in the buffer
 * x : IN : mpfr number
 * return pointer to a buffer storing the limb of the mpfr number x in a binary 
 *        format
 */
static unsigned char * mpfr_fpif_store_limbs
  (unsigned char * buffer, size_t * buffer_size, mpfr_t x)
{
  unsigned char * result;
  mpfr_prec_t precision;
  size_t nb_byte;
  size_t nb_limb, mp_bytes_per_limb;
  size_t nb_partial_byte;
  size_t i, j;  
  
  precision = mpfr_get_prec(x);
  nb_byte = (precision + 7) >> 3;
  mp_bytes_per_limb = mp_bits_per_limb >> 3;
  nb_partial_byte = nb_byte % mp_bytes_per_limb;
  nb_limb = (nb_byte + mp_bytes_per_limb - 1) / mp_bytes_per_limb;
  
  result = buffer;
  ALLOC_RESULT(result, buffer_size, nb_byte)

  putBigEndianData (result, (unsigned char*)x->_mpfr_d, 
      sizeof(*x->_mpfr_d), nb_partial_byte, 1);
  for(i=nb_partial_byte,  j=(nb_partial_byte == 0) ? 0 : 1; j<nb_limb; 
    i+=mp_bytes_per_limb, j++)
  {
    putLittleEndianData (&result[i], (unsigned char*)(&x->_mpfr_d[j]), 
      sizeof(*x->_mpfr_d), sizeof(*x->_mpfr_d), 1);
  }

  return result;
}

/*
 * x : OUT : mpfr number extracted from the binary buffer, should have the same 
 *           precision than the number in the binary format
 * buffer : IN : limb of the mpfr number x in a binary format, 
 * buffer_size : IN/OUT : size of the buffer => size used in the buffer
 * return 0 if successfull
 */
static int mpfr_fpif_read_limbs
  (mpfr_t x, unsigned char * buffer, size_t * buffer_size)
{
  mpfr_prec_t precision;
  size_t nb_byte;
  size_t nb_limb, mp_bytes_per_limb;
  size_t nb_partial_byte;
  size_t i, j;

  precision = mpfr_get_prec(x);
  nb_byte = (precision + 7) >> 3;
  mp_bytes_per_limb = mp_bits_per_limb >> 3;
  nb_partial_byte = nb_byte % mp_bytes_per_limb;
  nb_limb = (nb_byte + mp_bytes_per_limb - 1) / mp_bytes_per_limb;
  
  if ((buffer == NULL) || (*buffer_size < nb_byte))
  {
    *buffer_size = 0;
    return 1;
  }
  *buffer_size = nb_byte;
  
  if (nb_partial_byte > 0)
  {
    memset(x->_mpfr_d, 0, sizeof(*x->_mpfr_d));
    getBigEndianData ((unsigned char*)x->_mpfr_d, buffer, 
      sizeof(*x->_mpfr_d), nb_partial_byte, 1);
  }
  for(i=nb_partial_byte, j=(nb_partial_byte == 0) ? 0 : 1; i<nb_byte; 
    i+=mp_bytes_per_limb, j++)
  {
    getLittleEndianData ((unsigned char*)(&x->_mpfr_d[j]), &buffer[i], 
      sizeof(*x->_mpfr_d), sizeof(*x->_mpfr_d), 1);
  }
  
  return 0;
}

/* External Function */
/*
 * fh : IN : file hander
 * x : IN : mpfr number to put in the file
 * return 0 if successfull
 */
int mpfr_fpif_export_binary(struct FILE *fh, mpfr_t x)
{
  signed int status;
  unsigned char * buf;
  size_t used_size, buf_size;
  
  buf_size = MAX_VARIABLE_STORAGE(sizeof(mpfr_exp_t), mpfr_get_prec(x));
  buf = malloc (buf_size);
  if (buf == NULL)
  {
    return -1;
  }
  
  used_size = buf_size;
  buf = mpfr_fpif_store_precision(buf, &used_size, mpfr_get_prec(x));
  used_size > buf_size ? buf_size = used_size : 0;
  status = fwrite (buf, used_size, 1, fh);
  if (status != 1)
  {
    free(buf);
    return -1;
  }  
  used_size = buf_size;
  buf = mpfr_fpif_store_exponent (buf, &used_size, x);
  used_size > buf_size ? buf_size = used_size : 0;
  status = fwrite (buf, used_size, 1, fh);
  if (status != 1)
  {
    free(buf);
    return -1;
  }

  if (mpfr_regular_p(x) != 0)
  {
    used_size = buf_size;
    buf = mpfr_fpif_store_limbs (buf, &used_size, x);
    used_size > buf_size ? buf_size = used_size : 0;
    status = fwrite (buf, used_size, 1, fh);
    if (status != 1)
    {
      free(buf);
      return -1;
    }
  }
  
  free(buf);
  return 0;
}

/*
 * fh : IN : file hander
 * x : IN/OUT : mpfr number extracted from the file, his precision is reset to
 *              be able to hold the number,
 */
int mpfr_fpif_import_binary(struct FILE *fh, mpfr_t x)
{
  int status;
  mpfr_prec_t precision;
  unsigned char * buffer;
  size_t used_size;
  
  precision = mpfr_fpif_read_precision_from_file(fh);
  if (precision == 0)
  {
    return -1;
  }
  mpfr_set_prec(x, precision);

  status = mpfr_fpif_read_exponent_from_file(x, fh);
  if (status != 0)
  {
    return -1;
  }

  if (mpfr_regular_p(x) != 0)
  {
    used_size = (precision/8)+((precision%8) == 0 ? 0 : 1);
    buffer = malloc(used_size);
    status = fread(buffer, used_size, 1, fh);
    if (status != 1)
    {
      free(buffer);
      return -1;
    }
    status = mpfr_fpif_read_limbs(x, &buffer[0], &used_size);
    free(buffer);
    if (status != 0)
    {
      return -1;
    }
  }
  
  return 0;
}

