#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <mpfr.h>

int main(int argc, char *argv[])
{
  char *filenameCompressed = NULL;
  int delFile = 1;
  signed int status;
  FILE *fh;
  mpfr_t x[9];
  mpfr_t y;
  size_t limb_size, exponent_size;
  signed int i;
  unsigned char * buf;
  size_t buf_size, max_buf_size, buf_pos;
  mpfr_prec_t precision;

  if (argc == 2)
  {
    filenameCompressed = argv[1];
    delFile = 1;
  }

  mpfr_init2(x[0], 130);
  mpfr_init2(x[8], 130);
  mpfr_inits2(2048, x[1], x[2], x[3], x[4], x[5], x[6], x[7], NULL);
  mpfr_set_d(x[0], 45.2564215, MPFR_RNDN);
  mpfr_set_d(x[1], 45.2564215, MPFR_RNDN);
  mpfr_set_d(x[2], 45.2564215, MPFR_RNDN);
  mpfr_set_exp (x[2], -48000);
  mpfr_set_inf(x[3], -1);
  mpfr_set_zero(x[4], 0);
  mpfr_set_nan(x[5]);
  mpfr_set_d(x[6], 104348, MPFR_RNDN);
  mpfr_set_d(x[7], 33215, MPFR_RNDN);
  mpfr_div(x[8], x[6], x[7], MPFR_RNDN);
  mpfr_div(x[6], x[6], x[7], MPFR_RNDN);
    
  if (filenameCompressed == NULL)
  {
    filenameCompressed = "stream_data.test";
    fh = fopen(filenameCompressed, "w");  
    if (fh == NULL)
    {
      fprintf(stderr, "Failled to open for writing %s, exiting...\n", 
        filenameCompressed);
      return -1;
    }
  
    for(i=0; i<9; i++)
    {
      status = mpfr_fpif_export_binary(fh, x[i]);
      if (status != 0)
      {
        fclose(fh);
        if (delFile == 0)
        {
          unlink(filenameCompressed);
        }
        fprintf(stderr, "Failled to export the %d number, exiting...\n", i);
        return -1;
      }
    }
  
    fclose(fh);
  }

  fh = fopen(filenameCompressed, "r");  
  if (fh == NULL)
  {
    if (delFile == 0)
    {
      unlink(filenameCompressed);
    }
    fprintf(stderr, "Failled to open for reading %s, exiting...\n", 
      filenameCompressed);
    return -1;
  }

  for(i=0; i<9; i++)
  {
    mpfr_init2(y, 2);
    mpfr_fpif_import_binary(fh, y);
    if (mpfr_cmp(x[i], y) != 0)
    {
      fprintf(stderr, "mpfr_cmp failed on the %d number, exiting...\n", i);
      if (delFile == 0)
      {
        unlink(filenameCompressed);
      }
      return -1;
    }
    mpfr_clear(y);
  }
  fclose(fh);
  if (delFile == 0)
  {
    unlink(filenameCompressed);
  }

  return 0;
}

