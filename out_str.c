#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"

size_t 
#if __STDC__
mpfr_out_str (FILE *stream, int base, size_t n_digits, mpfr_srcptr op,
	      unsigned char rnd_mode)
#else
mpfr_out_str (stream, base, n_digits, op, rnd_mode)
     FILE *stream; 
     int base;
     size_t n_digits; 
     mpfr_srcptr op; 
     unsigned char rnd_mode; 
#endif
{
  char *s,*s0; size_t l; mp_exp_t e;

  if (FLAG_NAN(op)) { fprintf(stream, "NaN"); return 3; }
  if (!NOTZERO(op)) { fprintf(stream, "0"); return 1; }

  s = mpfr_get_str(NULL, &e, base, n_digits, op, rnd_mode);
  s0 = s;
  /* for op=3.1416 we have s = "31416" and e = 1 */
  
  l = strlen(s)+1;
  if (*s == '-') fputc(*s++, stream);

  fputc(*s++, stream); e--; /* writes leading digit */
  fputc('.', stream);       /* decimal point */
  fputs(s, stream);         /* rest of mantissa */
  if (e) {
    fputc((base>10) ? '@' : 'e', stream); l++;
    sprintf(s, "%ld", e);
    l += strlen(s);
    fprintf(stream, "%s", s);
  }

  (*_mp_free_func)(s0, l); 
  return l;
}
