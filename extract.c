#include <stdio.h>
#include <math.h>
#include <stdlib.h> 

#include "gmp.h"
#include "gmp-impl.h"
#include "mpfr.h"
#include "assert.h"


int
#if __STDC__
mpfr_extract(mpz_ptr y, mpfr_srcptr p, int  i)
#else
mpfr_extract(y,p,i)
mpz_ptr y; 
mpfr_srcptr p;
int  i;
#endif
{
  int two_i=1 << i;
  int two_i_2 = i ? two_i / 2 : 1;
  int size;
  int j;
  if (ABSSIZE(p) < two_i){
    int j;
    y->_mp_alloc=two_i_2 ;	
    y->_mp_size=two_i_2 ;	
    y->_mp_d = calloc(two_i_2,sizeof(  mp_limb_t));
    assert(y->_mp_d!=NULL);
    /* attention : initialiser a 0 si on utilise malloc */
    for(j= 0; j < ABSSIZE(p) -  two_i_2 ; j++){
      y->_mp_d[j + two_i - ABSSIZE(p)] = p->_mp_d[j];
    }
  } else
    {
      y->_mp_d = malloc(two_i_2  * sizeof(mp_limb_t));
      memcpy(y -> _mp_d, p->_mp_d+ABSSIZE(p) - two_i, two_i_2 * sizeof(mp_limb_t));
      y->_mp_alloc=two_i_2 ;	
      y->_mp_size=two_i_2 ;
    }
  
  size = ABSSIZE(y);
  for (j = 0; j < size; j++)
    {
      if (y->_mp_d[j])
         { 
	   if ISNEG(p)
		     mpz_neg(y,y);
            return 0;
         }
    }
  y->_mp_size=0;	
  
  return 0;
}



