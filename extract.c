#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#include <strings.h>
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

  /* TODO: gestion de l'infini. */
  
  if (MPFR_ABSSIZE(p) < two_i){
    int j;
    y->_mp_alloc=two_i_2 ;	
    y->_mp_size=two_i_2 ;	
    PTR(y) = (*_mp_allocate_func)(two_i_2 * sizeof(mp_limb_t));
    MPN_ZERO (PTR(y), two_i_2);
    assert(y->_mp_d!=NULL);
    for(j= 0; j < MPFR_ABSSIZE(p) -  two_i_2 ; j++){
      y->_mp_d[j + two_i - MPFR_ABSSIZE(p)] = p->_mp_d[j];
    }
  } else
    {
      PTR(y) = (*_mp_allocate_func)(two_i_2  * sizeof(mp_limb_t));
      memcpy(y -> _mp_d, p->_mp_d+MPFR_ABSSIZE(p) - two_i, two_i_2 * sizeof(mp_limb_t));
      y->_mp_alloc=two_i_2 ;	
      y->_mp_size=two_i_2 ;
    }
  
  size = MPFR_ABSSIZE(y);
  for (j = 0; j < size; j++)
    {
      if (y->_mp_d[j])
         { 
	   if MPFR_ISNEG(p)
		     mpz_neg(y,y);
            return 0;
         }
    }
  y->_mp_size=0;	
  
  return 0;
}



