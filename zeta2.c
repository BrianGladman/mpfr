#include <stdio.h>
#include <math.h>
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
#include "mpfr.h"

main()
{mpfr_t s,s2,x,y,u,b,v,nn,z,z2,p,result,res_p; int i,n,succes; size_t t;
mpfr_init(result);
mpfr_init(res_p);
mpfr_set_default_prec(67);
mpfr_init(x);
mpfr_init(y);
mpfr_init(s);
mpfr_init(s2);
mpfr_init(u);
mpfr_init(b);
mpfr_init(v);
mpfr_init(nn);
mpfr_init(z);
mpfr_init(z2);
mpfr_init(p);
mpfr_set_ui(u,1,GMP_RNDN);
mpfr_set_ui(s,0,GMP_RNDN);
/*s=Somme des 1/i^2 (i=100...2)*/
n=100;
for (i=n; i>1; i--)
{
mpfr_set_ui(x,i*i,GMP_RNDN);
mpfr_div(y,u,x,GMP_RNDN);
mpfr_add(s,s,y,GMP_RNDN);
};
/*mpfr_print_raw(s);printf("\n");
t=mpfr_out_str(stdout,10,0,s,GMP_RNDN);printf("\n");*/
/*Application d'Euler-Maclaurin, jusqu'au terme 1/n^7 - n=100)*/
mpfr_set_ui(nn,n,GMP_RNDN);
mpfr_div(z,u,nn,GMP_RNDN);
mpfr_set(s2,z,GMP_RNDN);
mpfr_mul(z2,z,z,GMP_RNDN);
mpfr_div_2exp(v,z2,1,GMP_RNDN);
mpfr_sub(s2,s2,v,GMP_RNDN);
mpfr_set_ui(b,6,GMP_RNDN);
mpfr_mul(z,z,z2,GMP_RNDN);
mpfr_div(v,z,b,GMP_RNDN);
mpfr_add(s2,s2,v,GMP_RNDN);
mpfr_set_si(b,-30,GMP_RNDN);
mpfr_mul(z,z,z2,GMP_RNDN);
mpfr_div(v,z,b,GMP_RNDN);
mpfr_add(s2,s2,v,GMP_RNDN);
mpfr_set_si(b,42,GMP_RNDN);
mpfr_mul(z,z,z2,GMP_RNDN);
mpfr_div(v,z,b,GMP_RNDN);
mpfr_add(s2,s2,v,GMP_RNDN);
/*mpfr_print_raw(s2);printf("\n");
t=mpfr_out_str(stdout,10,0,s2,GMP_RNDN);printf("\n");*/
mpfr_add(s,s,s2,GMP_RNDN);
/*mpfr_print_raw(s);printf("\n");
t=mpfr_out_str(stdout,10,0,s,GMP_RNDN);printf("\n");*/
mpfr_add(s,s,u,GMP_RNDN);
/*mpfr_print_raw(s);printf("\n");
t=mpfr_out_str(stdout,10,0,s,GMP_RNDN);printf("\n");*/
/*Peut-on arrondir ? La reponse est oui*/
succes=mpfr_can_round(s,57,GMP_RNDN,GMP_RNDN,53);
mpfr_set(result,s,GMP_RNDN);
/*printf("%d\n",succes);*/
printf("Valeur de zeta(2) avec prec=53 et arrondi au plus pres:\n");
mpfr_print_raw(result);printf("\n");
t=mpfr_out_str(stdout,10,0,result,GMP_RNDN);printf("\n");
/*Test de comparaison avec pi^2/6*/
mpfr_pi(p,GMP_RNDN);
mpfr_mul(p,p,p,GMP_RNDN);
mpfr_set_ui(b,6,GMP_RNDN);
mpfr_div(p,p,b,GMP_RNDN);
/*mpfr_print_raw(p);printf("\n");
t=mpfr_out_str(stdout,10,0,p,GMP_RNDN);printf("\n");*/
succes=mpfr_can_round(p,63,GMP_RNDN,GMP_RNDN,53);
mpfr_set(res_p,p,GMP_RNDN);
/*printf("%d\n",succes);*/
printf("Valeur de pi^2/6 avec prec=53 et arrondi au plus pres:\n");
mpfr_print_raw(res_p);printf("\n");
t=mpfr_out_str(stdout,10,0,res_p,GMP_RNDN);printf("\n");

mpfr_clear(x);
mpfr_clear(y);
mpfr_clear(s);
mpfr_clear(s2);
mpfr_clear(u);
mpfr_clear(b);
mpfr_clear(v);
mpfr_clear(nn);
mpfr_clear(z);
mpfr_clear(z2);
mpfr_clear(p);
mpfr_clear(result);
mpfr_clear(res_p);
}
