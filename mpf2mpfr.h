#define mpf_abs((x),(y)) mpfr_abs((x),(y),__gmp_default_rounding_mode)
#define mpf_add((x),(y),(z)) mpfr_add((x),(y),(z),__gmp_default_rounding_mode) 

#define mpf_add_ui((x),(y),(z)) \
             mpfr_add_ui((x),(y),(z),__gmp_default_rounding_mode) /* TODO */

#define mpf_ceil((x),(y)) mpfr_ceil((x),(y),__gmp_default_rounding_mode) /* TODO */

#define mpf_clear((x)) mpfr_clear(x)
#define mpf_cmp((x),(y)) mpfr_cmp((x),(y)) 
#define mpf_cmp_si((x),(y)) mpfr_cmp_si((x),(y)) 
#define mpf_cmp_ui((x),(y)) mpfr_cmp_ui((x),(y)) 
#define mpf_div((x),(y),(z)) mpfr_div((x),(y),(z),__gmp_default_rounding_mode) 
#define mpf_div_ui((x),(y),(z)) \
                          mpfr_div_ui((x),(y),(z),__gmp_default_rounding_mode) 
#define mpf_div_2exp((x),(y),(z)) \
                         mpfr_div2exp((x),(y),(z),__gmp_default_rounding_mode) 

#define mpf_dump((x),(y),(z)) \
                mpfr_dump((x),(y),(z),__gmp_default_rounding_mode)
#define mpf_eq((x),(y),(z)) mpfr_eq((x),(y),(z),__gmp_default_rounding_mode)

#define mpf_floor((x),(y),(z)) \
                 mpfr_floor((x),(y),(z),__gmp_default_rounding_mode) /* TODO */

#define mpf_get_d((x)) mpfr_get_d((x))
#define mpf_get_prec((x)) mpfr_get_prec((x))
#define mpf_get_str((x),(y),(z),(t),(u)) \
               mpfr_get_str((x),(y),(z),(t),(u),__gmp_default_rounding_mode) 
#define mpf_init((x)) mpfr_init((x)) 
#define mpf_init2((x),(y)) mpfr_init2((x),(y)) 

#define mpf_inp_str((x),(y),(z)) mpfr_inp_str((x),(y),(z),__gmp_default_rounding_mode) /* TODO */

#define mpf_init_set((x),(y)) mpfr_init_set((x),(y),__gmp_default_rounding_mode) 
#define mpf_init_set_d((x),(y)) mpfr_init_set_d((x),(y),__gmp_default_rounding_mode) 
#define mpf_init_set_si((x),(y)) mpfr_init_set_si((x),(y),__gmp_default_rounding_mode) 
#define mpf_init_set_str((x),(y)) mpfr_init_set_str((x),(y),__gmp_default_rounding_mode) 
#define mpf_init_set_ui((x),(y)) mpfr_init_set_ui((x),(y),__gmp_default_rounding_mode) 
#define mpf_mul((x),(y),(z)) mpfr_mul((x),(y),(z),__gmp_default_rounding_mode) 
#define mpf_mul_2exp((x),(y),(z)) mpfr_mul_2exp((x),(y),(z),__gmp_default_rounding_mode) 
#define mpf_mul_ui((x),(y),(z)) mpfr_mul_ui((x),(y),(z),__gmp_default_rounding_mode) 
#define mpf_neg((x),(y)) mpfr_neg((x),(y),__gmp_default_rounding_mode) 
#define mpf_out_str((x),(y),(z),(t)) mpfr_out_str((x),(y),(z),(t),__gmp_default_rounding_mode) 
#define mpf_pow_ui((x),(y),(z)) mpfr_pow_ui((x),(y),(z),__gmp_default_rounding_mode) 
#define mpf_random2((x),(y),(z)) mpfr_random2((x),(y),(z),__gmp_default_rounding_mode)

#define mpf_reldiff((x),(y),(z)) mpfr_reldiff((x),(y),(z),__gmp_default_rounding_mode) /* TODO */

#define mpf_set((x),(y)) mpfr_set((x),(y),__gmp_default_rounding_mode) 
#define mpf_set_d((x),(y)) mpfr_set_d((x),(y),__gmp_default_rounding_mode) 
#define mpf_set_default_prec((x)) mpfr_set_default_prec((x)) 
#define mpf_set_prec((x),(y)) mpfr_set_prec((x),(y)) 
#define mpf_set_prec_raw((x),(y)) mpfr_set_prec_raw((x),(y),__gmp_default_rounding_mode)
#define mpf_set_q((x),(y)) mpfr_set_q((x),(y),__gmp_default_rounding_mode) 
#define mpf_set_si((x),(y)) mpfr_set_si((x),(y),__gmp_default_rounding_mode) 
#define mpf_set_ui((x),(y),(z)) mpfr_set_ui((x),(y),__gmp_default_rounding_mode) 
#define mpf_set_z((x),(y)) mpfr_set_z((x),(y),__gmp_default_rounding_mode) 
#define mpf_size((x)) mpfr_size((x),__gmp_default_rounding_mode)
#define mpf_sqrt((x),(y)) mpfr_sqrt((x),(y),__gmp_default_rounding_mode) 

#define mpf_sqrt_ui((x),(y)) mpfr_sqrt_ui((x),(y),__gmp_default_rounding_mode) /* TODO */

#define mpf_sub((x),(y),(z)) mpfr_sub_ui((x),(y),(z),__gmp_default_rounding_mode) 

#define mpf_trunc((x),(y)) mpfr_trunc((x),(y),__gmp_default_rounding_mode) /* TODO */

#define mpf_ui_div((x),(y),(z)) mpfr_ui_div((x),(y),(z),__gmp_default_rounding_mode) /* TODO */

#define mpf_ui_sub((x),(y),(z)) mpfr_ui_sub((x),(y),(z),__gmp_default_rounding_mode) /* TODO */

#define mpf_urandomb((x),(y),(z)) mpfr_urandomb((x),(y),(z),__gmp_default_rounding_mode)
#define mpf_sgn((x)) mpfr_sgn((x),__gmp_default_rounding_mode)


