/* functions which don't take as argument the rounding mode */
#define mpf_ceil mpfr_ceil
#define mpf_clear mpfr_clear
#define mpf_cmp mpfr_cmp
#define mpf_cmp_si mpfr_cmp_si
#define mpf_cmp_ui mpfr_cmp_ui
#define mpf_eq mpfr_eq
#define mpf_floor mpfr_floor
#define mpf_get_d mpfr_get_d
#define mpf_get_prec mpfr_get_prec
#define mpf_init mpfr_init
#define mpf_init2 mpfr_init2
#define mpf_random2 mpfr_random2
#define mpf_reldiff mpfr_reldiff
#define mpf_set_default_prec mpfr_set_default_prec
#define mpf_set_prec mpfr_set_prec
#define mpf_set_prec_raw mpfr_set_prec_raw
#define mpf_trunc mpfr_trunc
#define mpf_urandomb mpfr_urandomb
#define mpf_sgn mpfr_sgn

/* functions which take as argument the rounding mode */
#define mpf_abs((x),(y)) mpfr_abs((x),(y),__gmp_default_rounding_mode)
#define mpf_add((x),(y),(z)) mpfr_add((x),(y),(z),__gmp_default_rounding_mode) 
#define mpf_add_ui((x),(y),(z)) \
             mpfr_add_ui((x),(y),(z),__gmp_default_rounding_mode)
#define mpf_div((x),(y),(z)) mpfr_div((x),(y),(z),__gmp_default_rounding_mode) 
#define mpf_div_ui((x),(y),(z)) \
                          mpfr_div_ui((x),(y),(z),__gmp_default_rounding_mode) 
#define mpf_div_2exp((x),(y),(z)) \
                         mpfr_div2exp((x),(y),(z),__gmp_default_rounding_mode) 
#define mpf_dump((x),(y),(z)) \
                mpfr_dump((x),(y),(z),__gmp_default_rounding_mode)
#define mpf_get_str((x),(y),(z),(t),(u)) \
               mpfr_get_str((x),(y),(z),(t),(u),__gmp_default_rounding_mode) 
#define mpf_inp_str((x),(y),(z)) mpfr_inp_str((x),(y),(z),__gmp_default_rounding_mode)
#define mpf_set_str((x),(y),(z)) mpfr_set_str((x),(y),(z),__gmp_default_rounding_mode)
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
#define mpf_set((x),(y)) mpfr_set((x),(y),__gmp_default_rounding_mode) 
#define mpf_set_d((x),(y)) mpfr_set_d((x),(y),__gmp_default_rounding_mode) 
#define mpf_set_q((x),(y)) mpfr_set_q((x),(y),__gmp_default_rounding_mode) 
#define mpf_set_si((x),(y)) mpfr_set_si((x),(y),__gmp_default_rounding_mode) 
#define mpf_set_ui((x),(y),(z)) mpfr_set_ui((x),(y),__gmp_default_rounding_mode) 
#define mpf_set_z((x),(y)) mpfr_set_z((x),(y),__gmp_default_rounding_mode) 
#define mpf_sqrt((x),(y)) mpfr_sqrt((x),(y),__gmp_default_rounding_mode) 
#define mpf_sqrt_ui((x),(y)) mpfr_sqrt_ui((x),(y),__gmp_default_rounding_mode)
#define mpf_sub((x),(y),(z)) mpfr_sub((x),(y),(z),__gmp_default_rounding_mode) 
#define mpf_sub_ui((x),(y),(z)) mpfr_sub_ui((x),(y),(z),__gmp_default_rounding_mode) 
#define mpf_ui_div((x),(y),(z)) mpfr_ui_div((x),(y),(z),__gmp_default_rounding_mode)
#define mpf_ui_sub((x),(y),(z)) mpfr_ui_sub((x),(y),(z),__gmp_default_rounding_mode)


