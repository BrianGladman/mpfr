#define mpfr_init_set_si(x, i, p, rnd) {\
 mpfr_init2((x), (p)); \
 mpfr_set_si((x), (i), (rnd)); \
}

#define mpfr_init_set_ui(x, i, p, rnd) {\
 mpfr_init2((x), (p)); \
 mpfr_set_ui((x), (i), (rnd)); \
}

#define mpfr_init_set_d(x, d, p, rnd) {\
 mpfr_init2((x), (p)); \
 mpfr_set_d((x), (d), (rnd)); \
}

#define mpfr_init_set(x, y, p, rnd) {\
 mpfr_init2((x), (p)); \
 mpfr_set((x), (y), (rnd)); \
}

#define mpfr_init_set_f(x, y, p, rnd) {\
 mpfr_init2((x), (p)); \
 mpfr_set_f((x), (y), (rnd)); \
}

#define mpfr_init_set_str(x, y, p, rnd) {\
 mpfr_init2((x), (p)); \
 mpfr_set_str((x), (y), (rnd)); \
}

#define mpfr_init_set_str_raw(x, y, p, rnd) {\
 mpfr_init2((x), (p)); \
 mpfr_set_str_raw((x), (y), (rnd)); \
}

