SUBDIRS=tests
EXTRA_DIST=generic.c

include_HEADERS = mpfr.h mpf2mpfr.h mpfr-test.h

lib_LIBRARIES = libmpfr.a

libmpfr_a_SOURCES= cputime.h mpfr.h mpfr-impl.h sin_cos.c extract.c add.c add_ui.c add_ulp.c agm.c clear.c cmp.c cmp_ui.c div_2exp.c div.c div_ui.c dump.c eq.c exp2.c exp3.c exp.c get_str.c init.c inp_str.c log2.c log.c mul_2exp.c mul.c mul_ui.c neg.c out_str.c pi.c pow.c print_raw.c print_rnd_mode.c random2.c random.c reldiff.c rnd_mode.c round.c set.c set_d.c set_dfl_prec.c set_rnd.c set_f.c set_prc_raw.c set_prec.c set_q.c set_si.c set_str.c set_str_raw.c set_z.c sqrt.c sqrt_ui.c sub.c sub_ui.c trunc.c ui_div.c ui_sub.c urandomb.c sqrtrem.c mpz_set_fr.c swap.c
libmpfr_a_LIBADD= ceil.o floor.o


ceil.o: trunc.c mpfr.h
	$(COMPILE) -DMPFR_CEIL -c $(srcdir)/trunc.c -o ceil.o
floor.o: trunc.c  mpfr.h
	$(COMPILE) -DMPFR_FLOOR -c $(srcdir)/trunc.c -o floor.o
trunc.o: trunc.c  mpfr.h
	$(COMPILE) -DMPFR_TRUNC -c @MISCFLAGS@

info_TEXINFOS = mpfr.texi

tests: 
	cd tests ; make tests
