SUBDIRS=tests
EXTRA_DIST=generic.c COPYING.LIB

include_HEADERS = mpfr.h mpf2mpfr.h mpfr-test.h

lib_LIBRARIES = libmpfr.a

libmpfr_a_SOURCES= cputime.h mpfr.h mpfr-impl.h exceptions.c save_expo.c extract.c add.c add1.c add_one_ulp.c add_ui.c agm.c clear.c cmp.c cmp_abs.c cmp_ui.c div_2exp.c div_2si.c div_2ui.c div.c div_ui.c dump.c eq.c exp2.c exp3.c exp.c get_str.c init.c inp_str.c isinteger.c isinf.c isnan.c isnum.c const_log2.c log.c mul_2exp.c mul_2si.c mul_2ui.c mul.c mul_ui.c neg.c out_str.c const_pi.c pow.c pow_si.c pow_ui.c print_raw.c print_rnd_mode.c random2.c random.c reldiff.c rnd_mode.c round.c set.c set_d.c set_dfl_prec.c set_rnd.c set_f.c set_prc_raw.c set_prec.c set_q.c set_si.c set_str.c set_str_raw.c set_ui.c set_z.c sqrt.c sqrt_ui.c sub.c sub1.c sub_one_ulp.c sub_ui.c trunc.c ui_div.c ui_sub.c urandomb.c sqrtrem.c mpz_set_fr.c swap.c  factorial.c cosh.c sinh.c tanh.c acosh.c asinh.c atanh.c atan.c cmp2.c exp_2.c asin.c const_euler.c cos.c sin.c tan.c fma.c hypot.c log1p.c expm1.c log2.c log10.c ui_pow.c ui_pow_ui.c minmax.c dim.c copysign.c gmp_op.c init2.c acos.c

libmpfr_a_LIBADD= ceil.o floor.o

ceil.o: trunc.c mpfr.h
	$(COMPILE) -DMPFR_CEIL -c $(srcdir)/trunc.c -o ceil.o
floor.o: trunc.c  mpfr.h
	$(COMPILE) -DMPFR_FLOOR -c $(srcdir)/trunc.c -o floor.o
trunc.o: trunc.c  mpfr.h
	$(COMPILE) -DMPFR_TRUNC -c @MISCFLAGS@

info_TEXINFOS = mpfr.texi

DISTCLEANFILES = mpfr.info

tests: 
	cd tests ; make tests




