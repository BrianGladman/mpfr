AUTOMAKE_OPTIONS = gnu

check_PROGRAMS = reuse tabs tadd tagm tcan_round tcmp tcmp2 tcmp_ui tdiv tdiv_ui tdump teq texp tfrac tget_str tlog tconst_log2 tmul tmul_2exp tmul_ui tout_str tconst_pi tpow trandom tround_prec tset_d tset_f tset_q tset_si tset_str tset_z tsqrt tsqrt_ui tui_div tui_sub tswap ttrunc trint tisnan tget_d tatan tcosh tsinh ttanh tacosh tasinh tatanh thyperbolic texp2 tfactorial tsub tasin tconst_euler tcos tsin ttan tsub_ui tset tlog1p texpm1 tlog2 tlog10 tui_pow tpow3 tadd_ui texceptions tfma thypot tacos tgamma tset_ld tcbrt tsin_cos
EXTRA_DIST = tgeneric.c

LDADD = libfrtests.a $(MPFR_LIBM) $(top_builddir)/libmpfr.a  @LDADD@
INCLUDES = -I$(top_srcdir) @INCLUDES@

check_LIBRARIES = libfrtests.a
libfrtests_a_SOURCES = memory.c rnd_mode.c tests.c
EXTRA_libfrtests_a_SOURCES = x86.asm
libfrtests_a_DEPENDENCIES = $(TESTS_ASM_OBJECTS)
libfrtests_a_LIBADD = $(libfrtests_a_DEPENDENCIES)

$(top_builddir)/libmpfr.a:
	cd $(top_builddir); make libmpfr.a

TESTS = $(check_PROGRAMS)


# This $(srcdir) and $(CYGPATH_W) stuff follows automake 1.7.2
SUFFIXES = .asm
.asm.o:
	grep -v '^#' `test -f '$<' || echo '$(srcdir)/'`$< >tmp-$*.s
	$(COMPILE) -c tmp-$*.s -o $@
	rm -f tmp-$*.s
.asm.obj:
	grep -v '^#' `if test -f '$<'; then $(CYGPATH_W) '$<'; else $(CYGPATH_W) '$(srcdir)/$<' >tmp-$*.s
	$(COMPILE) -c tmp-$*.s -o $@
	rm -f tmp-$*.s
