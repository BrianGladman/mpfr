AUTOMAKE_OPTIONS = gnu

check_PROGRAMS = tversion tinits tsgn tcheck tisnan texceptions tset tabs tset_d tset_f tset_q tset_si tset_str tset_z tset_ld tadd tsub tmul tdiv tcmp tcmp2 tcan_round tround_prec tswap reuse tsqrt tnext tcomparisons teq tadd_ui tsub_ui tcmp_ui tdiv_ui tmul_ui tsqrt_ui tui_div tui_sub tcmp_d tmul_2exp tfma tfrac tget_str tout_str tget_d tget_d_2exp tconst_log2 tconst_pi tconst_euler trandom ttrunc trint texp texp2 texpm1 tlog tlog2 tlog10 tlog1p tpow tui_pow tpow3 tcosh tsinh ttanh tacosh tasinh tatanh thyperbolic tasin tacos tcos tatan tsin ttan tsin_cos tagm thypot tfactorial tgamma terf tcbrt mpf_compat mpfr_compat tzeta tsub1sp tgmpop tsum

EXTRA_DIST = tgeneric.c mpf_compat.h

LDADD = libfrtests.a $(MPFR_LIBM) $(top_builddir)/libmpfr.a  @LDADD@
INCLUDES = -I$(top_srcdir) @INCLUDES@

check_LIBRARIES = libfrtests.a
libfrtests_a_SOURCES = memory.c rnd_mode.c tests.c cmp_str.c
EXTRA_libfrtests_a_SOURCES = amd64.asm x86.asm
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
