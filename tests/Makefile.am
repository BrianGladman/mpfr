AUTOMAKE_OPTIONS = gnu

check_PROGRAMS = tversion tinits tsgn tcheck tisnan texceptions tset_exp tset tabs tset_d tset_f tset_q tset_si tset_str tset_z tset_ld tset_sj tswap tcopysign tcmp tcmp2 tcmpabs tcmp_d tcmp_ld tcomparisons teq tadd tsub tmul tdiv tsub1sp tadd1sp tadd_ui tsub_ui tcmp_ui tdiv_ui tmul_ui tsqrt_ui tui_div tui_sub tgmpop tsi_op tmul_2exp tfma tsum tdim tminmax tnext tfits tget_d tget_d_2exp tget_z tget_str tget_sj tout_str tinp_str toutimpl tcan_round tround_prec tsqrt tconst_log2 tconst_pi tconst_euler trandom ttrunc trint tfrac texp texp2 texpm1 tlog tlog2 tlog10 tlog1p tpow tui_pow tpow3 tcosh tsinh ttanh tacosh tasinh tatanh thyperbolic tasin tacos tcos tatan tsin ttan tsin_cos tagm thypot tfactorial tgamma terf tcbrt tzeta mpf_compat mpfr_compat reuse tsqr tstrtofr tpow_z

EXTRA_DIST = tgeneric.c mpf_compat.h inp_str.data

LDADD = libfrtests.a $(MPFR_LIBM) $(top_builddir)/libmpfr.a
INCLUDES = -I$(top_srcdir) @INCLUDES@
LOADLIBES=$(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS) $(LDADD) $(LIBS)

check_LIBRARIES = libfrtests.a
libfrtests_a_SOURCES = memory.c rnd_mode.c tests.c cmp_str.c

$(top_builddir)/libmpfr.a:
	cd $(top_builddir); make libmpfr.a

TESTS = $(check_PROGRAMS)
