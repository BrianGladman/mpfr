check_PROGRAMS= reuse tabs tadd tagm tcan_round tcmp tcmp2 tcmp_ui tdiv tdiv_ui tdump teq texp tget_str tlog tconst_log2 tmul tmul_2exp tmul_ui tout_str tconst_pi tpow trandom tround_prec tset_d tset_f tset_q tset_si tset_str tset_z tsqrt tsqrt_ui tui_div tui_sub tswap ttrunc trint tisnan tget_d tatan tcosh tsinh ttanh tacosh tasinh tatanh thyperbolic texp2 tfactorial tsub tasin tconst_euler tcos tsin ttan tsub_ui tset tlog1p texpm1 tlog2 tlog10 tui_pow tpow3 tadd_ui texceptions tfma thypot tacos
EXTRA_DIST= tgeneric.c

LDADD= $(top_builddir)/libmpfr.a  @LDADD@
INCLUDES=-I$(top_srcdir) @INCLUDES@

../libmpfr.a:
	cd $(top_srcdir); make libmpfr.a  

check: 	all
	for i in $(check_PROGRAMS); do echo Testing $$i ; ./$$i || echo "### $$i failed ###"; done

