check_PROGRAMS= reuse tabs tadd tagm tcan_round tcmp tcmp2 tcmp_ui tdiv tdiv_ui tdump teq texp tget_str tlog tlog2 tmul tmul_2exp tmul_ui tout_str tpi tpow trandom tround tset_d tset_f tset_q tset_si tset_str tset_z tsqrt tsqrt_ui tui_div tui_sub tsin_cos tswap ttrunc tisnan tget_d tarctan
 
LDADD= $(top_builddir)/libmpfr.a  @LDADD@
INCLUDES=@INCLUDES@ -I$(top_srcdir)

../libmpfr.a:
	cd $(top_srcdir); make libmpfr.a  

check: 	all
	for i in $(check_PROGRAMS); do echo Testing $$i ; ./$$i || echo "### $$i failed ###"; done





