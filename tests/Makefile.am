bin_PROGRAMS=  tabs tadd tagm tcan_round tcmp tcmp_ui tdiv tdiv_ui texp  tget_str tlog tlog2 tmul tmul_2exp tmul_ui tout_str tpi tround tset_d tset_f tset_q tset_si tset_str tset_z  tsqrt tsqrt_ui tui_div tui_sub tzeta
 
LDADD=  $(top_srcdir)/libmpfr.a  @LDADD@
INCLUDES=@INCLUDES@ -I$(top_srcdir)

../libmpfr.a:
	cd $(top_srcdir); make libmpfr.a  

tests: 	all
	for i in $(bin_PROGRAMS); do echo Testing $$i ; ./$$i; done





