bin_PROGRAMS=  tabs.c tadd.c tagm.c tcan_round.c tcmp.c tcmp_ui.c 
tdiv.c tdiv_ui.c texp.c  tget_str.c tlog.c 
tlog2.c tmul.c tmul_2exp.c tmul_ui.c tout_str.c tpi.c tround.c tset_d.c
tset_f.c tset_q.c tset_si.c tset_str.c tset_z.c tsin_cos.c tsqrt.c
tsqrt_ui.c tui_div.c tui_sub.c tzeta.c
 
LDADD=  ../libmpfr.a  @LDADD@
INCLUDES=@INCLUDES@ -I$(top_srcdir)

tests: 	all
	for i in $(bin_PROGRAMS); do echo Testing $$i ; ./$$i; done





