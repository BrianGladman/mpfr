CC=gcc
INCFLAGS=-I$(GMPDIR)/include -I..
LDFLAGS=libmpfr.a $(GMPDIR)/lib/libgmp.a 
CFLAGS=-g -DMsb

all_files=tmul tround tset_f tcmp tcmp2 tdiv tmul_2exp tmul_ui tset_d tset_i tset_str tadd tsqrt tagm

all: $(all_files)
	for f in $(all_files) ; do \
	   echo Testing $$f; \
	   ./$$f; \
	done

../$(ODIR)/libmpfr.a: ../mpfr.h
	cd .. ; ./mmpfr

tadd: tadd.c ../add.c
	cd .. ; ./mmpfr tadd
	-chmod g+w tadd

tcmp: tcmp.c ../cmp.c
	cd .. ; ./mmpfr tcmp
	-chmod g+w tcmp

tcmp2: tcmp2.c ../cmp.c
	cd .. ; ./mmpfr tcmp2
	-chmod g+w tcmp2

tmul: tmul.c ../mul.c
	cd .. ; ./mmpfr tmul
	-chmod g+w tmul

tdiv: tdiv.c ../div.c
	cd .. ; ./mmpfr tdiv
	-chmod g+w tdiv

tround: tround.c ../round.c
	cd .. ; ./mmpfr tround
	-chmod g+w tround

tset_f: tset_f.c ../set_f.c
	cd .. ; ./mmpfr tset_f
	-chmod g+w tset_f

tset_d: tset_d.c ../set_d.c
	cd .. ; ./mmpfr tset_d
	-chmod g+w tset_d

tset_si: tset_si.c ../set_si.c
	cd .. ; ./mmpfr tset_si
	-chmod g+w tset_si

tset_str: tset_str.c ../set_str_raw.c
	cd .. ; ./mmpfr tset_str
	-chmod g+w tset_str

tmul_ui: tmul_ui.c ../mul_ui.c
	cd .. ; ./mmpfr tmul_ui
	-chmod g+w tmul_ui

tmul_2exp: tmul_2exp.c ../mul_2exp.c
	cd .. ; ./mmpfr tmul_2exp
	-chmod g+w tmul_2exp

tagm: tagm.c ../agm.c
	cd .. ; ./mmpfr tagm
	-chmod g+w tagm

clean: 
	-rm *~ $(all_files)

