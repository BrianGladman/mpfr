CC=gcc
INCFLAGS=-I$(GMPDIR)/include
LDFLAGS=$(GMPDIR)/lib/libgmp.a -lm
CFLAGS=-g -Wall -DMsb $(special) -D$(dirname)
ODIR=o.$(dirname)

MPFR_OBJECTS = o.$(dirname)/add.o o.$(dirname)/clear.o o.$(dirname)/cmp.o o.$(dirname)/init.o o.$(dirname)/mul.o o.$(dirname)/print_raw.o o.$(dirname)/rnd_mode.o o.$(dirname)/round.o o.$(dirname)/set_d.o o.$(dirname)/set_prec.o o.$(dirname)/set_dfl_prec.o o.$(dirname)/set_dfl_rnd.o o.$(dirname)/neg.o o.$(dirname)/set_f.o o.$(dirname)/set_si.o o.$(dirname)/set_str_raw.o o.$(dirname)/mul_ui.o o.$(dirname)/mul_2exp.o o.$(dirname)/div_2exp.o o.$(dirname)/sub.o o.$(dirname)/set.o o.$(dirname)/div.o o.$(dirname)/sqrt.o o.$(dirname)/agm.o o.$(dirname)/get_str.o o.$(dirname)/cmp_ui.o

dft target::
	@echo "Utiliser soit make all, soit ./mmpfr"
	@echo "./mmpfr peut prendre les options tests ou clean"

libmpfr.a: mpfr.h $(MPFR_OBJECTS)
	$(AR) cr o.$(dirname)/libmpfr.a $(MPFR_OBJECTS)
	-chmod g+w o.$(dirname)/*.o o.$(dirname)/libmpfr.a
	-ranlib o.$(dirname)/libmpfr.a

clean: 
	-rm o.$(dirname)/*.o *~ o.$(dirname)/*.a core

dist:
	-rm mpfr.tgz
	tar cvf mpfr.tar Makefile *.c *.h tests/*.c tests/Makefile
	gzip -9 mpfr.tar

o.$(dirname)/add.o: mpfr.h add.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c add.c -o o.$(dirname)/add.o

o.$(dirname)/sub.o: mpfr.h sub.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c sub.c -o o.$(dirname)/sub.o

o.$(dirname)/set.o: mpfr.h set.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c set.c -o o.$(dirname)/set.o

o.$(dirname)/neg.o: mpfr.h neg.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c neg.c -o o.$(dirname)/neg.o

o.$(dirname)/clear.o: mpfr.h clear.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c clear.c -o o.$(dirname)/clear.o

o.$(dirname)/cmp.o: mpfr.h cmp.c
	-$(CC) $(CFLAGS) $(INCFLAGS) -c cmp.c -o o.$(dirname)/cmp.o

o.$(dirname)/cmp_ui.o: mpfr.h cmp_ui.c
	-$(CC) $(CFLAGS) $(INCFLAGS) -c cmp_ui.c -o o.$(dirname)/cmp_ui.o

o.$(dirname)/init.o: mpfr.h init.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c init.c -o o.$(dirname)/init.o

o.$(dirname)/mul.o: mpfr.h mul.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c mul.c -o o.$(dirname)/mul.o

o.$(dirname)/print_raw.o: mpfr.h print_raw.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c print_raw.c -o o.$(dirname)/print_raw.o

o.$(dirname)/rnd_mode.o: mpfr.h rnd_mode.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c rnd_mode.c -o o.$(dirname)/rnd_mode.o

o.$(dirname)/set_f.o: mpfr.h set_f.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c set_f.c -o o.$(dirname)/set_f.o

o.$(dirname)/set_si.o: mpfr.h set_si.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c set_si.c -o o.$(dirname)/set_si.o

o.$(dirname)/round.o: mpfr.h round.c
	-$(CC) $(CFLAGS) $(INCFLAGS) -c round.c -o o.$(dirname)/round.o

o.$(dirname)/set_d.o: mpfr.h set_d.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c set_d.c -o o.$(dirname)/set_d.o

o.$(dirname)/set_prec.o: mpfr.h set_prec.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c set_prec.c -o o.$(dirname)/set_prec.o

o.$(dirname)/set_str_raw.o: mpfr.h set_str_raw.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c set_str_raw.c -o o.$(dirname)/set_str_raw.o

o.$(dirname)/get_str.o: mpfr.h get_str.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c get_str.c -o o.$(dirname)/get_str.o

o.$(dirname)/mul_2exp.o: mpfr.h mul_2exp.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c mul_2exp.c -o o.$(dirname)/mul_2exp.o

o.$(dirname)/div_2exp.o: mpfr.h div_2exp.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c div_2exp.c -o o.$(dirname)/div_2exp.o

o.$(dirname)/div.o: mpfr.h div.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c div.c -o o.$(dirname)/div.o

o.$(dirname)/sqrt.o: mpfr.h sqrt.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c sqrt.c -o o.$(dirname)/sqrt.o

o.$(dirname)/set_dfl_prec.o: mpfr.h set_dfl_prec.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c set_dfl_prec.c -o o.$(dirname)/set_dfl_prec.o

o.$(dirname)/set_dfl_rnd.o: mpfr.h set_dfl_rnd.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c set_dfl_rnd.c -o o.$(dirname)/set_dfl_rnd.o

o.$(dirname)/mul_ui.o: mpfr.h mul_ui.c
	$(CC) $(CFLAGS) $(INCFLAGS) -c mul_ui.c -o o.$(dirname)/mul_ui.o

o.$(dirname)/agm.o: mpfr.h agm.c o.$(dirname)/div.o
	$(CC) $(CFLAGS) $(INCFLAGS) -c agm.c -o o.$(dirname)/agm.o

tadd: tests/tadd.c o.$(dirname)/add.o libmpfr.a
	$(CC) $(CFLAGS) $(INCFLAGS) -I. -o tests/tadd tests/tadd.c $(ODIR)/libmpfr.a $(LDFLAGS)

tdiv: tests/tdiv.c o.$(dirname)/div.o libmpfr.a
	$(CC) $(CFLAGS) $(INCFLAGS) -I. -o tests/tdiv tests/tdiv.c $(ODIR)/libmpfr.a $(LDFLAGS)

tsqrt: tests/tsqrt.c o.$(dirname)/sqrt.o libmpfr.a
	$(CC) $(CFLAGS) $(INCFLAGS) -I. -o tests/tsqrt tests/tsqrt.c $(ODIR)/libmpfr.a $(LDFLAGS)

tcmp: tests/tcmp.c o.$(dirname)/add.o libmpfr.a
	$(CC) $(CFLAGS) $(INCFLAGS) -I. -o tests/tcmp tests/tcmp.c $(ODIR)/libmpfr.a $(LDFLAGS)

tcmp_ui: tests/tcmp_ui.c libmpfr.a
	$(CC) $(CFLAGS) $(INCFLAGS) -I. -o tests/tcmp_ui tests/tcmp_ui.c $(ODIR)/libmpfr.a $(LDFLAGS)

tcmp2: tests/tcmp2.c o.$(dirname)/add.o libmpfr.a
	$(CC) $(CFLAGS) $(INCFLAGS) -I. -o tests/tcmp2 tests/tcmp2.c $(ODIR)/libmpfr.a $(LDFLAGS)

tmul: tests/tmul.c o.$(dirname)/mul.o mul.c libmpfr.a
	$(CC) $(CFLAGS) $(INCFLAGS) -I. -o tests/tmul tests/tmul.c $(ODIR)/libmpfr.a $(LDFLAGS)

tround: tests/tround.c round.c $(ODIR)/libmpfr.a
	$(CC) $(CFLAGS) $(INCFLAGS) -I. -o tests/tround tests/tround.c $(ODIR)/libmpfr.a $(LDFLAGS)

tset_f: tests/tset_f.c set_f.c $(ODIR)/libmpfr.a
	$(CC) $(CFLAGS) $(INCFLAGS) -I. -o tests/tset_f tests/tset_f.c $(ODIR)/libmpfr.a $(LDFLAGS)

tset_d: tests/tset_d.c set_d.c libmpfr.a
	$(CC) $(CFLAGS) $(INCFLAGS) -I. -o tests/tset_d tests/tset_d.c $(ODIR)/libmpfr.a $(LDFLAGS)

tset_si: tests/tset_si.c set_si.c $(ODIR)/libmpfr.a
	$(CC) $(CFLAGS) $(INCFLAGS) -I. -o tests/tset_si tests/tset_si.c $(ODIR)/libmpfr.a $(LDFLAGS)

tmul_ui: tests/tmul_ui.c mul_ui.c $(ODIR)/libmpfr.a
	$(CC) $(CFLAGS) $(INCFLAGS) -I. -o tests/tmul_ui tests/tmul_ui.c $(ODIR)/libmpfr.a $(LDFLAGS)

tset_str: tests/tset_str.c set_str_raw.c $(ODIR)/libmpfr.a
	$(CC) $(CFLAGS) $(INCFLAGS) -I. -o tests/tset_str tests/tset_str.c $(ODIR)/libmpfr.a $(LDFLAGS)

tget_str: tests/tget_str.c o.$(dirname)/get_str.o libmpfr.a
	$(CC) $(CFLAGS) $(INCFLAGS) -I. -o tests/tget_str tests/tget_str.c $(ODIR)/libmpfr.a $(LDFLAGS)

tmul_2exp: tests/tmul_2exp.c mul_2exp.c $(ODIR)/libmpfr.a
	$(CC) $(CFLAGS) $(INCFLAGS) -I. -o tests/tmul_2exp tests/tmul_2exp.c $(ODIR)/libmpfr.a $(LDFLAGS)

tagm: tests/tagm.c agm.c libmpfr.a
	$(CC) $(CFLAGS) $(INCFLAGS) -I. -o tests/tagm tests/tagm.c $(ODIR)/libmpfr.a $(LDFLAGS)


all: 
	./mmpfr; 

