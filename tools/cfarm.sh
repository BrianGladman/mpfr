#!/bin/bash
# this script helps to check MPFR on the GCC compile farm
# 1) update the GMP version if needed
# 2) update the MPFR release candidate
# 3) ssh gcc10 < cfarm.sh
GMP=gmp-6.1.2
MPFR=mpfr-3.1.6
RC=rc1
/bin/rm -fr gmp*
if [ ! -d $GMP ]; then
   wget ftp://ftp.gmplib.org/pub/$GMP/$GMP.tar.bz2
   tar jxf $GMP.tar.bz2
   cd $GMP
   ./configure --prefix=$HOME
   make -j4
   make install
   cd $HOME
fi
/bin/rm -fr mpfr*
wget http://www.mpfr.org/$MPFR/$MPFR-$RC.tar.gz
gunzip $MPFR-$RC.tar.gz
tar xf $MPFR-$RC.tar
cd $MPFR-$RC
if [ "`hostname`" = "power-aix" ]; then # gcc111
   export OBJECT_MODE=64
   # or ./configure AR="ar -X64" NM="nm -B -X64"
fi
./configure --with-gmp=$HOME
make -j4
make check -j4

# results with mpfr-3.1.6-rc1.tar.gz (160 tests)
# gcc10 # PASS:  159 # SKIP:  1
# gcc11 Connection refused (asks for a password)
# gcc12 # PASS:  159 # SKIP:  1
# gcc13 Connection timed out
# gcc14 Connection timed out
# gcc15 # PASS:  159 # SKIP:  1
# gcc16 # PASS:  159 # SKIP:  1
# gcc17 Connection timed out
# gcc20 Connection timed out
# gcc21 Network is unreachable
# gcc22 # PASS:  159 # SKIP:  1
# gcc23 # PASS:  159 # SKIP:  1
# gcc24 # PASS:  159 # SKIP:  1
# gcc33 Connection timed out
# gcc34 Connection timed out
# gcc35 Connection timed out
# gcc36 Connection timed out
# gcc37 Connection timed out
# gcc38 Connection timed out
# gcc40 Connection timed out
# gcc41 Connection timed out
# gcc42 Connection timed out
# gcc43 Connection timed out
# gcc45 # PASS:  159 # SKIP:  1
# gcc46 Connection timed out
# gcc47 Connection timed out
# gcc49 Name or service not known
# gcc50 Connection timed out
# gcc51 Connection timed out
# gcc52 Connection timed out
# gcc53 Connection timed out
# gcc54 Connection timed out
# gcc55 Connection timed out
# gcc56 Connection timed out
# gcc57 Connection timed out
# gcc60 Connection timed out
# gcc61 Connection timed out
# gcc62 Connection timed out
# gcc63 Connection timed out
# gcc64 Connection timed out
# gcc66 Connection timed out
# gcc67 # PASS:  159 # SKIP:  1
# gcc68 # PASS:  159 # SKIP:  1
# gcc70 Connection timed out
# gcc75 # PASS:  159 # SKIP:  1 
# gcc76 No route to host
# gcc100 Name or service not known
# gcc101 Name or service not known
# gcc110 # PASS:  159 # SKIP:  1
# gcc111 # PASS:  159 # SKIP:  1
# gcc112 # PASS:  159 # SKIP:  1
# gcc113 # PASS:  159 # SKIP:  1
# gcc114 # PASS:  159 # SKIP:  1 
# gcc115 # PASS:  159 # SKIP:  1
# gcc116 # PASS:  159 # SKIP:  1
# gcc117 Connection timed out
# gcc118 # PASS:  159 # SKIP:  1
# gcc119 # PASS:  159 # SKIP:  1 (manual compilation)
# gcc200 Connection timed out
# gcc201 Connection timed out
# gcc202 # PASS:  159 # SKIP:  1 (gmp-6.1.2 configured with --disable-assembly)
