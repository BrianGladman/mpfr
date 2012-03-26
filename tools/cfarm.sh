#!/bin/bash
# this script helps to check MPFR is running on the GCC compile farm
# 1) update the GMP version if needed
# 2) update the MPFR release candidate
# 3) ssh gcc100 < cfarm.sh
if [ ! -d gmp-5.0.4 ]; then
   /bin/rm -fr gmp*
   wget ftp://ftp.gmplib.org/pub/gmp-5.0.4/gmp-5.0.4.tar.bz2
   tar jxf gmp-5.0.4.tar.bz2
   cd gmp-5.0.4
   if [ "`hostname`" = "dingo" ]; then
      # gcc 4.3.2 miscompiles GMP
      ./configure --prefix=$HOME CC=/opt/cfarm/release/4.4.1/bin/gcc
   else
      ./configure --prefix=$HOME
   fi
   make
   make install
   cd $HOME
fi
/bin/rm -fr mpfr*
wget http://www.loria.fr/~zimmerma/mpfr-3.2.0-dev.tar.bz2
tar jxf mpfr-3.2.0-dev.tar.bz2
cd  mpfr-3.2.0-dev
if [ "`hostname`" = "dingo" ]; then
   # http://websympa.loria.fr/wwsympa/arc/mpfr/2011-10/msg00048.html
   ./configure --with-gmp=$HOME --disable-thread-safe
else
   ./configure --with-gmp=$HOME
fi
make
make check

# results with mpfr-3.1.0-rc1.tar.bz2
# gcc10 All 160 tests passed (1 test was not run)
# gcc11 All 160 tests passed (1 test was not run)
# gcc12 All 160 tests passed (1 test was not run)
# gcc13 All 160 tests passed (1 test was not run)
# gcc14 All 160 tests passed (1 test was not run)
# gcc15 All 160 tests passed (1 test was not run)
# gcc16 All 160 tests passed (1 test was not run)
# gcc17 All 160 tests passed (1 test was not run)
# gcc20 All 160 tests passed (1 test was not run)
# gcc33 Connection refused
# gcc34 Connection refused
# gcc35 Connection refused
# gcc36 Connection refused
# gcc37 Connection refused
# gcc38 All 160 tests passed (1 test was not run)
# gcc40 Connection timed out
# gcc41 Connection refused
# gcc42 Connection timed out
# gcc43 Connection refused
# gcc45 Connection refused
# gcc46 All 160 tests passed (1 test was not run)
# gcc47 All 160 tests passed (1 test was not run)
# gcc48 All 160 tests passed (1 test was not run)
# gcc50 Connection refused
# gcc51 Connection refused
# gcc52 Connection timed out
# gcc53 Connection timed out
# gcc54 1 of 160 tests failed FAIL: tagm ok with --disable-thread-safe or -O1
# gcc55 Connection timed out
# gcc56 Connection refused
# gcc57 Connection refused
# gcc60 All 160 tests passed (1 test was not run)
# gcc61 All 160 tests passed (1 test was not run) ABI=1.0 GCC 4.4.1
# gcc62 1 of 160 tests failed FAIL: tagm
# gcc63 1 of 160 tests failed FAIL: tagm
# gcc64 All 160 tests passed (1 test was not run)
# gcc70 All 160 tests passed (1 test was not run)
# gcc100 Connection timed out
# gcc101 Connection timed out
# gcc200 Connection timed out
# gcc201 Connection timed out

# results with mpfr-3.1.0-rc2.tar.bz2
# gcc54 All 160 tests passed (1 test was not run)
# gcc61 All 160 tests passed (1 test was not run) ABI=1.0 GCC 4.4.1
# gcc62 All 160 tests passed (1 test was not run)
# gcc63 All 160 tests passed (1 test was not run)
