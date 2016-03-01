#!/bin/bash
# this script helps to check MPFR on the GCC compile farm
# 1) update the GMP version if needed
# 2) update the MPFR release candidate
# 3) ssh gcc10 < cfarm.sh
if [ ! -d gmp-6.1.0 ]; then
   /bin/rm -fr gmp*
   wget ftp://ftp.gmplib.org/pub/gmp-6.1.0/gmp-6.1.0.tar.bz2
   tar jxf gmp-6.1.0.tar.bz2
   cd gmp-6.1.0
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
wget http://www.mpfr.org/mpfr-3.1.4/mpfr-3.1.4-rc2.tar.gz
gunzip mpfr-3.1.4-rc2.tar.gz
tar xf mpfr-3.1.4-rc2.tar
cd mpfr-3.1.4-rc2
./configure --with-gmp=$HOME
make
make check

# results with mpfr-3.1.4-rc2.tar.gz (160 tests)
# gcc10 # PASS:  159 # SKIP:  1
# gcc11 Connection refused (asks for a password)
# gcc12 # PASS:  159 # SKIP:  1
# gcc13 # PASS:  159 # SKIP:  1
# gcc14 # PASS:  159 # SKIP:  1
# gcc15 # PASS:  159 # SKIP:  1
# gcc16 # PASS:  159 # SKIP:  1
# gcc17 Connection timed out
# gcc20 # PASS:  159 # SKIP:  1
# gcc21 # PASS:  159 # SKIP:  1
# gcc33 Connection timed out
# gcc34 Connection refused
# gcc35 Connection refused
# gcc36 Connection timed out
# gcc37 Connection refused
# gcc38 Connection refused
# gcc40 Connection refused
# gcc41 Connection refused
# gcc42 Connection refused
# gcc43 Connection refused
# gcc45 # PASS:  159 # SKIP:  1
# gcc46 Connection timed out
# gcc47 Connection timed out
# gcc49 Name or service not known
# gcc50 Connection timed out
# gcc51 Connection refused
# gcc52 Connection timed out
# gcc53 Connection timed out
# gcc54 Connection refused
# gcc55 Connection refused
# gcc56 Connection refused
# gcc57 Connection refused
# gcc60 Connection refused
# gcc61 Connection timed out
# gcc62 Connection refused
# gcc63 Connection timed out
# gcc64 Connection timed out
# gcc66 Connection refused
# gcc70 # PASS:  159 # SKIP:  1
# gcc75 # PASS:  159 # SKIP:  1 
# gcc76 # PASS:  159 # SKIP:  1
# gcc100 Connection timed out
# gcc101 Connection timed out
# gcc110 # PASS:  159 # SKIP:  1
# gcc111 configure: error: could not determine ar interface
# gcc112 # PASS:  159 # SKIP:  1
# gcc113 # PASS:  159 # SKIP:  1
# gcc114 # PASS:  159 # SKIP:  1 
# gcc115 # PASS:  159 # SKIP:  1
# gcc116 # PASS:  159 # SKIP:  1
# gcc200 Connection timed out
# gcc201 Connection timed out
