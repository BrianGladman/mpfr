# this script helps to check MPFR is running on the GCC compile farm
# 1) update the GMP version if needed
# 2) update the MPFR release candidate
# 3) ssh gcc100 < cfarm.sh
/bin/rm -fr gmp-5.0.1*
wget ftp://ftp.gmplib.org/pub/gmp-5.0.1/gmp-5.0.1.tar.bz2
tar jxf gmp-5.0.1.tar.bz2
cd gmp-5.0.1
./configure --prefix=$HOME --disable-shared
make
make install
cd $HOME
/bin/rm -fr mpfr-3.0.1*
wget http://www.mpfr.org/mpfr-3.0.1/mpfr-3.0.1-rc1.tar.bz2
tar jxf mpfr-3.0.1-rc1.tar.bz2
cd mpfr-3.0.1-rc1
./configure --with-gmp=$HOME --prefix=$HOME --disable-shared
make
make check

# results with mpfr-3.0.1-rc1.tar.bz2
# gcc10 All 156 tests passed
# gcc11 All 156 tests passed
# gcc12 All 156 tests passed
# gcc13 All 156 tests passed
# gcc14 All 156 tests passed
# gcc15 All 156 tests passed
# gcc16 All 156 tests passed
# gcc17 All 156 tests passed
# gcc20 All 156 tests passed
# gcc33 All 156 tests passed
# gcc34 All 156 tests passed
# gcc35 Connection timed out
# gcc36 Connection timed out
# gcc37 All 156 tests passed
# gcc38 Connection timed out
# gcc40 Connection timed out
# gcc42 All 156 tests passed
# gcc43 Connection timed out
# gcc45 All 156 tests passed
# gcc46 All 156 tests passed
# gcc47 All 156 tests passed
# gcc50 Connection timed out
# gcc51 All 156 tests passed
# gcc52 All 156 tests passed
# gcc53 Connection timed out
# gcc54 All 156 tests passed
# gcc55 Connection timed out
# gcc56 Connection timed out
# gcc57 Connection timed out
# gcc60 All 156 tests passed
# gcc61 All 156 tests passed (GMP configured with ABI=1.0)
# gcc62 All 156 tests passed
# gcc63 All 156 tests passed
# gcc64 All 156 tests passed
# gcc70 All 156 tests passed
# gcc100 Connection closed by remote host
# gcc101 Connection closed by remote host
# gcc200 Connection timed out
# gcc201 Connection timed out
