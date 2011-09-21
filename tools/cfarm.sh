# this script helps to check MPFR is running on the GCC compile farm
# 1) update the GMP version if needed
# 2) update the MPFR release candidate
# 3) ssh gcc100 < cfarm.sh
/bin/rm -fr gmp-5.0.2*
wget http://mirror.ibcp.fr/pub/gnu/gmp/gmp-5.0.2.tar.bz2
tar jxf gmp-5.0.2.tar.bz2
cd gmp-5.0.2
./configure --prefix=$HOME --disable-shared
make
make install
cd $HOME
/bin/rm -fr mpfr-3.1.0*
wget http://www.mpfr.org/mpfr-3.1.0/mpfr-3.1.0-rc2.tar.bz2
tar jxf mpfr-3.1.0-rc2.tar.bz2
cd mpfr-3.1.0-rc2
./configure --with-gmp=$HOME --prefix=$HOME --disable-shared
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
