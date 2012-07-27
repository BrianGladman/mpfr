#!/bin/sh

# "autoreconf -f" will clobber our INSTALL file with a generic one if we
# don't move it out of the way

mv -f INSTALL INSTALL.$$.tmp

autoreconf -v -f -i -W all

rm -f INSTALL
mv -f INSTALL.$$.tmp INSTALL

rm -rf autom4te.cache
