#!/bin/sh

# "autoreconf -f" will clobber our INSTALL file with a generic one if we
# don't move it out of the way.

# EXIT and signals that correspond to SIGHUP, SIGINT, SIGQUIT and SIGTERM.
signals="0 1 2 3 15"

cleanup()
{
  trap '' $signals
  if [ -f INSTALL.$$.tmp ]; then
    echo "$0: restoring the INSTALL file" >&2
    mv -f INSTALL.$$.tmp INSTALL
  fi
}

rm -f INSTALL.$$.tmp
trap cleanup $signals
mv -f INSTALL INSTALL.$$.tmp

autoreconf -v -f -i -W all
status=$?

rm -rf autom4te.cache

exit $status
