# Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 Free Software Foundation, Inc.
# This Makefile.am is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.


AUTOMAKE_OPTIONS = gnu
ACLOCAL_AMFLAGS = -I m4

SUBDIRS = doc src tests tune

nobase_dist_doc_DATA = AUTHORS BUGS COPYING COPYING.LESSER NEWS TODO \
  examples/ReadMe examples/divworst.c examples/rndo-add.c examples/sample.c \
  examples/version.c

EXTRA_DIST = PATCHES VERSION tools/ck-version-info tools/get_patches.sh

# In a "make dist", check consistency concerning -version-info
# and if the VERSION file doesn't end with "-dev", check that
# the -version-info value is up-to-date.
# Note: this is a heuristic, to detect some mistakes.
dist-hook:
	cd $(srcdir) && tools/ck-version-info

# The following needs to be removed once GNU Automake 1.11 support is added.
DIST_ARCHIVES = $(distdir).tar.gz $(distdir).tar.bz2 $(distdir).tar.xz \
	$(distdir).zip
dist dist-all: distdir
	tardir=$(distdir) && $(am__tar) | GZIP=$(GZIP_ENV) gzip -c >$(distdir).tar.gz
	tardir=$(distdir) && $(am__tar) | bzip2 -9 -c >$(distdir).tar.bz2
	tardir=$(distdir) && $(am__tar) | xz -c >$(distdir).tar.xz
	-rm -f $(distdir).zip
	zip -rq $(distdir).zip $(distdir)
	$(am__remove_distdir)
