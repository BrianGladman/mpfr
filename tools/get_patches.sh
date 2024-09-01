#!/bin/sh

set -e

if [ $# -ne 1 ]; then
  echo "Usage: $0 <PATCHES_file>" >&2
  exit 1
fi

patches=`cat "$1"`

cat <<EOF
/* mpfr_get_patches -- Patches that have been applied

Copyright 2007-2024 Free Software Foundation, Inc.
Contributed by the AriC and Caramba projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.
If not, see <https://www.gnu.org/licenses/>. */

#include "mpfr-impl.h"

const char *
mpfr_get_patches (void)
{
EOF

# $patches is written as a list of words separated by a space character.
echo '  return "'$patches'";'
echo '}'
