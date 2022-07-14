#!/bin/sh
set -e # errexit
aclocal
autoheader
autoconf
if ! test -f ltmain.sh ; then
	libtoolize -c ;
fi
automake -c -Woverride --add-missing

