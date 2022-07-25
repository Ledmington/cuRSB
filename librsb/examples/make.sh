#!/bin/bash
# Example script to build the librsb example programs.
# Uses librsb-config in the $PATH for build flags.
# Environment-provided $LIBRSB_CONFIG override that.
set -e
set -x

srcdir=${srcdir:-`pwd`}
builddir=${builddir:-`pwd`}
prefix="/usr/local"
exec_prefix="${prefix}"
bindir="${exec_prefix}/bin"
if test -z "${LIBRSB_CONFIG}"; then
	export PATH="${bindir}:$PATH"
fi
LIBRSB_CONFIG=${LIBRSB_CONFIG:-librsb-config}
PKG_CONFIG=pkg-config
WANT_PKGCONFIG="no"
WANT_CLEANUP=${WANT_CLEANUP:-false}

if test x"yes" == x"yes" ; then
CXX="`${LIBRSB_CONFIG} --cxx`"
if test x"rsblib" != x"" -a x"${CXX}" != x"" ; then
for s in ${srcdir}/*.cpp
do
	p=${builddir}/`basename ${s/.cpp/}`
	rm -f $p
	CXXFLAGS=`${LIBRSB_CONFIG} --cxxflags --I_opts`
	LDFLAGS=`${LIBRSB_CONFIG} --ldflags --extra_libs`
	LINK=`${LIBRSB_CONFIG} --link`
	o="${p}.o"
	ccmd="$CXX $CXXFLAGS -c $s -o $o"
	lcmd="$LINK $o $LDFLAGS -o $p"
	echo "$ccmd && $lcmd"
	( $ccmd && $lcmd )
	${WANT_CLEANUP} && rm -f "$p"
	# one may use a single command, but that's error-prone (may miss libraries):
	#cmd="$CXX $CXXFLAGS $s $LDFLAGS -o $p"
	#echo $cmd
	#$cmd
done
fi
fi

if test x"yes" == x"yes" ; then
for s in ${srcdir}/*.c
do
	p=`basename ${s/.c/}`
	if test $p == hello-spblas -a x"yes" != x"yes" ; then continue; fi
	if test $p ==    io-spblas -a x"yes" != x"yes" ; then continue; fi
	rm -f $p 
	CFLAGS=`${LIBRSB_CONFIG} --I_opts --cppflags`
       	LDFLAGS=`${LIBRSB_CONFIG} --ldflags --extra_libs`
	CC=`${LIBRSB_CONFIG} --cc`
	LINK=`${LIBRSB_CONFIG} --link`
	o="${p}.o"
	ccmd="$CC $CFLAGS -c $s -o $o"
	lcmd="$LINK $o $LDFLAGS -o $p"
	echo "$ccmd && $lcmd"
	( $ccmd && $lcmd )
	${WANT_CLEANUP} && rm -f "$p"
	# one may use a single command, but that's error-prone (may miss libraries):
	#cmd="$CC $CFLAGS $s $LDFLAGS -o $p"
	#echo $cmd
	#$cmd
	if test x"${WANT_PKGCONFIG}" != x"no" ; then
		CFLAGS=`${PKG_CONFIG} --cflags librsb`
		LIBS=`${PKG_CONFIG} --libs --static librsb`
		ccmd="$CC $CFLAGS -c $s -o $o"
		lcmd="$LINK $o $LIBS -o $p"
		${WANT_CLEANUP} && rm -f "$p"
		echo "$ccmd && $lcmd"
		( $ccmd && $lcmd )
	fi
done
fi

if test x"no" == x"yes" ; then
if test x"yes" = x"yes" ; then
FP=${srcdir}/fortran_rsb_fi.F90
if test x"yes" = x"yes" ; then
	FP+=\ ${srcdir}/fortran.F90
fi
# activated if you have built the Fortran modules and installed them in the right path.
for s in $FP
do
	p=`basename ${s/.F90/}`
	rm -f $p 
	FCFLAGS=`${LIBRSB_CONFIG} --I_opts --fcflags`
       	LDFLAGS=`${LIBRSB_CONFIG} --ldflags --extra_libs`
	FC=`${LIBRSB_CONFIG} --fc`
	LINK=`${LIBRSB_CONFIG} --link`
	o="${p}.o"
	FCLIBS=`${LIBRSB_CONFIG} --fclibs`
	ccmd="$FC $FCFLAGS -c $s -o $o"
	lcmd="$LINK $o $LDFLAGS $FCLIBS -o $p"
	echo "$ccmd && $lcmd"
	( $ccmd && $lcmd )
	${WANT_CLEANUP} && rm -f "$p"
done
fi
fi

echo " [*] done building examples!"
