#!/bin/bash
#
# Copyright (C) 2008-2021 Michele Martone
# 
# This file is part of librsb.
# 
# librsb is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# 
# librsb is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with librsb; see the file COPYING.
# If not, see <http://www.gnu.org/licenses/>.

# This script is intended for the librsb developer usage.

# This script gets a bunch of useful test matrices.

# our default matrix set is :
#bayer02.mtx  coater2.mtx  crystk03.mtx  ex11.mtx  lhr10.mtx  memplus.mtx  orani678.mtx  raefsky4.mtx  wang4.mtx

#http://sparse.tamu.edu/MM/Pothen/commanche_dual.tar.gz
#http://sparse.tamu.edu/MM/Simon/venkat01.tar.gz
#http://sparse.tamu.edu/MM/Simon/venkat25.tar.gz
#http://sparse.tamu.edu/MM/Simon/venkat50.tar.gz

MATRICES="\
http://sparse.tamu.edu/MM/Grund/bayer02.tar.gz		\
http://sparse.tamu.edu/MM/Brethour/coater2.tar.gz	\
http://sparse.tamu.edu/MM/Boeing/crystk03.tar.gz	\
http://sparse.tamu.edu/MM/FIDAP/ex11.tar.gz		\
http://sparse.tamu.edu/MM/Mallya/lhr10.tar.gz		\
http://sparse.tamu.edu/MM/Hamm/memplus.tar.gz		\
http://sparse.tamu.edu/MM/HB/orani678.tar.gz		\
http://sparse.tamu.edu/MM/Simon/raefsky4.tar.gz	\
http://sparse.tamu.edu/MM/Simon/raefsky3.tar.gz	\
http://sparse.tamu.edu/MM/Wang/wang4.tar.gz"

MATRICES_CSB="\
http://sparse.tamu.edu/MM/Sandia/ASIC_320k.tar.gz	\
http://sparse.tamu.edu/MM/FEMLAB/sme3Dc.tar.gz		\
http://sparse.tamu.edu/MM/Wissgott/parabolic_fem.tar.gz\
http://sparse.tamu.edu/MM/Mittelmann/cont11_l.tar.gz	\
http://sparse.tamu.edu/MM/Rucci/Rucci1.tar.gz		\
http://sparse.tamu.edu/MM/Norris/torso1.tar.gz		\
http://sparse.tamu.edu/MM/Zaoui/kkt_power.tar.gz	\
http://sparse.tamu.edu/MM/Rajat/rajat31.tar.gz		\
http://sparse.tamu.edu/MM/GHS_psdef/ldoor.tar.gz	\
http://sparse.tamu.edu/MM/Oberwolfach/bone010.tar.gz"


[[ -d "$1" ]] && { cd "$1" || exit -1 ; }

for m in $MATRICES
do
	mbn=`basename $m`
	mn=${mbn//.tar.gz/}
	mfn=$mn.mtx

#	file based
#	[ -f $mbn ] || wget $m
#	tar xzf $mbn $mn/$mfn -O > $mfn
	
	# pipe based
	[ -f $mfn ] || wget $m -O - | tar xzf - $mn/$mfn -O > $mfn || exit -1

done
