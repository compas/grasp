#!/usr/bin/env bash
export DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

BINARIES="
hf
jj2lsj
jjgen
lscomp.pl
rangular
rangular_mpi
rasfsplit
rbiotransform rbiotransform_mpi
rci rci_mpi
rcsfblock
rcsfgenerate
rcsfinteract
rcsfmr
rcsfsplit
rcsfzerofirst
rhfs
rhfs_lsj
rlevels
rlevelseV
rmcdhf rmcdhf_mpi
rmixaccumulate
rmixextract
rnucleus
rsave
rseqenergy
rseqhfs
rseqtrans
rsms
rtabhfs
rtablevels
rtabtrans1 rtabtrans2 rtabtransE1
rtransition rtransition_mpi
rwfnestimate
rwfnmchfmcdf
rwfnplot
rwfnrelabel
rwfnrotate
wfnplot
"

LIBRARIES="
lib9290.a
libdvd90.a
libmcp90.a
libmod.a
libmpiu90.a
librang90.a
"

success=true

BIN="${DIR}/../bin"
for p in ${BINARIES}; do
	if ! [ -f "${BIN}/$p" ]; then
		>&2 echo "ERROR: binary ${p} missing from bin/"
		success=false
	fi
done

LIB="${DIR}/../lib"
for lib in ${LIBRARIES}; do
	if ! [ -f "${LIB}/$lib" ]; then
		>&2 echo "ERROR: library ${lib} missing from lib/"
		success=false
	fi
done

if [ "$success" = "false" ]; then
	>&2 echo "FAIL: Verification failed, check the logs."
	>&2 echo "INFO: ls -Alh bin/"
	ls -Alh "${BIN}"
	>&2 echo "INFO: ls -Alh lib/"
	ls -Alh "${LIB}"
	exit 1
else
	>&2 echo "SUCCESS: Found all the binaries in all the right places."
	exit 0
fi
