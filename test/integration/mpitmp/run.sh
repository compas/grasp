#!/bin/bash
#
# ./run.sh [--preserve-tmp] casename
#
# Create a temporary directory in the current working directory and run the rnucleus tests.
#
export SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Check for mandatory environment variables
if [ -z ${GRASP_BINDIR+x} ]; then
	>&2 echo "ERROR: \$GRASP_BINDIR variable is unset."
	exit 127
fi
GRASP_BINDIR=`realpath "${GRASP_BINDIR}"`
if ! [ -d "${GRASP_BINDIR}" ]; then
	>&2 echo "ERROR: \$GRASP_BINDIR is not a diretory."
	>&2 echo "  GRASP_BINDIR=$GRASP_BINDIR"
	exit 127
fi

# If the user passes "--preserve-tmp" as the first argument, we'll note that down for later.
if [ "$1" == "--preserve-tmp" ]; then
	preserve_tmp=true
	shift
fi

# Get the case name. The user should pass only one argument.
if [ "$#" -ne 1 ]; then
	>&2 echo "ERROR: Invalid number of arguments passed ($#)"
	exit 127
fi
CASE=$1
if [ "$CASE" != "serial" ] && [ "$CASE" != "mpi" ] && [ "$CASE" != "mpi-longpath" ]; then
	>&2 echo "ERROR: Invalid configuration passed ($CASE)"
	>&2 echo "  expected one of: serial, mpi, mpi-longpath"
	exit 125
fi

# Determine the paths to the necessary binaries
function checkbinary {
	varname=$1
	path=$2
	if ! [ -f "${path}" ]; then
		>&2 echo "ERROR: Unable to find binary $varname"
		>&2 echo "  at ${path}"
		exit 127
	fi
	>&2 echo "INFO: $varname=${path}"
}
RNUCLEUS="${GRASP_BINDIR}/rnucleus"; checkbinary RNUCLEUS $RNUCLEUS
RCSFGENERATE="${GRASP_BINDIR}/rcsfgenerate"; checkbinary RCSFGENERATE $RCSFGENERATE
RWFNESTIMATE="${GRASP_BINDIR}/rwfnestimate"; checkbinary RWFNESTIMATE $RWFNESTIMATE
if [[ $CASE =~ ^mpi ]]; then
	checkbinary RANGULAR "${GRASP_BINDIR}/rangular_mpi"
	checkbinary RMCDHF "${GRASP_BINDIR}/rmcdhf_mpi"
	RANGULAR="mpirun -n 4 ${GRASP_BINDIR}/rangular_mpi"
	RMCDHF="mpirun -n 4 ${GRASP_BINDIR}/rmcdhf_mpi"
else
	RANGULAR="${GRASP_BINDIR}/rangular"; checkbinary RANGULAR $RANGULAR
	RMCDHF="${GRASP_BINDIR}/rmcdhf"; checkbinary RMCDHF $RMCDHF
fi

# Create a temporary directory to run GRASP in:
TMP=`mktemp -d grasp-test-mpitmp.XXXXXXXXX` || exit 120
TMP=`realpath "${TMP}"`
if ! [ -d "${TMP}" ]; then
	>&2 echo "ERROR: Temporary directory as not created."
	>&2 echo "  TMP=$TMP"
	exit 121
fi
# This will be called any time we exit, independent of whether it's an early exit due to a
# failure, or the final exit when tests pass.
function clean_up_tmp_directory {
	# Keep the temporary directory around if the user passed --preserve-tmp
	if [ "$preserve_tmp" == "true" ]; then
		>&2 echo "INFO: Keeping temporary directory ${TMP}"
	else
		>&2 echo "INFO: Removing temporary directory ${TMP}"
		rm -vR "${TMP}"
	fi
}
trap clean_up_tmp_directory EXIT
>&2 echo "INFO: switching to temporary directory: ${TMP}"
cd "${TMP}" || exit 122

# Function to test existence of a generated file:
function test_file_exists {
	if ! [ -f "$1" ]; then
		>&2 echo "ERROR: failed to generate file $1"
		exit 50
	fi
}

# Run rnucleus to generate a simple isodata file
${RNUCLEUS} <<-EOF
	92
	238
	n
	238.02891
	0
	0
	0
EOF
exitcode=$?
if ! [ $exitcode -eq 0 ]; then
	>&2 echo "ERROR: rnucleus failed with $exitcode"
	exit 1
fi
test_file_exists "isodata"

# Run rcsfgenerate to generate a simple CSL
${RCSFGENERATE} <<-EOF
	*
	0
	1s(2,*)2s(2,*)

	2s,2p
	0,2
	2
	n
EOF
exitcode=$?
if ! [ $exitcode -eq 0 ]; then
	>&2 echo "ERROR: rcsfgenerate failed with $exitcode"
	exit 1
fi
test_file_exists "rcsf.out"
mv rcsf.out rcsf.inp || exit 2
test_file_exists "rcsf.inp"

# Run rwfnestimate to generate basic orbitals
${RWFNESTIMATE} <<-EOF
	y
	2
	*
EOF
exitcode=$?
if ! [ $exitcode -eq 0 ]; then
	>&2 echo "ERROR: rwfnestimate failed with $exitcode"
	exit 1
fi
test_file_exists "rwfn.inp"

# Set up MPI_TMP on MPI cases
if [[ $CASE =~ ^mpi ]]; then
	if [ "$CASE" == "mpi-longpath" ]; then
		export MPI_TMP="$TMP/mpitmp"
	else
		export MPI_TMP=`mktemp -d`
		function clean_up_mpitmp {
			rm -Rv ${MPI_TMP}
		}
		trap clean_up_mpitmp EXIT
	fi
	nchars=`echo -n $MPI_TMP | wc -c`
	echo "$MPI_TMP ($nchars characters)"
fi

# Run rangular
echo "Running: ${RANGULAR}"
${RANGULAR} <<-EOF
	y
EOF
exitcode=$?
if ! [ $exitcode -eq 0 ]; then
	>&2 echo "ERROR: rangular failed with $exitcode"
	exit 1
fi

# Run rmcdhf
echo "Running: ${RMCDHF}"
${RMCDHF} <<-EOF
	y
	1
	1
	5
	*
	*
	20
EOF
exitcode=$?
if ! [ $exitcode -eq 0 ]; then
	>&2 echo "ERROR: rmcdhf failed with $exitcode"
	exit 1
fi
test_file_exists "rmcdhf.sum"
test_file_exists "rmix.out"

echo "INFO: Final directory contents:"
ls -Alh

# If we got this far, everything is a-ok
>&2 echo "TESTS SUCCEEDED"
exit 0
