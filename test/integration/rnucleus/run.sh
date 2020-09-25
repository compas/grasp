#!/bin/bash
#
# ./run.sh [--preserve-tmp] casename
#
# Create a temporary directory in the current working directory and run the rnucleus tests.
#
export SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Check for mandatory environement variables
if [ -z ${GRASP_BINDIR+x} ]; then
	>&2 echo "ERROR: \$GRASP_BINDIR variable is unset."
	exit 127
fi
GRASP_BINDIR=`realpath "${GRASP_BINDIR}"`
# Determine the path to the rnucleus binary
RNUCLEUS="${GRASP_BINDIR}/rnucleus"
if ! [ -f "${RNUCLEUS}" ]; then
	>&2 echo "ERROR: Unable to find rnucleus binary"
	>&2 echo "  at: ${RNUCLEUS}"
	exit 127
fi
>&2 echo "RNUCLEUS=${RNUCLEUS}"

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
if ! [[ "$CASE" =~ ^[A-Za-z0-9_-]+$ ]]; then
	>&2 echo "ERROR: Case name ($CASE) must match '^[A-Za-z0-9_-]+d$'"
	exit 125
fi
STDIN="${SCRIPTDIR}/${CASE}.stdin"
STDOUT="${CASE}.stdout"
STDERR="${CASE}.stderr"
if ! [ -f "${STDIN}" ]; then
	>&2 echo "ERROR: Missing input file for case '$CASE'"
	>&2 echo "  at: ${STDIN}"
	exit 126
fi

# Create a temporary directory to run GRASP in:
TMP=`mktemp -d grasp-test-rnucleus.XXXXXXXXX` || exit 120
TMP=`realpath "${TMP}"`
if ! [ -d "${TMP}" ]; then
	>&2 echo "ERROR: Temporary directory as not created."
	>&2 echo "  TMP=$TMP"
	exit 121
fi
>&2 echo "INFO: switching to temporary directory: ${TMP}"
cd "${TMP}" || exit 122
# Run rnucleus. If any error conditions are detected, we'll set `success` to false
success=true
>&2 echo "INFO: running rnucleus for case=${CASE}"
>&2 echo " PWD=${PWD}"
>&2 echo "--------------------------------------------------------------------------------"
>&2 echo "STDIN> ${STDIN}:"
>&2 cat ${STDIN}
>&2 echo "--------------------------------------------------------------------------------"
${RNUCLEUS} <${STDIN} 1>${STDOUT} 2>${STDERR}
exitcode=$?
>&2 echo "INFO: rnucleus exited with ${exitcode}"
>&2 echo "--------------------------------------------------------------------------------"
>&2 echo "STDOUT> ${STDOUT}:"
>&2 cat ${STDOUT}
>&2 echo "--------------------------------------------------------------------------------"
>&2 echo "STDERR> ${STDERR}:"
>&2 cat ${STDERR}
>&2 echo "--------------------------------------------------------------------------------"
# Check for the exit code:
if ! [ "${exitcode}" -eq 0 ]; then
	>&2 echo "TEST FAILED: rnucleus exit code ($exitcode) not zero"
	success=false
fi
# Check for the output isodata file:
if ! [ -f "isodata" ]; then
	>&2 echo "TEST FAILED: missing output isodata file in $PWD"
	success=false
else
	>&2 echo "--------------------------------------------------------------------------------"
	>&2 echo "isodata file> $PWD/isodata"
	>&2 cat $PWD/isodata
	>&2 echo "--------------------------------------------------------------------------------"
fi

# Keep the temporary directory around if the user passed --preserve-tmp
if [ "$preserve_tmp" == "true" ]; then
	>&2 echo "INFO: Keeping temporary directory ${TMP}"
else
	>&2 echo "INFO: Removing temporary directory ${TMP}"
	rm -vR "${TMP}"
fi

# Exit with code 1 if we detected any error conditions
if [ "$success" == "false" ]; then
	>&2 echo "TESTS FAILED, check the log for errors."
	exit 1
else
	>&2 echo "TESTS SUCCEEDED"
	exit 0
fi
