#!/bin/bash
export DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# cd into to root directory of the repository
cd "${DIR}/.." || exit

# Run Doxygen, which should create a doc/html directory
doxygen || {
	>&2 echo "ERROR: Doxygen failed with $?"
	exit 1
}
if ! [ -d "doc/html" ]; then
	>&2 echo "ERROR: doc/html missing -- documentation not built properly."
	ls -Alh doc/ # for debugging
	exit 1
fi
echo "Documentation built successfully."
