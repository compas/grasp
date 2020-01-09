#!/usr/bin/env bash
export DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
GENBUILD="${DIR}/genbuild.sh"
SRCDIR="${DIR}/../src"
if ! [ -f "${GENBUILD}" ]; then >&2 echo "ERROR: Unable to find genbuild.sh at ${GENBUILD}"; exit 1; fi
if ! [ -d "${SRCDIR}" ]; then >&2 echo "ERROR: Not a directory: ${SRCDIR}"; exit 1; fi
if ! [ -d "${SRCDIR}/appl" ]; then >&2 echo "ERROR: Not a directory: ${SRCDIR}/appl"; exit 1; fi
if ! [ -d "${SRCDIR}/lib" ]; then >&2 echo "ERROR: Not a directory: ${SRCDIR}/lib"; exit 1; fi
if ! [ -d "${SRCDIR}/tool" ]; then >&2 echo "ERROR: Not a directory: ${SRCDIR}/tool"; exit 1; fi

directories="
	$(find ${SRCDIR}/appl/ -mindepth 1 -maxdepth 1 -type d)
	$(find ${SRCDIR}/lib/ -mindepth 1 -maxdepth 1 -type d)
	${SRCDIR}/tool
"
for d in ${directories}; do
	if ! [ -f "${d}/BUILDCONF.sh" ]; then
		>&2 echo "> WARNING: BUILDCONF.sh missing in $(basename $d)"
		>&2 echo ">  in $d"
		continue
	fi
	echo "> Calling genbuild.sh for $(basename $d)"
	$GENBUILD "$d" || {
		>&2 echo "> WARNING: genbuild.sh failed with $? for $(basename $d)"
		>&2 echo ">  in $d"
	}
done
