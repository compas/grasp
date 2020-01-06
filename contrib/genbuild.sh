#!/usr/bin/env bash

function generate-makefile {
	ofile="Makefile"
	echo "Generating ${ofile}"

	cat <<-EOF > ${ofile}
		BIN=\${GRASP}/bin
		EXE=${EXE}
	EOF

	if ! [ -z ${LIBRARIES+x} ]; then
		makelibs_string=$(for lib in ${LIBRARIES}; do echo -n "-l${lib} "; done)
		echo "LIBS=${makelibs_string}" >> $ofile
		LIBS="\$(LIBS)"
	else
		LIBS=""
	fi

	echo >> $ofile
	echo -n "OBJS=" >> $ofile
	for file in ${FILES}; do
		echo " \\" >> $ofile
		echo -n "    ${file}" >> $ofile
	done
	echo >> $ofile; echo >> $ofile

	cat <<-EOF >> ${ofile}
		\$(BIN)/\$(EXE): \$(OBJS)
		    \$(FC) \$(FC_LD) -o \$@ ${LIBS} \$?

		%.o: %.f90
		    \$(FC) \$(FC_FLAGS) -c -o \$@ \$<

		clean:
		    -rm -f *.o *.mod
	EOF
	# Replace 4 spaces with tabs
	sed -i 's/    /\t/' Makefile
}

function generate-cmakelists {
	ofile="CMakeLists.txt"
	echo "Generating ${ofile}"
	if ! [ -z ${EXE+x} ]; then
		echo "add_executable(${EXE}" > $ofile
	elif ! [ -z ${LIB+x} ]; then
		echo "add_library(${LIB} STATIC" > $ofile
	else
		>&2 echo "ERROR: neither EXE nor LIB specified"
		exit 1
	fi
	for file in ${FILES}; do
		echo "    ${file}" >> $ofile
	done
	echo ")" >> $ofile
	if ! [ -z ${EXE+x} ]; then
		if ! [ -z ${LIBRARIES+x} ]; then
			echo "target_link_libraries_Fortran(jj2lsj PUBLIC ${LIBRARIES})" >> $ofile
		fi
		echo "install(TARGETS ${EXE} DESTINATION bin/)" >> $ofile
	elif ! [ -z ${LIB+x} ]; then
		echo "setup_fortran_modules($LIB)" >> $ofile
	fi
}

# The main script:
if [ "$#" -ne 1 ]; then
	>&2 echo "ERROR: Must provide a single argument (target directory)"
	exit 1
fi
target="${PWD}/$1"

if ! [ -d "${target}" ]; then
	>&2 echo "ERROR: Invalid directory ${target}"
	exit 1
fi

if ! [ -f "${target}/BUILDCONF.sh" ]; then
	>&2 echo "ERROR: Missing BUILDCONF.sh file in ${target}"
	exit 1
fi

if ! output=$(cd ${target} || exit; source "${target}/BUILDCONF.sh"); then
	>&2 echo "ERROR: BUILDCONF.sh script failed for $target"
	echo "Output:"
	echo $output
	exit 2
else
	echo "BUILDCONF.sh ran successfully in $target"
	echo "Output:"
	echo $output
	exit 0
fi
