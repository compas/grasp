#!/usr/bin/env bash

function generate-makefile {
	>&2 echo "Generating a Makefile"
	if ! [ -z ${EXE+x} ]; then
		_generate-makefile-binary
	elif ! [ -z ${LIB+x} ]; then
		_generate-makefile-library
	else
		>&2 echo "ERROR: neither EXE nor LIB specified"
		exit 1
	fi
}

function _generate-makefile-binary {
	cat <<-EOF
		BIN=\${GRASP}/bin
		EXE=${EXE}
	EOF

	if ! [ -z ${LIBRARIES+x} ]; then
		makelibs_string=$(for lib in ${LIBRARIES}; do echo -n "-l${lib} "; done)
		echo "LIBS=${makelibs_string}"
		LIBS="\$(LIBS)"
	else
		LIBS=""
	fi

	echo
	echo -n "OBJS="
	for file in ${FILES}; do
		if [[ "$file" =~ ^(.+)\.f90$ ]]; then
			echo " \\"
			echo -n "	${BASH_REMATCH[1]}.o"
		else
			>&2 echo "ERROR: Invalid file in FILES: $file"
			exit 1
		fi
	done
	echo; echo

	cat <<-EOF | sed 's/    /\t/'
		\$(BIN)/\$(EXE): \$(OBJS)
		    \$(FC) \$(FC_LD) -o \$@ ${LIBS} \$?

		%.o: %.f90
		    \$(FC) \$(FC_FLAGS) -c -o \$@ \$<

		clean:
		    -rm -f *.o *.mod
	EOF
}

function _generate-makefile-library {
	cat <<-EOF
		LIBA=\${GRASP}/lib/lib${LIB}.a
	EOF

	if ! [ -z ${LIBRARIES+x} ]; then
		makelibs_string=$(for lib in ${LIBRARIES}; do echo -n "-l${lib} "; done)
		echo "LIBS=${makelibs_string}"
		LIBS="\$(LIBS)"
	else
		LIBS=""
	fi

	echo
	echo -n "OBJS="
	for file in ${FILES}; do
		if [[ "$file" =~ ^(.+)\.f90$ ]]; then
			echo " \\"
			echo -n "	${BASH_REMATCH[1]}.o"
		else
			>&2 echo "ERROR: Invalid file in FILES: $file"
			exit 1
		fi
	done
	echo; echo

	cat <<-EOF | sed 's/    /\t/'
		\$(LIBA): \$(OBJS)
		    @echo "Installing \$@"
		    ar -curs \$@ \$?

		%.o: %.f90
		    \$(FC) \$(FC_FLAGS) -c -o \$@ \$<

		clean:
		    -rm -f \$(LIBA)
		    -rm -f *.o *.mod
	EOF
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
