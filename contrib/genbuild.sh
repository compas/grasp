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
		moddirs=$(for lib in ${LIBRARIES}; do echo -n "-I \${GRASP}/src/lib/$(libdir $lib) "; done)
		echo "FC_MODULES=${moddirs}"
		FC_MODS=" \$(FC_MODULES)" # Note: space is significant!
	else
		FC_MODS=""
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
		    \$(FC) \$(FC_FLAGS)${FC_MODS} -c -o \$@ \$<

		clean:
		    -rm -f \$(LIBA)
		    -rm -f *.o *.mod
	EOF
}

function libdir {
	if [ "$1" = "mod" ]; then
		echo "libmod"
	elif [ "$1" = "mpi" ]; then
		echo "mpi90"
	else
		echo "lib${1}90"
	fi
}

function generate-cmakelists {
	>&2 echo "Generating CMakeLists.txt"
	if ! [ -z ${EXE+x} ]; then
		TARGET=$EXE
		echo "add_executable(${EXE}"
	elif ! [ -z ${LIB+x} ]; then
		TARGET=$LIB
		echo "add_library(${LIB} STATIC"
	else
		>&2 echo "ERROR: neither EXE nor LIB specified"
		exit 1
	fi
	for file in ${FILES}; do
		echo "    ${file}"
	done
	echo ")"
	if ! [ -z ${LIB+x} ]; then
		echo "setup_fortran_modules($LIB)"
	fi
	if ! [ -z ${LIBRARIES+x} ]; then
		echo "target_link_libraries_Fortran(${TARGET} PRIVATE ${LIBRARIES})"
	fi
	# Add LAPACK and BLAS libraries
	if ! [ -z ${LAPACK} ]; then
		echo "target_link_libraries(${TARGET} PRIVATE \${BLAS_LIBRARIES} \${BLAS_LINKER_FLAGS})"
		echo "target_link_libraries(${TARGET} PRIVATE \${LAPACK_LIBRARIES} \${LAPACK_LINKER_FLAGS})"
	fi
	if ! [ -z ${EXE+x} ]; then
		echo "install(TARGETS ${EXE} DESTINATION bin/)"
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

if ! output=$(cd ${target} || exit; source "${target}/BUILDCONF.sh" 2>&1); then
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
