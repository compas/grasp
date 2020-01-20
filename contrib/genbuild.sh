#!/usr/bin/env bash
#
# Variables that affect the behaviour:
#
# - EXE & LIB: only one should be set, name of the target
# - FILES: List of .f90 source files, in the correct order (if there are between the files)
# - LIBRARIES: List of libraries that this library or binary depends on
#
# External dependencies:
# - LAPACK: Attach LAPACK and BLAS libraries
# - ISMPI: Attach MPI libraries
#
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
		EXE=\${GRASP}/bin/${EXE}
	EOF

	if ! [ -z ${ISMPI+x} ]; then
		FC="FC_MPI"
		FC_FLAGS="\$(FC_MPIFLAGS)"
		FC_LDFLAGS="\$(FC_MPILD)"
	else
		FC="FC"
		FC_FLAGS="\$(FC_FLAGS)"
		FC_LDFLAGS="\$(FC_LD)"
	fi

	if ! [ -z ${LIBRARIES+x} ]; then
		# TODO: trim space at the end of the string
		makelibs_string=$(for lib in ${LIBRARIES}; do echo -n " -l${lib}"; done)
		echo "LIBS=-L \${GRASP}/lib/${makelibs_string}"
		FC_LDFLAGS="${FC_LDFLAGS} \$(LIBS)"
		moddirs=$(for lib in ${LIBRARIES}; do echo -n " -I \${GRASP}/src/lib/$(libdir $lib)"; done)
		echo "FC_MODULES=${moddirs}"
		FC_FLAGS="${FC_FLAGS} \$(FC_MODULES)"
	fi

	if ! [ -z ${LAPACK} ]; then
		FC_LDFLAGS="${FC_LDFLAGS} \$(LAPACK_LIBS)"
	fi

	echo
	echo -n "OBJS="
	for file in $(echo "$FILES" | sed 's/#.*$//'); do
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
		\$(EXE): \$(OBJS)
		    \$($FC) -o \$@ \$? ${FC_LDFLAGS}

		%.o: %.f90
		    \$($FC) -c ${FC_FLAGS} -o \$@ \$<

		clean:
		    -rm -f \$(EXE)
		    -rm -f *.o *.mod
	EOF
}

function _generate-makefile-library {
	cat <<-EOF
		LIBA=\${GRASP}/lib/lib${LIB}.a
	EOF

	if ! [ -z ${ISMPI+x} ]; then
		FC="FC_MPI"
		FC_FLAGS="\$(FC_MPIFLAGS)"
	else
		FC="FC"
		FC_FLAGS="\$(FC_FLAGS)"
	fi

	if ! [ -z ${LIBRARIES+x} ]; then
		moddirs=$(for lib in ${LIBRARIES}; do echo -n " -I \${GRASP}/src/lib/$(libdir $lib)"; done)
		echo "FC_MODULES=${moddirs}"
		FC_FLAGS="${FC_FLAGS} \$(FC_MODULES)"
	fi

	echo
	echo -n "OBJS="
	for file in $(echo "$FILES" | sed 's/#.*$//'); do
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
		    \$($FC) -c ${FC_FLAGS} -o \$@ \$<

		clean:
		    -rm -f \$(LIBA)
		    -rm -f *.o *.mod
	EOF
}

function libdir {
	if [ "$1" = "mod" ]; then
		echo "libmod"
	elif [ "$1" = "9290" ]; then
		echo "lib9290"
	elif [ "$1" = "mpiu90" ]; then
		echo "mpi90"
	else
		echo "lib${1}"
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
	for file in $(echo "$FILES" | sed 's/#.*$//'); do
		echo "    ${file}"
	done
	echo ")"
	if ! [ -z ${LIB+x} ]; then
		echo "setup_fortran_modules($LIB)"
	fi
	if ! [ -z ${LIBRARIES+x} ]; then
		if ! [ -z ${EXE+x} ]; then
			echo "target_link_libraries_Fortran(${TARGET} PUBLIC ${LIBRARIES})"
		else
			echo "target_link_libraries_Fortran(${TARGET} PRIVATE ${LIBRARIES})"
		fi
	fi
	# Add LAPACK and BLAS libraries to libraries
	if [ -z ${EXE+x} ]; then
		if ! [ -z ${LAPACK} ]; then
			echo "target_link_libraries(${TARGET} PRIVATE \${BLAS_LIBRARIES} \${BLAS_LINKER_FLAGS})"
			echo "target_link_libraries(${TARGET} PRIVATE \${LAPACK_LIBRARIES} \${LAPACK_LINKER_FLAGS})"
		fi
		if ! [ -z ${ISMPI+x} ]; then
			cat <<-EOF
			target_include_directories(${TARGET} PRIVATE \${MPI_Fortran_INCLUDE_PATH})
			target_link_libraries(${TARGET} PRIVATE \${MPI_Fortran_LIBRARIES})
			set_target_properties(${TARGET} PROPERTIES
			  COMPILE_FLAGS "\${MPI_Fortran_COMPILE_FLAGS}"
			  LINK_FLAGS "\${MPI_Fortran_LINK_FLAGS}"
			)
			EOF
		fi
	fi
	if ! [ -z ${EXE+x} ]; then
		echo "install(TARGETS ${EXE} DESTINATION bin/)"
	else
		echo "install(TARGETS ${LIB} DESTINATION lib/)"
	fi
}

# The main script:
CMAKELISTSTXT=CMakeLists.txt
MAKEFILE=Makefile
for arg in $@; do
	if [ "$arg" = "--verify" ]; then
		VERIFY=true
		CMAKELISTSTXT=`tempfile`
		MAKEFILE=`tempfile`
		shift
	elif [[ $arg =~ ^- ]]; then
		>&2 echo "ERROR: Invalid argument $@"
		exit 1
	else
		break
	fi
done
if [ "$#" -ne 1 ]; then
	>&2 echo "ERROR: Must provide a single argument (target directory)"
	exit 1
fi
# Construct the path to the target directory, relative to $PWD if not an absolute path
if [[ "$1" =~ ^/(.+) ]]; then
	target=$1
else
	target="${PWD}/$1"
fi

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
	echo "$output"
	exit 2
else
	echo "BUILDCONF.sh ran successfully in $target"
	echo "Output:"
	echo "$output"

	if ! [ -z ${VERIFY+x} ]; then
		>&2 echo "INFO: Running in verification mode"
		>&2 echo " CMakeLists.txt = ${CMAKELISTSTXT}"
		>&2 echo " Makefile       = ${MAKEFILE}"

		diff ${CMAKELISTSTXT} "${target}/CMakeLists.txt" || {
			>&2 echo "ERROR: CMakeLists.txt differs"
			VERIFY=fail
		}
		diff ${MAKEFILE} "${target}/Makefile"

		rm -v "${CMAKELISTSTXT}" "${MAKEFILE}"

		if [ "$VERIFY" = "fail" ]; then
			exit 3
		fi
	fi
	exit 0
fi
