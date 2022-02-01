#!/usr/bin/env bash
export GRASP="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

build_directory="build" # we default to build/

cmake_args=""
for arg in $@; do
	if [ "$arg" == "--debug" ]; then
		echo "Creating a DEBUG build"
		build_directory="build-debug"
		cmake_args="-DCMAKE_BUILD_TYPE=Debug${cmake_args:+ $cmake_args}"
	fi
	if [ "$arg" == "--pic" ]; then
		echo "Force-enable position-independent code (e.g. -fPIC)"
		cmake_args="-DCMAKE_POSITION_INDEPENDENT_CODE=ON${cmake_args:+ $cmake_args}"
	fi
done

# We create an empty CMakeLists.user file, so that the user would not have to
# re-create it later.
touch "${GRASP}/CMakeLists.user" || exit

# Determine and check the build directory
build_abspath="${GRASP}/${build_directory}"
echo "Build directory: ${build_abspath}"
if [ -e "${build_abspath}" ]; then
	>&2 echo "ERROR: Build directory already exists."
	exit 1
fi

# We run the default setup for CMake's out-of-tree builds
#
#     mkdir build/
#     cd build/
#     cmake ..
#
mkdir "${build_abspath}" && cd "${build_abspath}" \
	&& cmake ${cmake_args} "${GRASP}" \
	|| exit

# Note: we need to use spaces, not tabs, to indent in the heredoc.
cat <<-EOF

 Build directory ${build_directory}/ created.

 To compile GRASP you need to cd into ${build_directory}/ and run make:

    cd ${build_directory}/
    make install

 which installs the GRASP binaries to the bin/ directory.

 Note that you also probably want to enable parallel build by passing -j to make:

    make -jN install

 where N is the number of cores you have available.

 Note: remove the build/ directory and rerun ./configure.sh to completely redo the configure step.

EOF
