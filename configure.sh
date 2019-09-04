#!/usr/bin/env bash
export GRASP="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

build_directory="build" # we default to build/

cmake_args=""
for arg in $@; do
	if [ "$arg" == "--debug" ]; then
		echo "Creating a DEBUG build"
		build_directory="build-debug"
		cmake_args="-DCMAKE_BUILD_TYPE=Debug"
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
To build GRASP run you need to cd into ${build_directory}/ and run make:

    cd ${build_directory}/
    make

Note that you probably want to also enable parallel builds by pass -j to make:

    make -jN

where N is the number of cores you have available.
EOF
