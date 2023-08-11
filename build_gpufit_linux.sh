#!/bin/sh
# build gpufit in a new directory prefixed with the current date/time. Assume that this script is run from the top directory of the gpufit source code. The code will be built in a folder one level above this

script_dir="$(dirname -- "$(readlink -f "${BASH_SOURCE}")")"

# #######################
# define location of cmake
# #######################
# input path to cmake command manually if desired
#cmake_path="/home/ptbrown/cmake-3.20.1-linux-x86_64/bin/cmake"

# otherwise use system cmake
cmake_path="cmake"

# #######################
# set build directory
# #######################
# directory where gpufit will be built
build_dir_rel="$script_dir""/../"$(date  +"%Y_%m_%d_%H_%M_%S")_gpufit_build
# resolve full path for convenience
build_dir="$(readlink -f "${build_dir_rel}")"

# build directory and move there
mkdir $build_dir
cd $build_dir

echo "building GPUfit from ${script_dir} to directory ${build_dir}"

# #######################
# build
# #######################

# generate build files
$cmake_path -DCMAKE_BUILD_TYPE=RELEASE ../Gpufit "$script_dir"

# build GPUfit
make

