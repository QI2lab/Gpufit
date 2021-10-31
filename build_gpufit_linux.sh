#!/bin/sh
# build gpufit in a new directory prefixed with the current date/time. Assume that this script is run from the top directory of the gpufit source code. The code will be built in a folder one level above this

# path to cmake command
cmake_path="/home/ptbrown/cmake-3.20.1-linux-x86_64/bin/cmake"

# directory where gpu
build_dir="../"$(date  +"%Y_%m_%d_%H_%M_%S")_gpufit_build

mkdir $build_dir
cd $build_dir
$cmake_path -DCMAKE_BUILD_TYPE=RELEASE ../Gpufit
make
