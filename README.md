# Notes on building the qi2lab fork of Gpufit:
First clone the repository
```
git clone https://github.com/QI2lab/Gpufit.git
cd Gpufit
```

## building on Windows
The build instructions for the Gpufit source suggest using Visual Studio 2013. However, this version is not compatible
with certain all versions of CUDA/nvcc. 

For example, with CUDA 11.2 the build files can be created with CMake using

`cmake -G "Visual Studio 16 2019" C:\Sources\Gpufit`

See my example batch files [for Windows](run_cmake_windows.bat). This script will generate Visual Studio project files. 

Then, navigate to the build directory and open ALL_BUILD.vcxproj with Visual Studio. Go to `build -> build solution` and build the project.

## building on Linux
The Gpufit source code can be built on Linux using cmake and make. See my [example shell script](build_gpufit_linux.sh)

# Installing pygpufit
After compiling, the python package is located in either `<build directory>\pyGpufit` or `<build directory>\Debug\pyGpufit` depending on your build options
Navigate to this folder, activate your conda or virtual environment, and install with `pip install .`

# Gpufit

Levenberg Marquardt curve fitting in CUDA.

Homepage: [github.com/gpufit/Gpufit](https://github.com/gpufit/Gpufit)

The manuscript describing Gpufit is now published in [Scientific Reports](https://www.nature.com/articles/s41598-017-15313-9).

## Quick start instructions

To verify that Gpufit is working correctly on the host computer, go to the folder gpufit_performance_test of the binary package and run Gpufit_Cpufit_Performance_Comparison.exe.  Further details of the test executable can be found in the documentation package.

## Binary distribution

The latest Gpufit binary release, supporting Windows 32-bit and 64-bit machines, can be found on the [release page](https://github.com/gpufit/Gpufit/releases).

## Documentation

[![Documentation Status](https://readthedocs.org/projects/gpufit/badge/?version=latest)](http://gpufit.readthedocs.io/en/latest/?badge=latest)

Documentation for the Gpufit library may be found online ([latest documentation](http://gpufit.readthedocs.io/en/latest/?badge=latest)), and also
as a PDF file in the binary distribution of Gpufit.

## Building Gpufit from source code

Instructions for building Gpufit are found in the documentation: [Building from source code](https://github.com/gpufit/Gpufit/blob/master/docs/installation.rst).

## Using the Gpufit binary distribution

Instructions for using the binary distribution may be found in the documentation.  The binary package contains:

- The Gpufit SDK, which consists of the 32-bit and 64-bit DLL files, and 
  the Gpufit header file which contains the function definitions.  The Gpufit
  SDK is intended to be used when calling Gpufit from an external application
  written in e.g. C code.
- Gpufit Performance test: A simple console application comparing the execution speed of curve fitting on the GPU and CPU.  This program also serves as a test to ensure the correct functioning of Gpufit.
- Matlab 32 bit and 64 bit bindings, with Matlab examples.
- Python version 2.x and version 3.x bindings (compiled as wheel files) and
  Python examples.
- Java binding, with Java examples.
- The Gpufit manual in PDF format

## Authors

Gpufit was created by Mark Bates, Adrian Przybylski, Björn Thiel, and Jan Keller-Findeisen at the Max Planck Institute for Biophysical Chemistry, in Göttingen, Germany.

## How to cite Gpufit

If you use Gpufit in your research, please cite our publication describing the software.  A paper describing the software was published in Scientific Reports.  The open-access manuscript is available from the Scientific Reports website, [here](https://www.nature.com/articles/s41598-017-15313-9).

  *  Gpufit: An open-source toolkit for GPU-accelerated curve fitting  
     Adrian Przybylski, Björn Thiel, Jan Keller-Findeisen, Bernd Stock, and Mark Bates  
     Scientific Reports, vol. 7, 15722 (2017); doi: https://doi.org/10.1038/s41598-017-15313-9 

## License

MIT License

Copyright (c) 2017 Mark Bates, Adrian Przybylski, Björn Thiel, and Jan Keller-Findeisen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
