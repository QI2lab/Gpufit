rem create build files for Gpufit on Windows for Microsoft Visual Studio (MSVS) using CMake
rem ensure that "Visual C++ tools for CMake" are installed when installing MSVS, or add later under the menu "Tools -> Get Tools and Features"
rem assumption: this script is in the top level Gpufit source code folder, and the build directory will be created one level above
rem construct directory in the same folder as this script with timestamp
set gpufitdir=%0\..
set bdir=%0\..\..\%date:~10,4%_%date:~4,2%_%date:~7,2%_%time:~0,2%;%time:~3,2%;%time:~6,2%_Gpufit_build_64
rem remote spaces. This is useful because the environmental variable time does not zero bad the hour
set bdir=%bdir: =%
echo %bdir%
mkdir "%bdir%"
chdir "%bdir%"
rem command recommended here: https://gpufit.readthedocs.io/en/latest/installation.html
rem cmake -G "Visual Studio 14 2015 Win64" C:\Users\ptbrown2\Documents\Gpufit\Gpufit
rem  cmake -G "Visual Studio 14 2015 Win64" "%gpufitdir%"
rem if using version of nvcc/CUDA which is not compatible with MSVS 2013, then install MSVS 2019 and run
cmake -G "Visual Studio 16 2019" "%gpufitdir%"
