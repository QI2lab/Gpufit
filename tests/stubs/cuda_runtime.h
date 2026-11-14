// Only the CUDA queries used by Info; no device memory is allocated by these tests.
#pragma once
#include <cstddef>

typedef int cudaError_t;
struct cudaDeviceProp
{
    int maxThreadsPerBlock;
    int maxGridSize[3];
    int warpSize;
};

cudaError_t cudaGetDeviceProperties(cudaDeviceProp *, int);
cudaError_t cudaMemGetInfo(std::size_t *, std::size_t *);
cudaError_t cudaGetDeviceCount(int *);
char const * cudaGetErrorString(cudaError_t);
