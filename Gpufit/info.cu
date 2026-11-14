#include "info.h"
#include <cuda_runtime.h>
#include <algorithm>

void Info::get_gpu_properties()
{
    cudaDeviceProp devProp;
    CUDA_CHECK_STATUS(cudaGetDeviceProperties(&devProp, 0));

    // Adding model functions to models.cuh can increase the number of registers
    // per thread used by the kernel cuda_calc_curve_values(). Exceeding the number
    // of available registers per thread block causes the error "too many resources
    // requested for launch". Reducing max_threads_ preserves the kernel cuda_calc_curve_values()
    // from reaching the maximum number of registers per thread block. If this error 
    // still occurs, comment out unused models in function calculate_model() in
    // file models.cuh.

    max_threads_ = std::min(devProp.maxThreadsPerBlock, 256);
    max_blocks_ = devProp.maxGridSize[0];
    warp_size_ = devProp.warpSize;

    std::size_t free_bytes;
    std::size_t total_bytes;
    CUDA_CHECK_STATUS(cudaMemGetInfo(&free_bytes, &total_bytes));
    available_gpu_memory_ = std::size_t(double(free_bytes) * 0.5);
    
    if (data_location_ == HOST && user_info_size_ > available_gpu_memory_)
    {
        throw std::runtime_error("0.5 free GPU memory budget cannot fit user info and one fit");
    }
}

int getDeviceCount()
{
	int deviceCount;
	CUDA_CHECK_STATUS(cudaGetDeviceCount(&deviceCount));
	return deviceCount;
}
