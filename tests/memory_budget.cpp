// Exercise the production planner with synthetic CUDA readings, without allocating VRAM.
#include "Gpufit/info.cu"
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>

static std::size_t free_bytes;
static std::size_t total_bytes = std::numeric_limits<std::size_t>::max();
static int max_blocks = 1024;

cudaError_t cudaGetDeviceProperties(cudaDeviceProp * properties, int)
{
    properties->maxThreadsPerBlock = 1024;
    properties->maxGridSize[0] = max_blocks;
    properties->warpSize = 32;
    return 0;
}
cudaError_t cudaMemGetInfo(std::size_t * available, std::size_t * total)
{
    *available = free_bytes;
    *total = total_bytes;
    return 0;
}
cudaError_t cudaGetDeviceCount(int * count) { *count = 1; return 0; }
char const * cudaGetErrorString(cudaError_t) { return "mock CUDA error"; }

static void check(bool condition)
{
    if (!condition)
        throw std::runtime_error("memory budget check failed");
}

static Info make_info(DataLocation location, bool constraints, bool weights)
{
    Info info;
    info.n_parameters_ = 4;
    info.n_parameters_to_fit_ = 2;
    info.n_points_ = weights ? 1024 : 16;
    info.n_fits_ = 1;
    info.data_location_ = location;
    info.use_constraints_ = constraints;
    info.use_weights_ = weights;
    return info;
}

static std::size_t capacity(Info info, std::size_t budget)
{
    free_bytes = 2 * budget;
    try
    {
        info.configure();
        return info.max_chunk_size_;
    }
    catch (std::runtime_error const & error)
    {
        if (std::string(error.what()).find("0.5 free GPU memory budget") == std::string::npos)
            throw;
        return 0;
    }
}

// Find the production planner's boundary rather than duplicating its per-fit formula.
static std::size_t minimum_budget(Info info)
{
    std::size_t lower = 0, upper = 1;
    while (capacity(info, upper) < info.n_fits_)
        upper *= 2;
    while (upper - lower > 1)
    {
        std::size_t const middle = lower + (upper - lower) / 2;
        if (capacity(info, middle) < info.n_fits_)
            lower = middle;
        else
            upper = middle;
    }
    return upper;
}

static void set_limit(char const * limit)
{
#ifdef _WIN32
    _putenv_s("GPUFIT_MAX_CHUNK_SIZE", limit);
#else
    if (*limit)
        setenv("GPUFIT_MAX_CHUNK_SIZE", limit, 1);
    else
        unsetenv("GPUFIT_MAX_CHUNK_SIZE");
#endif
}

int main()
try
{
    set_limit("");
    for (DataLocation location : {HOST, DEVICE})
    for (bool constraints : {false, true})
    for (bool weights : {false, true})
    {
        Info info = make_info(location, constraints, weights);
        std::size_t const one = minimum_budget(info);
        info.n_fits_ = 2;
        std::size_t const two = minimum_budget(info);
        // Fixed allocations are fitted indices, all_finished and optional constraint types.
        check(2 * one - two == sizeof(int) * (3 + (constraints ? 4 : 0)));
        info.n_fits_ = 1;
        check(capacity(info, 0) == 0);
        check(capacity(info, one - 1) == 0);
        check(capacity(info, one) == 1);
        info.user_info_size_ = 256;
        check(capacity(info, location == HOST ? one + 255 : one - 1) == 0);
        check(capacity(info, location == HOST ? one + 256 : one) == 1);
        info.user_info_size_ = std::numeric_limits<std::size_t>::max();
        check(capacity(info, one) == (location == HOST ? 0 : 1));
        info.user_info_size_ = one; // Exactly fills the budget, leaving no space for a fit.
        check(capacity(info, one) == (location == HOST ? 0 : 1));
    }

    Info info = make_info(HOST, false, false);
    std::size_t const one = minimum_budget(info);
    info.n_fits_ = 2;
    std::size_t const per_fit = minimum_budget(info) - one;
    info.n_fits_ = 100;
    std::size_t const budget = one + 26 * per_fit;
    check(capacity(info, budget) == 20); // Decimal rounding of capacity 27.
    total_bytes /= 2;
    check(capacity(info, budget) == 20); // Free memory, not total memory, sets the budget.
    info.n_fits_ = 7;
    check(capacity(info, budget) == 7);
    max_blocks = 6;
    check(capacity(info, budget) == 6);
    max_blocks = 1024;
    set_limit("3");
    check(capacity(info, budget) == 3);
    set_limit("1000");
    check(capacity(info, one) == 1); // An environment cap cannot raise memory capacity.
    set_limit("");
    std::puts("memory budget boundaries and fixed-buffer accounting passed");
    return 0;
}
catch (std::exception const & error)
{
    std::fprintf(stderr, "%s\n", error.what());
    return 1;
}
