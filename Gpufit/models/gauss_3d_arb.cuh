#ifndef GPUFIT_GAUSS2DARB_CUH_INCLUDED
#define GPUFIT_GAUSS2DARB_CUH_INCLUDED

__device__ void calculate_gauss3d_arb(
    float const * parameters,
    int const n_fits,
    int const n_points,
    float * value,
    float * derivative,
    int const point_index,
    int const fit_index,
    int const chunk_index,
    char * user_info,
    std::size_t const user_info_size)
{

    // read coordinates stored in user data
    // stored as arrays x1, y1, z1, x2, y2, z2, ...
    // stride between points is 3*fit_index*n_points
    float * user_info_float = (float *) user_info;
    float x, y, z;
    x = user_info_float[point_index + 3 * fit_index * n_points];
    y = user_info_float[point_index + n_points + 3 * fit_index * n_points];
    z = user_info_float[point_index + 2*n_points + 3 * fit_index * n_points];

    float const * p = parameters;

    // compute function value

    value[point_index] =  p[6] + p[0] * expf(-0.5f * ((x - p[1]) * (x - p[1]) / p[4] / p[4] + (y - p[2]) * (y - p[2]) / p[4] / p[4] + (z - p[3]) * (z - p[3]) / p[5] / p[5]));

    //compute function partial derivatives
    
    derivative[point_index + 0*n_points] = expf(-0.5f * ((x - p[1]) * (x - p[1]) / p[4] / p[4] + (y - p[2]) * (y - p[2]) / p[4] / p[4] + (z - p[3]) * (z - p[3]) / p[5] / p[5]));
    derivative[point_index + 1*n_points] = p[0] * derivative[point_index] * (x - p[1]) / p[4] / p[4];
    derivative[point_index + 2*n_points] = p[0] * derivative[point_index] * (y - p[2]) / p[4] / p[4];
    derivative[point_index + 3*n_points] = p[0] * derivative[point_index] * (z - p[3]) / p[5] / p[5];
    derivative[point_index + 4*n_points] = p[0] * derivative[point_index] * ((x - p[1]) * (x - p[1]) / p[4] / p[4] / p[4] + (y - p[2]) * (y - p[2]) / p[4] / p[4] / p[4]);
    derivative[point_index + 5*n_points] = p[0] * derivative[point_index] * (z - p[3]) * (z - p[3]) / p[5] / p[5] / p[5];
    derivative[point_index + 6*n_points] = 1.0;

}

#endif