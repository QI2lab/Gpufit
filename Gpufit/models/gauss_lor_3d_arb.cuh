#ifndef GPUFIT_GAUSS_LOR_ARB_CUH_INCLUDED
#define GPUFIT_GAUSS_LOR_ARB_CUH_INCLUDED

__device__ void calculate_gauss_lor_3d_arb(
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
    // stored as arrays x1, y1, z1, x2, y2, z2, ... xn, yn, zn, num_1, ..., num_n
    // the num_i allow the fits to be of different sizes. If so, image data and coordinates
    // for each fit must still be padded to n_points, but only the first num_i <= n_points will be non-zero.

    // extract coordinates
    // stride between points is 3*fit_index*n_points
    float * user_info_float = (float *) user_info;
    float x, y, z;
    x = user_info_float[point_index + 3 * fit_index * n_points];
    y = user_info_float[point_index + n_points + 3 * fit_index * n_points];
    z = user_info_float[point_index + 2*n_points + 3 * fit_index * n_points];

    // extract sizes
    // todo: probably better to store this as integers, but for a first test store as floats and cast
    //int * user_info_int = (int *) user_info_float[3 * n_points * n_fits];
    int nsize = (int) user_info_float[3 * n_points * n_fits + fit_index];

    // extract parameters
    float const * p = parameters;

    // compute function value
    if (point_index < nsize){
		float lor_factor = (1 + (z - p[3]) * (z - p[3]) / p[5] / p[5]);
        value[point_index] =  p[6] + p[0] * expf(-0.5f * ((x - p[1]) * (x - p[1]) + (y - p[2]) * (y - p[2])) / p[4] / p[4] / lor_factor) / lor_factor;
	}
    else {
    	 value[point_index] = 0;
	 }

    //compute function partial derivatives
    if (point_index < nsize){
		float lor_factor = (1 + (z - p[3]) * (z - p[3]) / p[5] / p[5]);
		float r_sqr = (x - p[1]) * (x - p[1]) + (y - p[2]) * (y - p[2]);
		
        derivative[point_index + 0*n_points] = expf(-0.5f * ((x - p[1]) * (x - p[1]) + (y - p[2]) * (y - p[2])) / p[4] / p[4] / lor_factor) / lor_factor;
    	derivative[point_index + 1*n_points] = p[0] * derivative[point_index] * (x - p[1]) / p[4] / p[4] / lor_factor;
        derivative[point_index + 2*n_points] = p[0] * derivative[point_index] * (y - p[2]) / p[4] / p[4] / lor_factor;
    	derivative[point_index + 3*n_points] = p[0] * derivative[point_index] * (-(z - p[3]) / p[5] / p[5] / lor_factor / lor_factor * 0.5f * r_sqr / p[4] / p[4] + 1.0f / lor_factor * 2.0f * (z - p[3]) / p[5] / p[5]);
		derivative[point_index + 4*n_points] = p[0] * derivative[point_index] * r_sqr / lor_factor / p[4] / p[4] / p[4];
		derivative[point_index + 5*n_points] = p[0] * derivative[point_index] * (r_sqr * 0.5f / p[4] / p[4] / lor_factor / lor_factor * 2.0f * (z - p[3]) * (z - p[3]) / p[5] / p[5] / p[5] + 1.0f / lor_factor * 2.0f * (z - p[3]) * (z - p[3]) / p[5] / p[5] / p[5]);
		derivative[point_index + 6*n_points] = 1.0;
	}
    else{
		derivative[point_index + 0*n_points] = 0;
		derivative[point_index + 1*n_points] = 0;
		derivative[point_index + 2*n_points] = 0;
		derivative[point_index + 3*n_points] = 0;
		derivative[point_index + 4*n_points] = 0;
		derivative[point_index + 5*n_points] = 0;
		derivative[point_index + 6*n_points] = 0;
    }

}

#endif