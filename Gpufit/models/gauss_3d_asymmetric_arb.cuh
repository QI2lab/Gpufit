#ifndef GPUFIT_GAUSS23ASYMARB_CUH_INCLUDED
#define GPUFIT_GAUSS23ASYMARB_CUH_INCLUDED

__device__ void calculate_gauss3d_asym_arb(
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
    // stored as arrays x1, y1, z1, x2, y2, z2, ... xn, yn, zn, num_1, ..., num_n, sxy_min, sz_min
    // the num_i allow the fits to be of different sizes. If so, image data and coordinates
    // for each fit must still be padded to n_points, but only the first num_i <= n_points will be non-zero.

    // cast user info to float
    float * user_info_float = (float *) user_info;

    // extract coordinates
    // stride between points is 3*fit_index*n_points
    float x, y, z;
    x = user_info_float[point_index + 3 * fit_index * n_points];
    y = user_info_float[point_index + n_points + 3 * fit_index * n_points];
    z = user_info_float[point_index + 2*n_points + 3 * fit_index * n_points];

    // extract sizes
    // todo: probably better to store this as integers, but for a first test store as floats and cast
    //int * user_info_int = (int *) user_info_float[3 * n_points * n_fits];
    int nsize = (int) user_info_float[3 * n_points * n_fits + fit_index];

    // extract minimum sigmas
    float sx_min = user_info_float[3 * n_points * n_fits + n_fits];
	float sy_min = user_info_float[3 * n_points * n_fits + n_fits + 1];
    float sz_min = user_info_float[3 * n_points * n_fits + n_fits + 2];

    // extract parameters
    float const * p = parameters;

    //compute function partial derivatives
    if (point_index < nsize){
        REAL const amp = p[0];
        REAL const cx = p[1];
        REAL const cy = p[2];
        REAL const cz = p[3];
        REAL const sx = p[4];
		REAL const sy = p[5];
        REAL const sz = p[6];
        REAL const bg = p[7];

       REAL const ex = expf(-0.5f * ((x - cx) * (x - cx) / (sx_min * sx_min + sx * sx) +
                                     (y - cy) * (y - cy) / (sy_min * sy_min + sy * sy) +
                                     (z - cz) * (z - cz) / (sz_min * sz_min + sz * sz)
                                     ));

        // compute function value
        value[point_index] =  bg + amp * ex;

        // compute derivative value
        derivative[point_index + 0*n_points] = ex;
    	derivative[point_index + 1*n_points] = amp * ex * (x - cx) / (sx_min * sx_min + sx * sx);
        derivative[point_index + 2*n_points] = amp * ex * (y - cy) / (sy_min * sy_min + sy * sy);
    	derivative[point_index + 3*n_points] = amp * ex * (z - cz) / (sz_min * sz_min + sz * sz);
	    derivative[point_index + 4*n_points] = amp * ex * sx / (sx_min * sx_min + sx * sx) / (sx_min * sx_min + sx * sx) * (x - cx) * (x - cx);           
		derivative[point_index + 5*n_points] = amp * ex * sy / (sy_min * sy_min + sy * sy) / (sy_min * sy_min + sy * sy) * (y - cy) * (y - cy);
	    derivative[point_index + 6*n_points] = amp * ex * sz / (sz_min * sz_min + sz * sz) / (sz_min * sz_min + sz * sz) * (z - cz) * (z - cz);
	    derivative[point_index + 7*n_points] = 1.0;
	}
    else{
        value[point_index] = 0;

	    derivative[point_index + 0*n_points] = 0;
	    derivative[point_index + 1*n_points] = 0;
	    derivative[point_index + 2*n_points] = 0;
	    derivative[point_index + 3*n_points] = 0;
	    derivative[point_index + 4*n_points] = 0;
	    derivative[point_index + 5*n_points] = 0;
	    derivative[point_index + 6*n_points] = 0;
		derivative[point_index + 7*n_points] = 0;
     	}

}

#endif