#ifndef GPUFIT_GAUSS3DARB_ROT_CUH_INCLUDED
#define GPUFIT_GAUSS3DARB_ROT_CUH_INCLUDED

__device__ void calculate_gauss3d_rot_arb(
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
    // read coordinates stored in user_info
    // stored as 1D raveled arrays: x1, y1, z1, x2, y2, z2, ... xn, yn, zn, num_1, ..., num_n, sx_min, sy_min, sz_min
    // the num_i allow the fits to be of different sizes. If so, image data and coordinates
    // for each fit must still be padded to n_points, but only the first num_i <= n_points will be non-zero.

    // cast user_info to float
    float * user_info_float = (float *) user_info;

    // extract coordinates
    float x, y, z;
    // stride between points is 3*fit_index*n_points
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

    // extract parameters [A, cx, cy, cz, sx, sy, sz, bg, phi, theta, psi]
    float const * p = parameters;
    // Define our Euler angles connecting the body frame to the space/lab frame by
    // r_lab = U_z(phi) * U_y(theta) * U_z(psi) * r_body
    // Consider the z-axis in the body frame. This axis is then orientated at [cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)]
    // in the space frame. i.e. phi, theta are the usual polar angles. psi represents a rotation of the object about its own axis.
    // U_z(phi) = [[cos(phi), -sin(phi), 0], [sin(phi), cos(phi), 0], [0, 0, 1]]
    // U_y(theta) = [[cos(theta), 0, sin(theta)], [0, 1, 0], [-sin(theta), 0, cos(theta)]]
    // Define our rotated gaussian function G_rot(r) = G(E^{-1} * r) where E^{-1} is the inverse Euler matrix
    // conveniently E(phi, theta, psi)^{-1} = E(-psi, -theta, -phi)


    if (point_index < nsize){
        REAL const amp = p[0];
        REAL const cx = p[1];
        REAL const cy = p[2];
        REAL const cz = p[3];
        REAL const sx = p[4];
        REAL const sy = p[5];
        REAL const sz = p[6];
        REAL const bg = p[7];
        REAL const phi = p[8];
        REAL const theta = p[9];
        REAL const psi = p[10];

        // effective sigma^2 values
        REAL const sx_sqr = sx_min * sx_min + sx * sx;
        REAL const sy_sqr = sy_min * sy_min + sy * sy;
        REAL const sz_sqr = sz_min * sz_min + sz * sz;

        // swapping phi/-psi, theta/-theta, psi/-phi will give us the inverse Euler matrix instead of the forward one, as desired

        // Inverse Euler rotation matrix = R^{-1}
        REAL const r00 =  cos(-psi) * cos(-theta) * cos(-phi) - sin(-psi) * sin(-phi);
        REAL const r01 = -cos(-psi) * cos(-theta) * sin(-phi) - sin(-psi) * cos(-phi);
        REAL const r02 =  cos(-psi) * sin(-theta);
        REAL const r10 =  sin(-psi) * cos(-theta) * cos(-phi) + cos(-psi) * sin(-phi);
        REAL const r11 = -sin(-psi) * cos(-theta) * sin(-phi) + cos(-psi) * cos(-phi);
        REAL const r12 =  sin(-psi) * sin(-theta);
        REAL const r20 = -sin(-theta) * cos(-phi);
        REAL const r21 =  sin(-theta) * sin(-phi);
        REAL const r22 =  cos(-theta);

        // dR^{-1}/dphi
        REAL const dp00 = -1 * (-cos(-psi) * cos(-theta) * sin(-phi) - sin(-psi) * cos(-phi));
        REAL const dp01 = -1 * (-cos(-psi) * cos(-theta) * cos(-phi) + sin(-psi) * sin(-phi));
        REAL const dp02 = -1 * 0;
        REAL const dp10 = -1 * (-sin(-psi) * cos(-theta) * sin(-phi) + cos(-psi) * cos(-phi));
        REAL const dp11 = -1 * (-sin(-psi) * cos(-theta) * cos(-phi) - cos(-psi) * sin(-phi));
        REAL const dp12 = -1 * 0;
        REAL const dp20 = -1 * sin(-theta) * sin(-phi);
        REAL const dp21 = -1 * sin(-theta) * cos(-phi);
        REAL const dp22 = -1 * 0;

        // dR^{-1}/dtheta
        REAL const dt00 = -1 * -cos(-psi) * sin(-theta) * cos(-phi);
        REAL const dt01 = -1 * cos(-psi) * sin(-theta) * sin(-phi);
        REAL const dt02 = -1 * cos(-psi) * cos(-theta);
        REAL const dt10 = -1 * -sin(-psi) * sin(-theta) * cos(-phi);
        REAL const dt11 = -1 * sin(-psi) * sin(-theta) * sin(-phi);
        REAL const dt12 = -1 * sin(-psi) * cos(-theta);
        REAL const dt20 = -1 * -cos(-theta) * cos(-phi);
        REAL const dt21 = -1 * cos(-theta) * sin(-phi);
        REAL const dt22 = - 1 * -sin(-theta);

        // dR^{-1}/dpsi
        REAL const dx00 = -1 * (-sin(-psi) * cos(-theta) * cos(-phi) - cos(-psi) * sin(-phi));
        REAL const dx01 = -1 * (sin(-psi) * cos(-theta) * sin(-phi) - cos(-psi) * cos(-phi));
        REAL const dx02 = -1 * -sin(-psi) * sin(-theta);
        REAL const dx10 = -1 * (cos(-psi) * cos(-theta) * cos(-phi) - sin(-psi) * sin(-phi));
        REAL const dx11 = -1 * (-cos(-psi) * cos(-theta) * sin(-phi) - sin(-psi) * cos(-phi));
        REAL const dx12 = -1 * cos(-psi) * sin(-theta);
        REAL const dx20 = -1 * 0;
        REAL const dx21 = -1 * 0;
        REAL const dx22 = -1 * 0;

        // r_rot = R^{-1} * r
        REAL const xrot = (x - cx) * r00 + (y - cy) * r01 + (z - cz) * r02;
        REAL const yrot = (x - cx) * r10 + (y - cy) * r11 + (z - cz) * r12;
        REAL const zrot = (x - cx) * r20 + (y - cy) * r21 + (z - cz) * r22;
        REAL const ex = expf(-0.5f * (xrot * xrot / sx_sqr +
                                      yrot * yrot / sy_sqr +
                                      zrot * zrot / sz_sqr
                                      ));

        // compute function value
        value[point_index] =  bg + amp * ex;

        //compute function partial derivatives
        derivative[point_index + 0*n_points] = ex;
        // df/dcx = df/dx * dxrot / dcx + df/dy * dyrot / dcx + df/dz * dzrot / dcx
    	derivative[point_index + 1*n_points] = amp * ex * (xrot / sx_sqr * r00 + yrot / sy_sqr * r10 + zrot / sz_sqr * r20);
        derivative[point_index + 2*n_points] = amp * ex * (xrot / sx_sqr * r01 + yrot / sy_sqr * r11 + zrot / sz_sqr * r21);
    	derivative[point_index + 3*n_points] = amp * ex * (xrot / sx_sqr * r02 + yrot / sy_sqr * r12 + zrot / sz_sqr * r22);
	    // df/dsigma
	    derivative[point_index + 4*n_points] = amp * ex * sx / sx_sqr * xrot * xrot;
	    derivative[point_index + 5*n_points] = amp * ex * sy / sy_sqr * yrot * yrot;
	    derivative[point_index + 6*n_points] = amp * ex * sz / sz_sqr * zrot * zrot;
	    // df/dbackground
	    derivative[point_index + 7*n_points] = 1;
	    // df / dphi = df/dx * dxrot / dphi + df/dy * dyrot / dphi + df/dz * dzrot / dphi
	    derivative[point_index + 8*n_points] =  -amp * ex * (xrot / sx_sqr * ((x - cx) * dp00 + (y - cy) * dp01  + (z - cz) * dp02) +
	                                                         yrot / sy_sqr * ((x - cx) * dp10 + (y - cy) * dp11  + (z - cz) * dp12) +
	                                                         zrot / sz_sqr * ((x - cx) * dp20 + (y - cy) * dp21  + (z - cz) * dp22));
	    derivative[point_index + 9*n_points] =  -amp * ex * (xrot / sx_sqr * ((x - cx) * dt00 + (y - cy) * dt01  + (z - cz) * dt02) +
	                                                         yrot / sy_sqr * ((x - cx) * dt10 + (y - cy) * dt11  + (z - cz) * dt12) +
	                                                         zrot / sz_sqr * ((x - cx) * dt20 + (y - cy) * dt21  + (z - cz) * dt22));
	    derivative[point_index + 10*n_points] = -amp * ex * (xrot / sx_sqr * ((x - cx) * dx00 + (y - cy) * dx01  + (z - cz) * dx02) +
	                                                         yrot / sy_sqr * ((x - cx) * dx10 + (y - cy) * dx11  + (z - cz) * dx12) +
	                                                         zrot / sz_sqr * ((x - cx) * dx20 + (y - cy) * dx21  + (z - cz) * dx22));

	}
    else {
    	 value[point_index] = 0;
    	 derivative[point_index + 0*n_points] = 0;
	     derivative[point_index + 1*n_points] = 0;
	     derivative[point_index + 2*n_points] = 0;
	     derivative[point_index + 3*n_points] = 0;
	     derivative[point_index + 4*n_points] = 0;
	     derivative[point_index + 5*n_points] = 0;
	     derivative[point_index + 6*n_points] = 0;
	     derivative[point_index + 7*n_points] = 0;
	     derivative[point_index + 8*n_points] = 0;
	     derivative[point_index + 9*n_points] = 0;
	     derivative[point_index + 10*n_points] = 0;
	 }
}

#endif