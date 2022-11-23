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
        REAL const phi = p[8];
        REAL const theta = p[9];
        REAL const psi = p[10];

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
        REAL const xrot = (x - p[1]) * r00 + (y - p[2]) * r01 + (z - p[3]) * r02;
        REAL const yrot = (x - p[1]) * r10 + (y - p[2]) * r11 + (z - p[3]) * r12;
        REAL const zrot = (x - p[1]) * r20 + (y - p[2]) * r21 + (z - p[3]) * r22;
        REAL const ex = expf(-0.5f * (xrot * xrot / p[4] / p[4] + yrot * yrot / p[5] / p[5] + zrot * zrot / p[6] / p[6]));

        // compute function value
        value[point_index] =  p[7] + p[0] * ex;

        //compute function partial derivatives
        derivative[point_index + 0*n_points] = ex;
        // df/dcx = df/dx * dxrot / dcx + df/dy * dyrot / dcx + df/dz * dzrot / dcx
    	derivative[point_index + 1*n_points] = p[0] * ex * (xrot / p[4] / p[4] * r00 + yrot / p[5] / p[5] * r10 + zrot / p[6] / p[6] * r20);
        derivative[point_index + 2*n_points] = p[0] * ex * (xrot / p[4] / p[4] * r01 + yrot / p[5] / p[5] * r11 + zrot / p[6] / p[6] * r21);
    	derivative[point_index + 3*n_points] = p[0] * ex * (xrot / p[4] / p[4] * r02 + yrot / p[5] / p[5] * r12 + zrot / p[6] / p[6] * r22);
	    derivative[point_index + 4*n_points] = p[0] * ex * xrot * xrot / p[4] / p[4] / p[4];
	    derivative[point_index + 5*n_points] = p[0] * ex * yrot * yrot / p[5] / p[5] / p[5];
	    derivative[point_index + 6*n_points] = p[0] * ex * zrot * zrot / p[6] / p[6] / p[6];
	    derivative[point_index + 7*n_points] = 1;
	    // df / dphi = df/dx * dxrot / dphi + df/dy * dyrot / dphi + df/dz * dzrot / dphi
	    derivative[point_index + 8*n_points] = -p[0] * ex * (xrot / p[4] / p[4] * ((x - p[1]) * dp00 + (y - p[2]) * dp01  + (z - p[3]) * dp02) +
	                                                         yrot / p[5] / p[5] * ((x - p[1]) * dp10 + (y - p[2]) * dp11  + (z - p[3]) * dp12) +
	                                                         zrot / p[6] / p[6] * ((x - p[1]) * dp20 + (y - p[2]) * dp21  + (z - p[3]) * dp22));
	    derivative[point_index + 9*n_points] = -p[0] * ex * (xrot / p[4] / p[4] * ((x - p[1]) * dt00 + (y - p[2]) * dt01  + (z - p[3]) * dt02) +
	                                                         yrot / p[5] / p[5] * ((x - p[1]) * dt10 + (y - p[2]) * dt11  + (z - p[3]) * dt12) +
	                                                         zrot / p[6] / p[6] * ((x - p[1]) * dt20 + (y - p[2]) * dt21  + (z - p[3]) * dt22));
	    derivative[point_index + 10*n_points] = -p[0] * ex * (xrot / p[4] / p[4] * ((x - p[1]) * dx00 + (y - p[2]) * dx01  + (z - p[3]) * dx02) +
	                                                         yrot / p[5] / p[5] * ((x - p[1]) * dx10 + (y - p[2]) * dx11  + (z - p[3]) * dx12) +
	                                                         zrot / p[6] / p[6] * ((x - p[1]) * dx20 + (y - p[2]) * dx21  + (z - p[3]) * dx22));

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