// Execute the existing microscopy model expressions on the host and compare
// their analytic derivatives with independent central finite differences.
#include <cmath>
#include <cstdio>
#include <cstddef>
#define __device__
#define REAL float
#include "Gpufit/models/gauss_3d_arb.cuh"
#include "Gpufit/models/gauss_lor_3d_arb.cuh"
#include "Gpufit/models/gauss_3d_rot_arb.cuh"
#include "Gpufit/models/gauss_3d_asymmetric_arb.cuh"

using Model = void (*)(float const *, int, int, float *, float *, int, int, int, char *, std::size_t);
int main() {
    Model models[] = {calculate_gauss3d_arb, calculate_gauss_lor_3d_arb,
                      calculate_gauss3d_rot_arb, calculate_gauss3d_asym_arb};
    int counts[] = {7,7,11,8};
    float starts[][11] = {{100,.1,.2,.3,.4,.8,5}, {100,.1,.2,.3,.4,.8,5},
                         {100,.1,.2,.3,.4,.5,.8,5,.2,.3,.4}, {100,.1,.2,.3,.4,.5,.8,5}};
    float info[] = {.72f,.42f,-.2f,1,.08f,.09f,.12f};
    int failed = 0;
    for(int m=0; m<4; ++m) {
        float value[1], derivative[11], scratch[11];
        models[m](starts[m],1,1,value,derivative,0,0,0,(char*)info,sizeof(info));
        for(int j=0; j<counts[m]; ++j) {
            float old = starts[m][j];
            float step = 0.001f * std::fmax(1.f,std::fabs(old));
            float upper[1], lower[1];
            starts[m][j] = old + step;
            models[m](starts[m],1,1,upper,scratch,0,0,0,(char*)info,sizeof(info));
            starts[m][j] = old - step;
            models[m](starts[m],1,1,lower,scratch,0,0,0,(char*)info,sizeof(info));
            starts[m][j] = old;
            float numerical = (upper[0]-lower[0])/(2*step);
            if(std::fabs(numerical-derivative[j]) > 0.01f + 0.005f*std::fabs(numerical)) {
                std::printf("model %d parameter %d: analytic %.6g numerical %.6g\n",15+m,j,derivative[j],numerical);
                ++failed;
            }
        }
    }
    return failed ? 1 : 0;
}
