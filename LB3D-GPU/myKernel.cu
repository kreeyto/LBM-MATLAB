#include <cuda_runtime.h>
#include <math.h>

// remember to generate ptx file to communicate with matlab
    /* nvcc -ptx myKernel.cu -o myKernel.ptx /*

// before simulation loop in matlab
    /* k = parallel.gpu.CUDAKernel('myKernel.ptx', 'myKernel.cu', 'updatePhi');
        k.ThreadBlockSize = [8, 8, 8];
        k.GridSize = [ceil(nx/8), ceil(ny/8), ceil(nz/8)];
        phi_gpu = gpuArray(phi);
        g_gpu = gpuArray(g); */

// in simulation loop
    /* phi_gpu = feval(k, phi_gpu, g_gpu, nx, ny, nz, gpoints);
        phi = gather(phi_gpu); */

__global__ void updatePhi(double *phi, double *g, int nx, int ny, int nz, int gpoints) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

    if (i < nx - 1 && j < ny - 1 && k < nz - 1) {
        double sumG = 0.0;
        for (int l = 0; l < gpoints; l++) {
            sumG += g[(i * ny * nz + j * nz + k) * gpoints + l];
        }
        phi[i * ny * nz + j * nz + k] = sumG;
    }
}
