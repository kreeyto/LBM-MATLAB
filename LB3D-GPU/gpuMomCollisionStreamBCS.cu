/* 
 * ----------------------------------------------------------------------------------------------------
 * TO FICANDO PUTO COM ESSE CODIGO
 * CUDA_ERROR_LAUNCH_FAILED
 * clc; !nvcc -arch=sm_86 -ptx gpuMomCollisionStreamBCS.cu -o gpuMomCollisionStreamBCS.ptx
 * ----------------------------------------------------------------------------------------------------
 */
 
#include "cuda_runtime.h"
#include <math.h>

__global__ void momCollisionStreamBCS(
    float *f, float *g, float *phi, float *rho, const float *w, const float *w_g,
    const float *cix, const float *ciy, const float *ciz,
    float *mod_grad, float *normx, float *normy, float *normz, float *indicator,
    float *curvature, float *ffx, float *ffy, float *ffz, 
    float *ux, float *uy, float *uz,
    float *pxx, float *pyy, float *pzz, float *pxy, float *pxz, float *pyz,
    int nx, int ny, int nz, int fpoints, int gpoints, 
    float sigma, float cssq, float omega, float sharp_c, int nsteps
) {
    
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    // int idx = i + nx * (j + ny * k);
    #define F_IDX(i,j,k,l) ((i) + nx * ((j) + ny * ((k) + nz * (l))))
    #define IDX3D(i,j,k) ((i) + nx * ((j) + ny * (k)))

    float fneq[19];
    float grad_fix, grad_fiy, grad_fiz, uu, udotc, HeF, feq, Hi;

    for (int t = 0; t < nsteps; t++) {

        // nested loops become hard to manage
        if (i > 0 && i < nx-1 && j > 0 && j < ny-1 && k > 0 && k < nz-1) {
            for (int l = 0; l < gpoints; l++) {
                phi[IDX3D(i,j,k)] += g[F_IDX(i, j, k, l)];
            }            
        }

        if (i > 0 && i < nx-1 && j > 0 && j < ny-1 && k > 0 && k < nz-1) {
            grad_fix = 0;
            grad_fiy = 0;
            grad_fiz = 0;
            for (int l = 0; l < fpoints; l++) {
                grad_fix += 3 * w[l] * cix[l] * phi[IDX3D(i + static_cast<int>(cix[l]),
                                                          j + static_cast<int>(ciy[l]),
                                                          k + static_cast<int>(ciz[l]))];
                grad_fiy += 3 * w[l] * ciy[l] * phi[IDX3D(i + static_cast<int>(cix[l]),
                                                          j + static_cast<int>(ciy[l]),
                                                          k + static_cast<int>(ciz[l]))];
                grad_fiz += 3 * w[l] * ciz[l] * phi[IDX3D(i + static_cast<int>(cix[l]),
                                                          j + static_cast<int>(ciy[l]),
                                                          k + static_cast<int>(ciz[l]))];
            }
            mod_grad[IDX3D(i,j,k)] = sqrt(pow(grad_fix,2) + pow(grad_fiy,2) + pow(grad_fiz,2));
            normx[IDX3D(i,j,k)] = grad_fix / (mod_grad[IDX3D(i,j,k)] + 1e-8);
            normy[IDX3D(i,j,k)] = grad_fiy / (mod_grad[IDX3D(i,j,k)] + 1e-8);
            normz[IDX3D(i,j,k)] = grad_fiz / (mod_grad[IDX3D(i,j,k)] + 1e-8);
            indicator[IDX3D(i,j,k)] = sqrt(pow(grad_fix,2) + pow(grad_fiy,2) + pow(grad_fiz,2));
        }

        if (i > 0 && i < nx-1 && j > 0 && j < ny-1 && k > 0 && k < nz-1) {
            curvature[IDX3D(i,j,k)] = 0;
            for (int l = 0; l < fpoints; l++) {
                curvature[IDX3D(i,j,k)] -= 3 * w[l] * 
                (cix[l] * normx[IDX3D(i + static_cast<int>(cix[l]),
                                      j + static_cast<int>(ciy[l]),
                                      k + static_cast<int>(ciz[l]))] + 
                 ciy[l] * normy[IDX3D(i + static_cast<int>(cix[l]),
                                      j + static_cast<int>(ciy[l]),
                                      k + static_cast<int>(ciz[l]))] + 
                 ciz[l] * normz[IDX3D(i + static_cast<int>(cix[l]),
                                      j + static_cast<int>(ciy[l]),
                                      k + static_cast<int>(ciz[l]))] 
                );
            }
            ffx[IDX3D(i,j,k)] = sigma * curvature[IDX3D(i,j,k)] * normx[IDX3D(i,j,k)] * indicator[IDX3D(i,j,k)];
            ffy[IDX3D(i,j,k)] = sigma * curvature[IDX3D(i,j,k)] * normy[IDX3D(i,j,k)] * indicator[IDX3D(i,j,k)];
            ffz[IDX3D(i,j,k)] = sigma * curvature[IDX3D(i,j,k)] * normz[IDX3D(i,j,k)] * indicator[IDX3D(i,j,k)];
        }

        if (i > 0 && i < nx-1 && j > 0 && j < ny-1 && k > 0 && k < nz-1) {
            ux[IDX3D(i,j,k)] = (
                (f[F_IDX(i,j,k,1)] + f[F_IDX(i,j,k,15)] + f[F_IDX(i,j,k,9)] + f[F_IDX(i,j,k,7)] + f[F_IDX(i,j,k,13)]) -
                (f[F_IDX(i,j,k,2)] + f[F_IDX(i,j,k,10)] + f[F_IDX(i,j,k,16)] + f[F_IDX(i,j,k,14)] + f[F_IDX(i,j,k,7)])
            ) / rho[IDX3D(i,j,k)] +
            ffx[IDX3D(i,j,k)] * 0.5 / rho[IDX3D(i,j,k)];

            uy[IDX3D(i,j,k)] = (
                (f[F_IDX(i,j,k,3)] + f[F_IDX(i,j,k,7)] + f[F_IDX(i,j,k,14)] + f[F_IDX(i,j,k,17)] + f[F_IDX(i,j,k,11)]) -
                (f[F_IDX(i,j,k,4)] + f[F_IDX(i,j,k,13)] + f[F_IDX(i,j,k,8)] + f[F_IDX(i,j,k,12)] + f[F_IDX(i,j,k,18)])
            ) / rho[IDX3D(i,j,k)] +
            ffy[IDX3D(i,j,k)] * 0.5 / rho[IDX3D(i,j,k)];

            uz[IDX3D(i,j,k)] = (
                (f[F_IDX(i,j,k,6)] + f[F_IDX(i,j,k,15)] + f[F_IDX(i,j,k,10)] + f[F_IDX(i,j,k,17)] + f[F_IDX(i,j,k,12)]) -
                (f[F_IDX(i,j,k,5)] + f[F_IDX(i,j,k,9)] + f[F_IDX(i,j,k,16)] + f[F_IDX(i,j,k,11)] + f[F_IDX(i,j,k,18)])
            ) / rho[IDX3D(i,j,k)] +
            ffz[IDX3D(i,j,k)] * 0.5 / rho[IDX3D(i,j,k)];

            uu = 0.5 * (pow(ux[IDX3D(i,j,k)],2) + pow(uy[IDX3D(i,j,k)],2) + pow(uz[IDX3D(i,j,k)],2)) / cssq;

            for (int l = 0; l < fpoints; l++) {
                rho[IDX3D(i,j,k)] += f[F_IDX(i, j, k, l)];
            }   

            for (int l = 0; l < fpoints; l++) {
                udotc = (ux[IDX3D(i,j,k)] * cix[l] + uy[IDX3D(i,j,k)] * ciy[l] + uz[IDX3D(i,j,k)] * ciz[l]) / cssq;
                 HeF = (w[l] * (rho[IDX3D(i,j,k)] + rho[IDX3D(i,j,k)] * (udotc + 0.5 * pow(udotc,2) - uu)))
                        * ((cix[l] - ux[IDX3D(i,j,k)]) * ffx[IDX3D(i,j,k)] + 
                            (ciy[l] - uy[IDX3D(i,j,k)]) * ffy[IDX3D(i,j,k)] + 
                            (ciz[l] - uz[IDX3D(i,j,k)]) * ffz[IDX3D(i,j,k)] 
                        ) / (rho[IDX3D(i,j,k)] * cssq);
                feq = w[l] * (rho[IDX3D(i,j,k)] + rho[IDX3D(i,j,k)] * (udotc + 0.5 * pow(udotc,2) - uu)) - 0.5 * HeF;
                fneq[l] = f[F_IDX(i,j,k,l)] - feq;
            }

            pxx[IDX3D(i,j,k)] = fneq[2] + fneq[3] + fneq[8] + fneq[9] + fneq[10] + fneq[11] + fneq[14] + fneq[15] + fneq[16] + fneq[17];
            pyy[IDX3D(i,j,k)] = fneq[4] + fneq[5] + fneq[8] + fneq[9] + fneq[12] + fneq[13] + fneq[14] + fneq[15] + fneq[18] + fneq[19];
            pzz[IDX3D(i,j,k)] = fneq[6] + fneq[7] + fneq[10] + fneq[11] + fneq[12] + fneq[13] + fneq[16] + fneq[17] + fneq[18] + fneq[19];
            pxy[IDX3D(i,j,k)] = fneq[8] + fneq[9] - fneq[14] - fneq[15];
            pxz[IDX3D(i,j,k)] = fneq[10] + fneq[11] - fneq[16] - fneq[17];
            pyz[IDX3D(i,j,k)] = fneq[12] + fneq[13] - fneq[18] - fneq[19];
        }

        if (i > 0 && i < nx-1 && j > 0 && j < ny-1 && k > 0 && k < nz-1) {

            float uu = 0.5 * (pow(ux[IDX3D(i,j,k)], 2) + pow(uy[IDX3D(i,j,k)], 2) + pow(uz[IDX3D(i,j,k)], 2)) / cssq;

            for (int l = 0; l < fpoints; l++) {
                udotc = (ux[IDX3D(i,j,k)] * cix[l] + uy[IDX3D(i,j,k)] * ciy[l] + uz[IDX3D(i,j,k)] * ciz[l]) / cssq;
                feq = w[l] * (rho[IDX3D(i,j,k)] + rho[IDX3D(i,j,k)] * (udotc + 0.5 * pow(udotc, 2) - uu));
                HeF = 0.5 * (w[l] * (rho[IDX3D(i,j,k)] + rho[IDX3D(i,j,k)] * (udotc + 0.5 * pow(udotc, 2) - uu)))
                        * ((cix[l] - ux[IDX3D(i,j,k)]) * ffx[IDX3D(i,j,k)] + 
                            (ciy[l] - uy[IDX3D(i,j,k)]) * ffy[IDX3D(i,j,k)] + 
                            (ciz[l] - uz[IDX3D(i,j,k)]) * ffz[IDX3D(i,j,k)] 
                        ) / (rho[IDX3D(i,j,k)] * cssq);
                float fneq = (cix[l] * cix[l] - cssq) * pxx[IDX3D(i,j,k)] + 
                            (ciy[l] * ciy[l] - cssq) * pyy[IDX3D(i,j,k)] + 
                            (ciz[l] * ciz[l] - cssq) * pzz[IDX3D(i,j,k)] + 
                            2 * cix[l] * ciy[l] * pxy[IDX3D(i,j,k)] + 
                            2 * cix[l] * ciz[l] * pxz[IDX3D(i,j,k)] + 
                            2 * ciy[l] * ciz[l] * pyz[IDX3D(i,j,k)]; // float fneq ~= fneq[l]
                f[F_IDX(i + static_cast<int>(cix[l]),
                        j + static_cast<int>(ciy[l]),
                        k + static_cast<int>(ciz[l]),
                        l)] = feq + (1 - omega) * (w[l] / (2 * pow(cssq, 2))) * fneq + HeF;
            }

            for (int l = 0; l < gpoints; l++) {
                udotc = (ux[IDX3D(i,j,k)] * cix[l] + uy[IDX3D(i,j,k)] * ciy[l] + uz[IDX3D(i,j,k)] * ciz[l]) / cssq;
                feq = w_g[l] * phi[IDX3D(i,j,k)] * (1 + udotc);
                Hi = sharp_c * phi[IDX3D(i,j,k)] * (1 - phi[IDX3D(i,j,k)]) * (cix[l] * normx[IDX3D(i,j,k)] + ciy[l] * normy[IDX3D(i,j,k)] + ciz[l] * normz[IDX3D(i,j,k)]); 
                g[F_IDX(i,j,k,l)] = feq + w_g[l] * Hi;
            }

        }

        if (i > 0 && i < nx-1 && j > 0 && j < ny-1 && k > 0 && k < nz-1) {
            for (int l = 0; l < gpoints; l++) {
                g[F_IDX(i,j,k,l)] = g[F_IDX(i + static_cast<int>(cix[l]),
                                            j + static_cast<int>(ciy[l]),
                                            k + static_cast<int>(ciz[l]),
                                            l)];
            }
        }

        if (i == 0 || i == nx-1 || j == 0 || j == ny-1 || k == 0 || k == ny-1) {
            for (int l = 0; l < fpoints; l++) {
                if (i + static_cast<int>(cix[l]) >= 0 && j + static_cast<int>(ciy[l]) >= 0 && k + static_cast<int>(ciz[l]) >= 0) {
                    f[F_IDX(i + static_cast<int>(cix[l]),
                            j + static_cast<int>(ciy[l]),
                            k + static_cast<int>(ciz[l]),
                            l)] = rho[IDX3D(i,j,k)] * w[l];
                }
            }
            for (int l = 0; l < gpoints; l++) {
                if (i + static_cast<int>(cix[l]) >= 0 && j + static_cast<int>(ciy[l]) >= 0 && k + static_cast<int>(ciz[l]) >= 0) {
                    g[F_IDX(i + static_cast<int>(cix[l]),
                            j + static_cast<int>(ciy[l]),
                            k + static_cast<int>(ciz[l]),
                            l)] = phi[IDX3D(i,j,k)] * w_g[l];
                }
            }
        }

        if (j == 0) {
            phi[IDX3D(i,j,k)] = phi[IDX3D(i,j+1,k)];
        }
        if (j == ny - 1) {
            phi[IDX3D(i,j,k)] = phi[IDX3D(i,j-1,k)]; 
        }
        if (k == 0) {
            phi[IDX3D(i,j,k)] = phi[IDX3D(i,j,k+1)];
        }
        if (k == nz - 1) {
            phi[IDX3D(i,j,k)] = phi[IDX3D(i,j,k-1)]; 
        }
        if (i == 0) {
            phi[IDX3D(i,j,k)] = phi[IDX3D(i+1,j,k)];
        }
        if (i == nx - 1) {
            phi[IDX3D(i,j,k)] = phi[IDX3D(i-1,j,k)]; 
        }

    } // end of time loop
}

