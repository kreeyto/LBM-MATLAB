#include <cuda_runtime.h>
#include <math.h>

// remember to generate ptx file to communicate with matlab
    /*  nvcc -ptx myKernel.cu -o myKernel.ptx  */

// before simulation loop in matlab
    /*  k = parallel.gpu.CUDAKernel('myKernel.ptx', 'myKernel.cu', 'updatePhi');
        k.ThreadBlockSize = [8, 8, 8];
        k.GridSize = [ceil(nx/8), ceil(ny/8), ceil(nz/8)];
        phi_gpu = gpuArray(phi);
        g_gpu = gpuArray(g);  */

// in simulation loop
    /*  phi_gpu = feval(k, phi_gpu, g_gpu, nx, ny, nz, gpoints);
        phi = gather(phi_gpu);  */

// updatePhi
    /* (double *phi, double *g, int nx, int ny, int nz, int gpoints) */

__global__ void momCollision(
        double *rho, 
        double *ux, 
        double *uy, 
        double *uz,  
        double *ffx,   
        double *ffy,    
        double *ffz,    
        double *f,         
        double *pxx,     
        double *pyy,        
        double *pzz,        
        double *pxy,        
        double *pxz,         
        double *pyz,        
        double *cix,         
        double *ciy,          
        double *ciz,          
        double *w,          
        double cssq,          
        int nx,              
        int ny,               
        int nz,              
        int fpoints           
    ) {
    
    // start at 1st index to match matlab 2nd index
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i > 0 && i < nx-1 && j > 0 && j < ny-1 && k > 0 && k < nz-1) { 

        ux[i][j][k] = (
            (f[i][j][k][1] + f[i][j][k][15] + f[i][j][k][9] + f[i][j][k][7] + f[i][j][k][13]) -
            (f[i][j][k][2] + f[i][j][k][10] + f[i][j][k][16] + f[i][j][k][14] + f[i][j][k][7])
        ) / rho[i][j][k] +
        ffx[i][j][k] * 0.5 / rho[i][j][k];
    
        uy[i][j][k] = (
            (f[i][j][k][3] + f[i][j][k][7] + f[i][j][k][14] + f[i][j][k][17] + f[i][j][k][11]) -
            (f[i][j][k][4] + f[i][j][k][13] + f[i][j][k][8] + f[i][j][k][12] + f[i][j][k][18])
        ) / rho[i][j][k] +
        ffy[i][j][k] * 0.5 / rho[i][j][k];
    
        uz[i][j][k] = (
            (f[i][j][k][6] + f[i][j][k][15] + f[i][j][k][10] + f[i][j][k][17] + f[i][j][k][12]) -
            (f[i][j][k][5] + f[i][j][k][9] + f[i][j][k][16] + f[i][j][k][11] + f[i][j][k][18])
        ) / rho[i][j][k] +
        ffz[i][j][k] * 0.5 / rho[i][j][k];
    
        uu[i][j][k] = 0.5 * (ux[i][j][k] * ux[i][j][k] + uy[i][j][k] * uy[i][j][k] + uz[i][j][k] * uz[i][j][k]) / cssq;

        for (int l = 0; l < fpoints; l++) {
            rho[i][j][k] += f[i][j][k][l];
        }

        for (int l = 0; l < fpoints; ++l) {

            double udotc = (ux[i][j][k] * cix[l] + uy[i][j][k] * ciy[l] + uz[i][j][k] * ciz[l]) / cssq;
        
            double HeF = (
                (w[l] * (rho[i][j][k] + rho[i][j][k] * (udotc + 0.5 * udotc * udotc - uu[i][j][k]))) *
                ((cix[l] - ux[i][j][k]) * ffx[i][j][k] +
                 (ciy[l] - uy[i][j][k]) * ffy[i][j][k] +
                 (ciz[l] - uz[i][j][k]) * ffz[i][j][k])
            ) / (rho[i][j][k] * cssq);
        
            double feq = w[l] * (rho[i][j][k] + rho[i][j][k] * (udotc + 0.5 * udotc * udotc - uu[i][j][k])) - 0.5 * HeF;
        
            fneq[l] = f[i][j][k][l] - feq;
        }

        pxx[i][j][k] = fneq[1] + fneq[2] + fneq[7] + fneq[8] + fneq[9] + fneq[10] + fneq[13] + fneq[14] + fneq[15] + fneq[16];
        pyy[i][j][k] = fneq[3] + fneq[4] + fneq[7] + fneq[8] + fneq[11] + fneq[12] + fneq[13] + fneq[14] + fneq[17] + fneq[18];
        pzz[i][j][k] = fneq[5] + fneq[6] + fneq[9] + fneq[10] + fneq[11] + fneq[12] + fneq[15] + fneq[16] + fneq[17] + fneq[18];
        pxy[i][j][k] = fneq[7] + fneq[8] - fneq[13] - fneq[14];
        pxz[i][j][k] = fneq[9] + fneq[10] - fneq[15] - fneq[16];
        pyz[i][j][k] = fneq[11] + fneq[12] - fneq[17] - fneq[18];

    }

    /* if (i <= nx-1 && j <= ny-1 && k <= nz-1) {

        // COLLISION
            uu = 0.5 * (ux(i,j,k).^2 + uy(i,j,k).^2 + uz(i,j,k).^2) / cssq;
            for l = 1:fpoints
                udotc = (ux(i,j,k) * cix(l) + uy(i,j,k) * ciy(l) + uz(i,j,k) * ciz(l)) / cssq;
                feq = w(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5.*udotc.^2 - uu));
                HeF = 0.5 .* (w(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5.*udotc.^2 - uu))) .* ...
                    ((cix(l) - ux(i,j,k)) .* ffx(i,j,k) + ...
                     (ciy(l) - uy(i,j,k)) .* ffy(i,j,k) + ...
                     (ciz(l) - uz(i,j,k)) .* ffz(i,j,k) ...
                    ) ./ (rho(i,j,k) .* cssq);
                fneq = (cix(l) .* cix(l) - cssq) * pxx(i,j,k) + ...
                       (ciy(l) .* ciy(l) - cssq) * pyy(i,j,k) + ...
                       (ciz(l) .* ciz(l) - cssq) * pzz(i,j,k) + ...
                       2 * cix(l) .* ciy(l) .* pxy(i,j,k) + ...
                       2 * cix(l) .* ciz(l) .* pxz(i,j,k) + ...
                       2 * ciy(l) .* ciz(l) .* pyz(i,j,k);
                f(i+cix(l),j+ciy(l),k+ciz(l),l) = feq + (1-omega) * (w(l) / (2*cssq^2)) * fneq + HeF;
            end
            for l = 1:gpoints
                udotc = (ux(i,j,k) * cix(l) + uy(i,j,k) * ciy(l) + uz(i,j,k) * ciz(l)) / cssq;
                feq = w_g(l) .* phi(i,j,k) .* (1 + udotc);
                Hi = sharp_c .* phi(i,j,k) .* (1 - phi(i,j,k)) .* (cix(l) .* normx(i,j,k) + ciy(l) .* normy(i,j,k) + ciz(l) .* normz(i,j,k)); 
                g(i,j,k,l) = feq + w_g(l) .* Hi;
            end  

    } */
}
