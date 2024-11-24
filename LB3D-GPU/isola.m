clc; clearvars; close all

[nx, ny, nz] = deal(128);

k = parallel.gpu.CUDAKernel('myKernel.ptx', 'myKernel.cu', 'momCollision');
k.ThreadBlockSize = [8, 8, 8];
k.GridSize = [ceil(nx/8), ceil(ny/8), ceil(nz/8)];

rho = rand(nx, ny, nz, 'gpuArray'); 
ux = zeros(nx, ny, nz, 'gpuArray');
uy = zeros(nx, ny, nz, 'gpuArray');
uz = zeros(nx, ny, nz, 'gpuArray');
ffx = rand(nx, ny, nz, 'gpuArray');
ffy = rand(nx, ny, nz, 'gpuArray');
ffz = rand(nx, ny, nz, 'gpuArray');
f = rand(nx, ny, nz, 19, 'gpuArray');  

[ux, uy, uz] = feval(k, rho, ux, uy, uz, ffx, ffy, ffz, f, nx, ny, nz);
