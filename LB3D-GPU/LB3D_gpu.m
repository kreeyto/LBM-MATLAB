clc; clearvars; close all

% parameters

radius = 20;
tau = 0.505;
cssq = 1/3;
omega = 1/tau;
sharp_c = 0.15*3;
sigma = 0.1;

[nx, ny, nz] = deal(128);

% kernel

sz = 8;
timeLoop = parallel.gpu.CUDAKernel('gpuMomCollisionStreamBCS.ptx', 'gpuMomCollisionStreamBCS.cu', 'momCollisionStreamBCS');
timeLoop.ThreadBlockSize = [sz, sz, sz];
timeLoop.GridSize = [ceil(nx/sz), ceil(ny/sz), ceil(nz/sz)]; 

% remaining pars

nsteps = 10000; 

fpoints = 19; 
gpoints = 15;
f = zeros(nx,ny,nz,fpoints,'gpuArray'); 
g = zeros(nx,ny,nz,gpoints,'gpuArray'); 

% index pre-alloc

ix = 2:nx-1; iy = 2:ny-1; iz = 2:nz-1;

% arrays and variables

[rho, ux, uy, uz, ...
 ffx, ffy, ffz, phi, ...
 normx, normy, normz, ...
 mod_grad, curvature, indicator] = deal(zeros(nx,ny,nz,'gpuArray'));

[pxx, pyy, pzz, ...
 pxy, pxz, pyz] = deal(ones(nx,ny,nz,'gpuArray'));

w = zeros(1,fpoints,'gpuArray');
w_g = zeros(1,gpoints,'gpuArray');
rho(:,:,:) = 1;

% velocity set properties

w(1) = 1/3;
w(2:6) = 1/18;
w(7:19) = 1/36;

w_g(1) = 2/9;
w_g(2:7) = 1/9;
w_g(8:15) = 1/72;

cix = gpuArray([0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0]);
ciy = gpuArray([0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1]);
ciz = gpuArray([0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1]);

% phase field init

for i = ix
    for j = iy
        for k = iz
            Ri = sqrt((i-nx/2)^2/2.^2 + (j-ny/2)^2 + (k-nz/2)^2);
            phi(i,j,k) = 0.5 + 0.5 * tanh(2*(radius-Ri)/3);
        end
    end
end

% distribution function init

for i = 1:fpoints
    f(:,:,:,i) = w(i) * rho(:,:,:);
end

for i = 1:gpoints
    g(:,:,:,i) = w_g(i) * phi(:,:,:);
end

%% simulation loop

try
    feval(timeLoop, ...
        f, g, phi, rho, w, w_g, ...
        cix, ciy, ciz, ...
        mod_grad, normx, normy, normz, indicator, ...
        curvature, ffx, ffy, ffz, ...
        ux, uy, uz, ...
        pxx, pyy, pzz, pxy, pxz, pyz, ...
        nx, ny, nz, fpoints, gpoints, ...
        sigma, cssq, omega, sharp_c, nsteps);
    gpuDeviceSynchronize();
catch ME
    disp("erro na exec do kernel:");
    disp(ME.message);
end

%% iterate through phi saved values

% insert code

