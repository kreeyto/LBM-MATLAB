% D3Q19 
clc; clearvars; close all
%% Parâmetros Gerais
% campo de velocidade do campo de fase
pf = "D3Q15";
actingForces = 1;

slicebool = 1;
tau = 0.505;
cssq = 1/3;
omega = 1/tau;
sharp_c = 0.15*3;
sigma = 0.1;

radius = 20;
res = 1;

[nx, ny, nz] = deal(150*res);
nsteps = 10000; 

fpoints = 19; 
if pf == "D3Q19"
    gpoints = 19;
elseif pf == "D3Q15"
    gpoints = 15;
elseif pf == "D3Q7"
    gpoints = 7;
end
f = zeros(nx,ny,nz,fpoints); 
g = zeros(nx,ny,nz,gpoints); 

%% Matrizes e Variáveis

[rho, ux, uy, uz, ...
 phi, normx, normy, normz, ...
 curvature, indicator, ...
 ffx, ffy, ffz, ...
 mod_grad, isfluid] = deal(zeros(nx,ny,nz));

[pxx, pyy, pzz, ...
 pxy, pxz, pyz] = deal(ones(nx,ny,nz));

w = zeros(1,fpoints);
w_g = zeros(1,gpoints);

fneq = zeros(fpoints,1,1); 

isfluid(2:nx-1,2:ny-1,2:nz-1) = 1;
rho(:,:,:) = 1;

%% Propriedades do Modelo

w(1) = 1/3;
w(2:6) = 1/18;
w(7:19) = 1/36;

if pf == "D3Q19"
    w_g = w;
elseif pf == "D3Q15"
    w_g(1) = 2/9;
    w_g(2:7) = 1/9;
    w_g(8:15) = 1/72;
elseif pf == "D3Q7"
    w_g(1) = 1/4;
    w_g(2:7) = 1/8;
end

% opp = [1, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16, 19, 18];

cix = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0];
ciy = [0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1];
ciz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1];

%% Cálculo da Função de Distribuição em Função da Distância Radial

nx2 = nx/2; ny2 = ny/2; nz2 = nz/2; 
for i = 2:nx-1
    for j = 2:ny-1
        for k = 2:nz-1
            if actingForces == 1
                Ri = sqrt((i-nx2)^2/2.^2 + (j-ny2)^2 + (k-nz2)^2);
            else
                Ri = sqrt((i-nx2)^2 + (j-ny2)^2 + (k-nz2)^2);
            end
            % phi = 0.5 + 0.5 * tanh((2*gamma(x))/W) 
            % with W being the interface width 
            % and gamma(x) the coordinate perpendicular to the interface
            phi(i,j,k) = 0.5 + 0.5 * tanh(2*(radius*res-Ri)/(3*res));
        end
    end
end

%% Inicialização de Funções de Distribuição

for i = 1:fpoints
    f(:,:,:,i) = w(i) * rho(:,:,:);
end

for i = 1:gpoints
    g(:,:,:,i) = w_g(i) * phi(:,:,:);
end

%% Visualização

if slicebool == 3
    hVol = volshow(phi, 'RenderingStyle', 'Isosurface');
    viewer = hVol.Parent;
    hFig = viewer.Parent;
end

%% Kernels

k = parallel.gpu.CUDAKernel('myKernel.ptx', 'myKernel.cu', 'updatePhi');
k.ThreadBlockSize = [8, 8, 8];
k.GridSize = [ceil(nx/8), ceil(ny/8), ceil(nz/8)];

% Alocação de variáveis na GPU
rho_gpu = gpuArray(rho);
ux_gpu = gpuArray(ux);
uy_gpu = gpuArray(uy);
uz_gpu = gpuArray(uz);
ffx_gpu = gpuArray(ffx);
ffy_gpu = gpuArray(ffy);
ffz_gpu = gpuArray(ffz);
f_gpu = gpuArray(f);
pxx_gpu = gpuArray(pxx);
pyy_gpu = gpuArray(pyy);
pzz_gpu = gpuArray(pzz);
pxy_gpu = gpuArray(pxy);
pxz_gpu = gpuArray(pxz);
pyz_gpu = gpuArray(pyz);
cix_gpu = gpuArray(cix);
ciy_gpu = gpuArray(ciy);
ciz_gpu = gpuArray(ciz);
w_gpu = gpuArray(w);

% Constantes na GPU
cssq_gpu = gpuArray(cssq);
fpoints_gpu = gpuArray(fpoints);
nx_gpu = gpuArray(nx);
ny_gpu = gpuArray(ny);
nz_gpu = gpuArray(nz);

%% Loop de Simulação

for t = 1:nsteps

    % phase field calc
    phi(ix,iy,iz) = sum(g(ix,iy,iz,:),4);

    % normal and arrays
    [grad_fix, grad_fiy, grad_fiz] = deal(0);
    for l = 1:fpoints
        grad_fix = grad_fix + 3 * w(l) .* cix(l) .* ((phi((ix)+cix(l),(iy)+ciy(l),(iz)+ciz(l))));
        grad_fiy = grad_fiy + 3 * w(l) .* ciy(l) .* ((phi((ix)+cix(l),(iy)+ciy(l),(iz)+ciz(l))));
        grad_fiz = grad_fiz + 3 * w(l) .* ciz(l) .* ((phi((ix)+cix(l),(iy)+ciy(l),(iz)+ciz(l))));
    end
    mod_grad(ix,iy,iz) = sqrt(grad_fix.^2 + grad_fiy.^2 + grad_fiz.^2);
    normx(ix,iy,iz) = grad_fix ./ (mod_grad(ix,iy,iz) + 1e-9);
    normy(ix,iy,iz) = grad_fiy ./ (mod_grad(ix,iy,iz) + 1e-9);
    normz(ix,iy,iz) = grad_fiz ./ (mod_grad(ix,iy,iz) + 1e-9);
    indicator(ix,iy,iz) = sqrt(grad_fix.^2 + grad_fiy.^2 + grad_fiz.^2);

    % curvature
    curvature(ix,iy,iz) = 0;
    for l = 1:fpoints
        curvature(ix,iy,iz) = curvature(ix,iy,iz) - 3 .* w(l) .* ...
        (cix(l) .* (normx((ix)+cix(l),(iy)+ciy(l),(iz)+ciz(l))) + ...
         ciy(l) .* (normy((ix)+cix(l),(iy)+ciy(l),(iz)+ciz(l))) + ...
         ciz(l) .* (normz((ix)+cix(l),(iy)+ciy(l),(iz)+ciz(l))) ...
        );
    end
    ffx(ix,iy,iz) = sigma .* curvature(ix,iy,iz) .* normx(ix,iy,iz) .* indicator(ix,iy,iz);
    ffy(ix,iy,iz) = sigma .* curvature(ix,iy,iz) .* normy(ix,iy,iz) .* indicator(ix,iy,iz);
    ffz(ix,iy,iz) = sigma .* curvature(ix,iy,iz) .* normz(ix,iy,iz) .* indicator(ix,iy,iz);

    % EXEMPLO FEVAL
    [rho_gpu, ux_gpu, uy_gpu, uz_gpu, pxx_gpu, pyy_gpu, pzz_gpu, ...
     pxy_gpu, pxz_gpu, pyz_gpu] = feval(kCollision, ...
        rho_gpu, ux_gpu, uy_gpu, uz_gpu, ...
        ffx_gpu, ffy_gpu, ffz_gpu, f_gpu, ...
        pxx_gpu, pyy_gpu, pzz_gpu, pxy_gpu, ...
        pxz_gpu, pyz_gpu, cix_gpu, ciy_gpu, ciz_gpu, ...
        w_gpu, cssq_gpu, nx_gpu, ny_gpu, nz_gpu, fpoints_gpu);

    rho = gather(rho_gpu);
    ux = gather(ux_gpu);
    uy = gather(uy_gpu);
    uz = gather(uz_gpu);
    pxx = gather(pxx_gpu);
    pyy = gather(pyy_gpu);
    pzz = gather(pzz_gpu);
    pxy = gather(pxy_gpu);
    pxz = gather(pxz_gpu);
    pyz = gather(pyz_gpu);

end




