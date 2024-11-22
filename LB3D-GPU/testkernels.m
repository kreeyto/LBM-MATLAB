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
phi_gpu = gpuArray(phi);
g_gpu = gpuArray(g);

%% Loop de Simulação

for t = 1:nsteps

    phi_gpu = feval(k, phi_gpu, g_gpu, nx, ny, nz, gpoints);    
    phi = gather(phi_gpu);

end




