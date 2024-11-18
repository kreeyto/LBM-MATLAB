%% D3Q19 
clc; clearvars; close all
%% Parâmetros Gerais
tic
% campo de velocidade do campo de fase
pf = "D3Q15";

slicebool = 1;
tau = 1;
cssq = 1/3;
omega = 1/tau;
sharp_c = 0.1;
sigma = 0.024;

[nx, ny, nz] = deal(65);
nsteps = 10000; 

fpoints = 19; 
if pf == "D3Q19"
    gpoints = 19;
elseif pf == "D3Q15"
    gpoints = 15;
elseif pf == "D3Q7"
    gpoints = 7;
end
f = gpuArray.zeros(nx,ny,nz,fpoints); 
g = gpuArray.zeros(nx,ny,nz,gpoints); 

%% Matrizes e Variáveis

[rho, ux, uy, uz, ...
 phi, normx, normy, normz, ...
 curvature, indicator, ...
 ffx, ffy, ffz, ...
 mod_grad, isfluid] = deal(gpuArray.zeros(nx,ny,nz));

[pxx, pyy, pzz, ...
 pxy, pxz, pyz] = deal(gpuArray.ones(nx,ny,nz));

w = gpuArray.zeros(1,fpoints);
w_g = gpuArray.zeros(1,gpoints);

fneq = gpuArray.zeros(fpoints,1,1); 

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

cix = gpuArray([0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0]);
ciy = gpuArray([0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1]);
ciz = gpuArray([0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1]);

%% Cálculo da Função de Distribuição em Função da Distância Radial

for i = 2:nx-1
    for j = 2:ny-1
        for k = 2:nz-1
            Ri = sqrt( ...
                    (i-(nx/2))^2 + ...
                    (j-(ny/2))^2 + ...
                    (k-(nz/2))^2 ...
                );
            phi(i,j,k) = 0.5 + 0.5 * tanh(2*(20-Ri)/3);
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

for i = 2:nx-1
    for j = 2:ny-1
        for k = 2:nz-1
            [grad_fix, grad_fiy, grad_fiz] = deal(0);
            for l = 1:fpoints
                grad_fix = grad_fix + 3 * w(l) .* cix(l) .* ((phi(i+cix(l),j+ciy(l),k+ciz(l))));
                grad_fiy = grad_fiy + 3 * w(l) .* ciy(l) .* ((phi(i+cix(l),j+ciy(l),k+ciz(l))));
                grad_fiz = grad_fiz + 3 * w(l) .* ciz(l) .* ((phi(i+cix(l),j+ciy(l),k+ciz(l))));
            end
            mod_grad(i,j,k) = sqrt(grad_fix.^2 + grad_fiy.^2 + grad_fiz.^2);
            normx(i,j,k) = grad_fix ./ (mod_grad(i,j,k) + 1e-9);
            normy(i,j,k) = grad_fiy ./ (mod_grad(i,j,k) + 1e-9);
            normz(i,j,k) = grad_fiz ./ (mod_grad(i,j,k) + 1e-9);
            indicator(i,j,k) = sqrt(grad_fix.^2 + grad_fiy.^2 + grad_fiz.^2);
        end
    end
end