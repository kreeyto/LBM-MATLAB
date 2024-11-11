%% D3Q19 
clc; clearvars; close all

%% Parâmetros Gerais

% campo de velocidade do campo de fase
pf = "D3Q15";

slicebool = 1;
nlinks = 19;
tau = 0.6;
cssq = 1/3;
omega = 1/tau;
sharp_c = 1;
sigma = 0.024;

[nx, ny, nz] = deal(50);
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

[rho, u, v, w, ...
 fi, normx, normy, normz, ...
 curvature, indicator, ...
 ffx, ffy, ffz, ...
 mod_grad, isfluid] = deal(zeros(nx,ny,nz));

[pxx, pyy, pzz, ...
 pxy, pxz, pyz] = deal(ones(nx,ny,nz));

p = zeros(1,fpoints);
p_g = zeros(1,gpoints);

fneq = zeros(fpoints,1,1); 
isfluid(2:nx-1,2:ny-1,2:nz-1) = 1;
rho(:,:,:) = 1;

%% Propriedades do Modelo

p(1) = 1/3;
p(2:6) = 1/18;
p(7:19) = 1/36;

if pf == "D3Q19"
    p_g = p;
elseif pf == "D3Q15"
    p_g(1) = 2/9;
    p_g(2:7) = 1/9;
    p_g(8:15) = 1/72;
elseif pf == "D3Q7"
    p_g(1) = 1/4;
    p_g(2:7) = 1/8;
end

ex = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0];
ey = [0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1];
ez = [0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1];

%% Cálculo da Função de Distribuição em Função da Distância Radial

for i = 2:nx-1
    for j = 2:ny-1
        for k = 2:nz-1
            Ri = sqrt( ...
                    (i-(nx/2))^2 + ...
                    (j-(ny/2))^2 + ...
                    (k-(nz/2))^2 ...
                );
            fi(i,j,k) = 0.5 + 0.5 * tanh(2*(20-Ri)/3);
        end
    end
end

%% Inicialização de Funções de Distribuição

for i = 1:fpoints
    f(:,:,:,i) = p(i) * rho(:,:,:);
end

for i = 1:gpoints
    g(:,:,:,i) = p_g(i) * fi(:,:,:);
end

%% Visualização

if slicebool == 3
    hVol = volshow(fi, 'RenderingStyle', 'Isosurface');
    viewer = hVol.Parent;
    hFig = viewer.Parent;
end

%% pre-alloc

uindices1 = [2,16,10,8,14]; uindices2 = [3,11,17,15,8];
vindices1 = [4,8,15,18,12]; vindices2 = [5,14,9,13,19];
windices1 = [7,16,11,18,13]; windices2 = [6,10,17,12,19];

%% Loop de Simulação

for t = 1:nsteps

    mask = (isfluid == 1);

    % calculo do campo de fase
    fisums = sum(g,4);
    fi(mask) = fisums(mask);

    % gradientes
    [grad_fix, grad_fiy, grad_fiz] = deal(zeros(nx,ny,nz));
    for l = 1:fpoints
        circfi = circshift(fi,[-ex(l),-ey(l),-ez(l)]);
        grad_fix(mask) = grad_fix(mask) + 3 * p(l) .* ex(l) .* circfi(mask);
        grad_fiy(mask) = grad_fiy(mask) + 3 * p(l) .* ey(l) .* circfi(mask);
        grad_fiz(mask) = grad_fiz(mask) + 3 * p(l) .* ez(l) .* circfi(mask);
    end
    mod_grad(mask) = sqrt(grad_fix(mask).^2 + grad_fiy(mask).^2 + grad_fiz(mask).^2);
    normx(mask) = grad_fix(mask) ./ (mod_grad(mask) + 1e-9);
    normy(mask) = grad_fiy(mask) ./ (mod_grad(mask) + 1e-9);
    normz(mask) = grad_fiz(mask) ./ (mod_grad(mask) + 1e-9);
    indicator(mask) = mod_grad(mask);

    % curvatura
    curvature(mask) = 0;
    for l = 1:fpoints
        circnormx = circshift(normx,[-ex(l),-ey(l),-ez(l)]);
        circnormy = circshift(normy,[-ex(l),-ey(l),-ez(l)]);
        circnormz = circshift(normz,[-ex(l),-ey(l),-ez(l)]);
        curvature(mask) = curvature(mask) - 3 * p(l) .* ( ...
            ex(l) .* circnormx(mask) + ...
            ey(l) .* circnormy(mask) + ...
            ez(l) .* circnormz(mask) ...
        );
    end
    ffx(mask) = sigma .* curvature(mask) .* normx(mask) .* indicator(mask);
    ffy(mask) = sigma .* curvature(mask) .* normy(mask) .* indicator(mask);
    ffz(mask) = sigma .* curvature(mask) .* normz(mask) .* indicator(mask);

    % momentos
    usums = sum(f(:,:,:,uindices1),4) - sum(f(:,:,:,uindices2),4);
    vsums = sum(f(:,:,:,vindices1),4) - sum(f(:,:,:,vindices2),4);
    wsums = sum(f(:,:,:,windices1),4) - sum(f(:,:,:,windices2),4);
    u(mask) = usums(mask) ./ rho(mask) + ffx(mask) * 0.5 ./ rho(mask);
    v(mask) = vsums(mask) ./ rho(mask) + ffy(mask) * 0.5 ./ rho(mask);
    w(mask) = wsums(mask) ./ rho(mask) + ffz(mask) * 0.5 ./ rho(mask);
    uu = 0.5 * (u(mask).^2 + v(mask).^2 + w(mask).^2) / cssq;

    rhosums = sum(f,4);
    rho(mask) = rhosums(mask);

    disp(['Passo de tempo: ', num2str(t)]);
end
