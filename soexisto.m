%% D3Q19 - Código Vetorizado
clc; clearvars; close all

%% Parâmetros Gerais

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
gpoints = 19;
if pf == "D3Q15"
    gpoints = 15;
elseif pf == "D3Q7"
    gpoints = 7;
end
f = zeros(nx, ny, nz, fpoints);
g = zeros(nx, ny, nz, gpoints);

%% Matrizes e Variáveis

[rho, u, v, w, ...
 fi, normx, normy, normz, ...
 curvature, indicator, ...
 ffx, ffy, ffz, ...
 mod_grad, isfluid] = deal(zeros(nx, ny, nz));

[pxx, pyy, pzz, ...
 pxy, pxz, pyz] = deal(ones(nx, ny, nz));

p = zeros(1, fpoints);
p_g = zeros(1, gpoints);

fneq = zeros(fpoints, 1, 1);
isfluid(2:nx-1, 2:ny-1, 2:nz-1) = 1;
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

[i_grid, j_grid, k_grid] = ndgrid(2:nx-1, 2:ny-1, 2:nz-1);
Ri = sqrt((i_grid - (nx/2)).^2 + (j_grid - (ny/2)).^2 + (k_grid - (nz/2)).^2);
fi(2:nx-1, 2:ny-1, 2:nz-1) = 0.5 + 0.5 * tanh(2 * (20 - Ri) / 3);

%% Inicialização de Funções de Distribuição

f(:,:,:,1:fpoints) = p .* rho;
g(:,:,:,1:gpoints) = p_g .* fi;

%% Loop de Simulação Vetorizado

for t = 1:nsteps
    % Campo de fase vetorizado
    fi(isfluid == 1) = sum(g(isfluid == 1, :, :, :), 4);

    % Gradientes vetorizados
    [grad_fix, grad_fiy, grad_fiz] = deal(zeros(size(fi)));
    for l = 1:fpoints
        shifted_fi = circshift(fi, [-ex(l), -ey(l), -ez(l)]);
        grad_fix = grad_fix + 3 * p(l) * ex(l) .* shifted_fi;
        grad_fiy = grad_fiy + 3 * p(l) * ey(l) .* shifted_fi;
        grad_fiz = grad_fiz + 3 * p(l) * ez(l) .* shifted_fi;
    end
    mod_grad(isfluid == 1) = sqrt(grad_fix(isfluid == 1).^2 + grad_fiy(isfluid == 1).^2 + grad_fiz(isfluid == 1).^2);
    normx(isfluid == 1) = grad_fix(isfluid == 1) ./ (mod_grad(isfluid == 1) + 1e-9);
    normy(isfluid == 1) = grad_fiy(isfluid == 1) ./ (mod_grad(isfluid == 1) + 1e-9);
    normz(isfluid == 1) = grad_fiz(isfluid == 1) ./ (mod_grad(isfluid == 1) + 1e-9);
    indicator(isfluid == 1) = mod_grad(isfluid == 1);

    % Curvatura vetorizada
    curvature(isfluid == 1) = -sum(3 * p .* ...
        (ex .* circshift(normx, [-ex; -ey; -ez]) + ...
         ey .* circshift(normy, [-ex; -ey; -ez]) + ...
         ez .* circshift(normz, [-ex; -ey; -ez])), 1);
    
    % Atualização de forças
    ffx(isfluid == 1) = sigma * curvature(isfluid == 1) .* normx(isfluid == 1) .* indicator(isfluid == 1);
    ffy(isfluid == 1) = sigma * curvature(isfluid == 1) .* normy(isfluid == 1) .* indicator(isfluid == 1);
    ffz(isfluid == 1) = sigma * curvature(isfluid == 1) .* normz(isfluid == 1) .* indicator(isfluid == 1);

    % Visualização a cada passo
    if mod(t, 1) == 0
        if slicebool == 1
            clf; x = 1:nx; y = 1:ny; z = 1:nz;
            h = slice(x, y, z, fi, [], ny/2, []);
            shading interp; colorbar;
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title(['t = ', num2str(t)]);
            view(3); drawnow;
        elseif slicebool == 2
            clf; x = 1:nx; y = 1:ny; z = 1:nz;
            surfpatch = patch(isosurface(x, y, z, fi));
            set(surfpatch, 'FaceColor', 'red', 'EdgeColor', 'none');
            xlabel('X'); ylabel('Y'); zlabel('Z');
            axis equal;
            camlight; lighting phong;
            title(['t = ', num2str(t)]);
            view(3); drawnow;
        else
            hVol.Data = fi;
            drawnow;
        end
    end

    disp(['Passo de tempo: ', num2str(t)]);
end
