%% D3Q19 
clc; clearvars; close all

%% Parâmetros Gerais

HeGuo = 2;

slicebool = 2;
nlinks = 19;
tau = 0.8s;
cssq = 1/3;
omega = 1/tau;
sharp_c = 0.1;
sigma = 0.024;

[nx, ny, nz] = deal(50);
nsteps = 20000; 

gpoints = 15;
f = zeros(nx,ny,nz,19); 
g = zeros(nx,ny,nz,gpoints); 

%% Matrizes e Variáveis

[rho, u, v, w, ...
 fi, normx, normy, normz, ...
 curvature, indicator, ...
 ffx, ffy, ffz, ...
 mod_grad, isfluid] = deal(zeros(nx,ny,nz));

[pxx, pyy, pzz, ...
 pxy, pxz, pyz] = deal(ones(nx,ny,nz));

p = zeros(1,19);
p_g = zeros(1,gpoints);

fneq = zeros(19,1,1); 
isfluid(2:nx-1,2:ny-1,2:nz-1) = 1;
rho(:,:,:) = 1;

%% Propriedades do Modelo

p(1) = 1/3;
p(2:6) = 1/18;
p(7:19) = 1/36;

p_g(1) = 2/9;
p_g(2:7) = 1/9;
p_g(8:15) = 1/72;

opp = [1, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16, 19, 18];

ex = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0];
ey = [0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1];
ez = [0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1];

%% Cálculo da Função de Distribuição em Função da Distância Radial

% R0 = 10; 
for i = 1:nx
    for j = 1:ny
        for k = 1:nz
            Ri = sqrt((i - nx/2)^2 + (j - ny/2)^2 + (k - nz/2)^2);
            % fi(i, j, k) = 0.5 * (1 - tanh((Ri - R0) / (2 * sqrt(2) * sharp_c)));
            fi(i,j,k) = 0.5 + 0.5 * tanh(10*(20-Ri)/3);
        end
    end
end

%% Inicialização de Funções de Distribuição

for i = 1:19
    f(:,:,:,i) = p(i) * rho(:,:,:);
end

for i = 1:gpoints
    g(:,:,:,i) = p_g(i) * fi(:,:,:);
end

%% Loop de Simulação

for t = 1:nsteps

    % Campo de fase
    fi = sum(g, 4);
    rho = sum(f, 4);

    % Gradiente e normais
    [grad_fix, grad_fiy, grad_fiz] = deal(zeros(nx, ny, nz));
    for l = 1:19
        grad_fix = grad_fix + 3 * p(l) .* ex(l) .* circshift(fi,[-ex(l),-ey(l),-ez(l)]);
        grad_fiy = grad_fiy + 3 * p(l) .* ey(l) .* circshift(fi,[-ex(l),-ey(l),-ez(l)]);
        grad_fiz = grad_fiz + 3 * p(l) .* ez(l) .* circshift(fi,[-ex(l),-ey(l),-ez(l)]);
    end
    mod_grad = sqrt(grad_fix.^2 + grad_fiy.^2 + grad_fiz.^2) + 1e-8;
    normx = grad_fix ./ mod_grad;
    normy = grad_fiy ./ mod_grad;
    normz = grad_fiz ./ mod_grad;
    indicator = mod_grad;

    % Curvatura e forças de tensão superficial
    curvature = zeros(nx, ny, nz);
    for l = 1:19
        circnormx = circshift(normx,[-ex(l),-ey(l),-ez(l)]);
        circnormy = circshift(normy,[-ex(l),-ey(l),-ez(l)]);
        circnormz = circshift(normz,[-ex(l),-ey(l),-ez(l)]);
        curvature = curvature - 3 * p(l) .* (ex(l) .* circnormx + ey(l) .* circnormy + ez(l) .* circnormz);
    end
    ffx = sigma .* curvature .* normx .* indicator;
    ffy = sigma .* curvature .* normy .* indicator;
    ffz = sigma .* curvature .* normz .* indicator;

    % Momentos
    [u, v, w] = deal(zeros(nx,ny,nz));
    for l = 1:19

        % CHECAR VERACIDADE !!!
        u = u + f(:,:,:,l) * ex(l);
        v = v + f(:,:,:,l) * ey(l);
        w = w + f(:,:,:,l) * ez(l);

    end
    u = u ./ rho + ffx * 0.5 ./ rho;
    v = v ./ rho + ffy * 0.5 ./ rho;
    w = w ./ rho + ffz * 0.5 ./ rho;
     
    uu = 0.5 * (u.^2 + v.^2 + w.^2) / cssq;
    [pxx, pyy, pzz, pxy, pxz, pyz] = deal(zeros(nx,ny,nz));
    for l = 1:19
        udotc = (u * ex(l) + v * ey(l) + w * ez(l)) / cssq;
        feq = p(l) * (rho + rho .* (udotc + 0.5 .* udotc.^2 - uu));
        fneq = f(:,:,:,l) - feq;
        pxx = pxx + ex(l)^2 * fneq;
		pyy = pyy + ey(l)^2 * fneq;
		pzz = pzz + ez(l)^2 * fneq;
		pxy = pxy + ex(l) * ey(l) * fneq;
		pxz = pxz + ex(l) * ez(l) * fneq;
		pyz = pyz + ey(l) * ez(l) * fneq;
    end

    % Colisão
    uu = 0.5 * (u.^2 + v.^2 + w.^2) / cssq;
    for l = 1:19
        udotc = (u * ex(l) + v * ey(l) + w * ez(l)) / cssq;
        feq = p(l) * (rho + rho .* (udotc + 0.5 .* udotc.^2 - uu));
        GuoF = (1 - 0.5 * omega) * p(l) .* ( ( (ex(l)-u).*ffx + (ey(l)-v).*ffy + (ez(l)-w).*ffz )/cssq + ( ((ex(l)-u).*ex(l)).*ffx + ((ey(l)-v).*ey(l)).*ffy + ((ez(l)-w).*ez(l)).*ffz )/cssq^2 );
        f(:,:,:,l) = f(:,:,:,l) - omega * (f(:,:,:,l) - feq) + GuoF;
    end
    for l = 1:gpoints
        udotc = (u * ex(l) + v * ey(l) + w * ez(l)) / cssq;
        feq = p_g(l) .* fi .* (1 + udotc);
        Hi = sharp_c .* fi .* (1 - fi) .* (ex(l) * normx + ey(l) * normy + ez(l) * normz);
        g(:,:,:,l) = feq + p_g(l) .* Hi + (1 - omega) .* (g(:,:,:,l) - feq);
    end

    % Streaming de g
    for l = 1:gpoints
        f(:,:,:,l) = circshift(f(:,:,:,l),[ex(l),ey(l),ez(l)]);
    end
    for l = 1:gpoints
        g(:,:,:,l) = circshift(g(:,:,:,l),[ex(l),ey(l),ez(l)]);
    end

    % Condições de contorno periódicas
    f = apd(f);
    g = apd(g);

    % Visualização 
    if mod(t, 1) == 0
        if slicebool == 1
            hFig = figure(1); clf;
            x = 1:nx; y = 1:ny; z = 1:nz;
            h = slice(x, y, z, fi, [], ny/2, []); 
            shading interp; colorbar; 
            xlabel('X'); ylabel('Y'); zlabel('Z'); 
            title(['t = ', num2str(t)]);
            view(3); drawnow; 
        else
            hFig = figure(1); clf;
            x = 1:nx; y = 1:ny; z = 1:nz;
            surfpatch = patch(isosurface(x, y, z, fi, 0.5));
            set(surfpatch, 'FaceColor', 'red', 'EdgeColor', 'none'); 
            xlabel('X'); ylabel('Y'); zlabel('Z');
            axis([1 nx 1 ny 1 nz]);
            axis equal;
            camlight; lighting phong; 
            title(['t = ', num2str(t)]);
            view(3); drawnow;
        end
    end

    disp(['Passo de tempo: ', num2str(t)]);
end

% Função condições de contorno periódicas
function f = apd(f)
    f(1,:,:,:) = f(end-1,:,:,:);
    f(end,:,:,:) = f(2,:,:,:);
    f(:,1,:,:) = f(:,end-1,:,:);
    f(:,end,:,:) = f(:,2,:,:);
    f(:,:,1,:) = f(:,:,end-1,:);
    f(:,:,end,:) = f(:,:,2,:);
end
