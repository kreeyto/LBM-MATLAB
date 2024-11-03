%% D3Q19 
clc; clearvars; close all

%% Parâmetros Gerais

op = 2;

slicebool = 2;
nlinks = 19;
tau = 0.8;
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
    rho = sum(f, 4);
    if op == 1

        [u, v, w] = deal(zeros(nx,ny,nz));
        for l = 1:19
            u = u + f(:,:,:,l) * ex(l);
            v = v + f(:,:,:,l) * ey(l);
            w = w + f(:,:,:,l) * ez(l);
        end

    elseif op == 2

        u = (f(:,:,:,2) + f(:,:,:,16) + f(:,:,:,10) + f(:,:,:,8) + f(:,:,:,14)) - (f(:,:,:,3) + f(:,:,:,11) + f(:,:,:,17) + f(:,:,:,15) + f(:,:,:,9));
        v = (f(:,:,:,4) + f(:,:,:,8) + f(:,:,:,15) + f(:,:,:,18) + f(:,:,:,12)) - (f(:,:,:,5) + f(:,:,:,14) + f(:,:,:,9) + f(:,:,:,13) + f(:,:,:,19));
        w = (f(:,:,:,7) + f(:,:,:,16) + f(:,:,:,11) + f(:,:,:,18) + f(:,:,:,13)) - (f(:,:,:,6) + f(:,:,:,10) + f(:,:,:,17) + f(:,:,:,12) + f(:,:,:,19));

    end
    u = u ./ rho + ffx * 0.5 ./ rho;
    v = v ./ rho + ffy * 0.5 ./ rho;
    w = w ./ rho + ffz * 0.5 ./ rho;
 
    % non existant
    if op == 2

        uu = 0.5 * (u.^2 + v.^2 + w.^2) / cssq;

    end
    
    for l = 1:19

        if op == 1

            [pxx, pyy, pzz, pxy, pxz, pyz] = deal(zeros(nx,ny,nz));
            
            udotc = ex(l) * u + ey(l) * v + ez(l) * w;
            feq = p(l) .* rho .* (1 + udotc / cssq + (udotc.^2) / (2 * cssq^2) - (u.^2 + v.^2 + w.^2) / (2 * cssq));
            fneq = f(:,:,:,l) - feq;

            pxx = pxx + ex(l)^2 * fneq;
            pyy = pyy + ey(l)^2 * fneq;
            pzz = pzz + ez(l)^2 * fneq;
            pxy = pxy + ex(l) * ey(l) * fneq;
            pxz = pxz + ex(l) * ez(l) * fneq;
            pyz = pyz + ey(l) * ez(l) * fneq;

        elseif op == 2

            udotc = (u * ex(l) + v * ey(l) + w * ez(l)) / cssq;
            HeF = (p(l) * (rho + rho .* (udotc + 0.5 .* udotc.^2 - uu))) ...
                    .* ((ex(l) - u) .* ffx + ...
                        (ey(l) - v) .* ffy + ...
                        (ez(l) - w) .* ffz ...
                       ) ./ (rho .* cssq);
            feq = p(l) * (rho + rho .* (udotc + 0.5 .* udotc.^2 - uu)) - 0.5 .* HeF;
            fneq = f(:,:,:,l) - feq;
 
            pxx = sum(fneq([2, 3, 8, 9, 10, 11, 14, 15, 16, 17]));
            pyy = sum(fneq([4, 5, 8, 9, 12, 13, 14, 15, 18, 19]));
            pzz = sum(fneq([6, 7, 10, 11, 12, 13, 16, 17, 18, 19]));
            pxy = sum(fneq([8, 9])) - sum(fneq([14, 15]));
            pxz = sum(fneq([10, 11])) - sum(fneq([16, 17]));
            pyz = sum(fneq([12, 13])) - sum(fneq([18, 19]));

        end

    end

    % Colisão

    if op == 2

        uu = 0.5 * (u.^2 + v.^2 + w.^2) / cssq;

    end
    
    for l = 1:19

        if op == 1

            udotc = ex(l) * u + ey(l) * v + ez(l) * w;
            feq = p(l) .* rho .* (1 + udotc / cssq + (udotc.^2) / (2 * cssq^2) - (u.^2 + v.^2 + w.^2) / (2 * cssq));
            GuoF = (1 - 0.5 * omega) * p(l) .* ((ex(l) - u) .* ffx + (ey(l) - v) .* ffy + (ez(l) - w) .* ffz) / cssq;
            f(:,:,:,l) = f(:,:,:,l) - omega * (f(:,:,:,l) - feq) + GuoF;

        elseif op == 2

            udotc = (u * ex(l) + v * ey(l) + w * ez(l)) / cssq;
            feq = p(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5 .* udotc.^2 - uu));
            HeF = 0.5 * (p(l) * (rho + rho .* (udotc + 0.5 .* udotc.^2 - uu))) ...
                    .* ((ex(l) - u) .* ffx + ...
                        (ey(l) - v) .* ffy + ...
                        (ez(l) - w) .* ffz ...
                       ) ./ (rho .* cssq);
            fneq = (ex(l) .* ex(l) - cssq) * pxx + ...
                   (ey(l) .* ey(l) - cssq) * pyy + ...
                   (ez(l) .* ez(l) - cssq) * pzz + ...
                    2 * ex(l) .* ey(l) .* pxy + ...
                    2 * ex(l) .* ez(l) .* pxz + ...
                    2 * ey(l) .* ez(l) .* pyz;
            fhold = circshift(f,[-ex(l),-ey(l),-ez(l),0]);
            fhold = feq + (1-omega) * (p(l) / (2*cssq^2)) * fneq + HeF;

        end

    end
    for l = 1:gpoints
        
        if op == 1

            udotc = ex(l) * u + ey(l) * v + ez(l) * w;
            feq = p_g(l) .* fi .* (1 + udotc / cssq);
            Hi = sharp_c .* fi .* (1 - fi) .* (ex(l) * normx + ey(l) * normy + ez(l) * normz);
            g(:,:,:,l) = (1 - omega) * (g(:,:,:,l) - feq) + feq + p_g(l) .* Hi;

        elseif op == 2

            udotc = (u * ex(l) + v * ey(l) + w * ez(l)) / cssq;
            feq = p_g(l) .* fi .* (1 + udotc);
            Hi = sharp_c .* fi .* (1 - fi) .* (ex(l) .* normx + ey(l) .* normy + ez(l) .* normz); 
            g(:,:,:,l) = feq + p_g(l) .* Hi;

        end
        
    end

    % Streaming

    if op == 1

        for l = 1:19
            f(:,:,:,l) = circshift(f(:,:,:,l),[ex(l),ey(l),ez(l)]);
        end
        for l = 1:gpoints
            g(:,:,:,l) = circshift(g(:,:,:,l),[ex(l),ey(l),ez(l)]);
        end

    elseif op == 2

        for l = 1:gpoints
            g(:,:,:,l) = circshift(g(:,:,:,l),[ex(l),ey(l),ez(l),0]);
        end

    end

    % Condições de contorno periódicas
    
    if op == 1

        f = apd(f);
        g = apd(g);

    elseif op == 2

        for i = 1:nx
            for j = 1:ny
                for k = 1:nz
                    if isfluid(i,j,k) == 1
                        for l = 1:19
                            if (i + ex(l) > 0 && j + ey(l) > 0 && k + ez(l) > 0)
                                f(i + ex(l), j + ey(l), k + ez(l), l) = rho(i,j,k) .* p(l); 
                            end
                        end
                        for l = 1:gpoints
                            if (i + ex(l) > 0 && j + ey(l) > 0 && k + ez(l) > 0)
                                g(i + ex(l), j + ey(l), k + ez(l), l) = fi(i,j,k) .* p_g(l);
                            end
                        end
                    end
                end
            end
        end

        fi(:, :, 1) = fi(:, :, 2);  
        fi(:, :, nz) = fi(:, :, nz-1); 
        fi(1, :, :) = fi(2, :, :); 
        fi(nx, :, :) = fi(nx-1, :, :); 
        fi(:, 1, :) = fi(:, 2, :); 
        fi(:, ny, :) = fi(:, ny-1, :); 

    end

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
