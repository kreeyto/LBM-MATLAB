%% D3Q19 
clc; clearvars; close all
%% Parâmetros Gerais
% campo de velocidade do campo de fase
pf = "D3Q15";

slicebool = 1;
tau = 1;
cssq = 1/3;
omega = 1/tau;
sharp_c = 0.1;
sigma = 0.024;

[nx, ny, nz] = deal(64);
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
            Ri = sqrt((i-nx2)^2 + (j-ny2)^2 + (k-nz2)^2);
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

%% Visualização

if slicebool == 3
    hVol = volshow(phi, 'RenderingStyle', 'Isosurface');
    viewer = hVol.Parent;
    hFig = viewer.Parent;
end

%% Loop de Simulação

for t = 1:nsteps

    % Campo de fase
    phi(2:nx-1,2:ny-1,2:nz-1) = sum(g(2:nx-1,2:ny-1,2:nz-1,:),4);

    % Normal e arrays
    [grad_fix, grad_fiy, grad_fiz] = deal(0);
    for l = 1:fpoints
        grad_fix = grad_fix + 3 * w(l) .* cix(l) .* ((phi((2:nx-1)+cix(l),(2:ny-1)+ciy(l),(2:nz-1)+ciz(l))));
        grad_fiy = grad_fiy + 3 * w(l) .* ciy(l) .* ((phi((2:nx-1)+cix(l),(2:ny-1)+ciy(l),(2:nz-1)+ciz(l))));
        grad_fiz = grad_fiz + 3 * w(l) .* ciz(l) .* ((phi((2:nx-1)+cix(l),(2:ny-1)+ciy(l),(2:nz-1)+ciz(l))));
    end
    mod_grad(2:nx-1,2:ny-1,2:nz-1) = sqrt(grad_fix.^2 + grad_fiy.^2 + grad_fiz.^2);
    normx(2:nx-1,2:ny-1,2:nz-1) = grad_fix ./ (mod_grad(2:nx-1,2:ny-1,2:nz-1) + 1e-9);
    normy(2:nx-1,2:ny-1,2:nz-1) = grad_fiy ./ (mod_grad(2:nx-1,2:ny-1,2:nz-1) + 1e-9);
    normz(2:nx-1,2:ny-1,2:nz-1) = grad_fiz ./ (mod_grad(2:nx-1,2:ny-1,2:nz-1) + 1e-9);
    indicator(2:nx-1,2:ny-1,2:nz-1) = sqrt(grad_fix.^2 + grad_fiy.^2 + grad_fiz.^2);

    % Curvatura
    curvature(2:nx-1,2:ny-1,2:nz-1) = 0;
    for l = 1:fpoints
        curvature(2:nx-1,2:ny-1,2:nz-1) = curvature(2:nx-1,2:ny-1,2:nz-1) - 3 .* w(l) .* ...
        (cix(l) .* (normx((2:nx-1)+cix(l),(2:ny-1)+ciy(l),(2:nz-1)+ciz(l))) + ...
         ciy(l) .* (normy((2:nx-1)+cix(l),(2:ny-1)+ciy(l),(2:nz-1)+ciz(l))) + ...
         ciz(l) .* (normz((2:nx-1)+cix(l),(2:ny-1)+ciy(l),(2:nz-1)+ciz(l))) ...
        );
    end
    ffx(2:nx-1,2:ny-1,2:nz-1) = sigma .* curvature(2:nx-1,2:ny-1,2:nz-1) .* normx(2:nx-1,2:ny-1,2:nz-1) .* indicator(2:nx-1,2:ny-1,2:nz-1);
    ffy(2:nx-1,2:ny-1,2:nz-1) = sigma .* curvature(2:nx-1,2:ny-1,2:nz-1) .* normy(2:nx-1,2:ny-1,2:nz-1) .* indicator(2:nx-1,2:ny-1,2:nz-1);
    ffz(2:nx-1,2:ny-1,2:nz-1) = sigma .* curvature(2:nx-1,2:ny-1,2:nz-1) .* normz(2:nx-1,2:ny-1,2:nz-1) .* indicator(2:nx-1,2:ny-1,2:nz-1);

    % Momentos
    for i = 2:nx-1
        for j = 2:ny-1
            for k = 2:nz-1
                ux(i,j,k) = sum(f(i,j,k,[2,16,10,8,14])) - sum(f(i,j,k,[3,11,17,15,8]));
                uy(i,j,k) = sum(f(i,j,k,[4,8,15,18,12])) - sum(f(i,j,k,[5,14,9,13,19]));
                uz(i,j,k) = sum(f(i,j,k,[7,16,11,18,13])) - sum(f(i,j,k,[6,10,17,12,19]));
                ux(i,j,k) = ux(i,j,k) ./ rho(i,j,k) + ffx(i,j,k) * 0.5 ./ rho(i,j,k);
                uy(i,j,k) = uy(i,j,k) ./ rho(i,j,k) + ffy(i,j,k) * 0.5 ./ rho(i,j,k);
                uz(i,j,k) = uz(i,j,k) ./ rho(i,j,k) + ffz(i,j,k) * 0.5 ./ rho(i,j,k);
                uu = 0.5 * (ux(i,j,k).^2 + uy(i,j,k).^2 + uz(i,j,k).^2) / cssq;
                rho(i,j,k) = sum(f(i,j,k,:),4);
                for l = 1:fpoints
                    udotc = (ux(i,j,k) * cix(l) + uy(i,j,k) * ciy(l) + uz(i,j,k) * ciz(l)) / cssq;
                    HeF = 0.5 .* (w(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5.*udotc.^2 - uu))) ...
                        .* ((cix(l) - ux(i,j,k)) .* ffx(i,j,k) + ...
                            (ciy(l) - uy(i,j,k)) .* ffy(i,j,k) + ...
                            (ciz(l) - uz(i,j,k)) .* ffz(i,j,k) ...
                           ) ./ (rho(i,j,k) .* cssq);
                    feq = w(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5.*udotc.^2 - uu)) - HeF;
                    fneq(l) = f(i,j,k,l) - feq;
                end
                pxx(i,j,k) = sum(fneq([2,3,8,9,10,11,14,15,16,17]));
                pyy(i,j,k) = sum(fneq([4,5,8,9,12,13,14,15,18,19]));
                pzz(i,j,k) = sum(fneq([6,7,10,11,12,13,16,17,18,19]));
                pxy(i,j,k) = fneq(8) + fneq(9) - fneq(14) - fneq(15);
                pxz(i,j,k) = fneq(10) + fneq(11) - fneq(16) - fneq(17);
                pyz(i,j,k) = fneq(12) + fneq(13) - fneq(18) - fneq(19);
            end
        end
    end

    % Colisão
    for i = 2:nx-1
        for j = 2:ny-1
            for k = 2:nz-1
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
                           (cix(l) .* ciy(l) - cssq) * pxy(i,j,k) + ...
                           (cix(l) .* ciz(l) - cssq) * pxz(i,j,k) + ...
                           (ciy(l) .* ciz(l) - cssq) * pyz(i,j,k);
                    f(i+cix(l),j+ciy(l),k+ciz(l),l) = feq + (1-omega) * (w(l) / (2*cssq^2)) * fneq + HeF;
                end
                for l = 1:gpoints
                    udotc = (ux(i,j,k) * cix(l) + uy(i,j,k) * ciy(l) + uz(i,j,k) * ciz(l)) / cssq;
                    feq = w_g(l) .* phi(i,j,k) .* (1 + udotc);
                    Hi = sharp_c .* phi(i,j,k) .* (1 - phi(i,j,k)) .* (cix(l) .* normx(i,j,k) + ciy(l) .* normy(i,j,k) + ciz(l) .* normz(i,j,k)); 
                    g(i,j,k,l) = feq + w_g(l) .* Hi;
                end
            end
        end
    end

    for l = 1:gpoints
        g(:,:,:,l) = circshift(g(:,:,:,l),[cix(l),ciy(l),ciz(l)]);
    end

    % Condições de contorno
    for i = [1,nx]
        for j = [1,ny]
            for k = [1,nz]
                for l = 1:fpoints
                    if (i+cix(l)>0 && j+ciy(l)>0 && k+ciz(l)>0)
                        f(i+cix(l),j+ciy(l),k+ciz(l),l) = rho(i,j,k) .* w(l); 
                    end
                end
                for l = 1:gpoints
                    if (i+cix(l)>0 && j+ciy(l)>0 && k+ciz(l)>0)
                        g(i+cix(l),j+ciy(l),k+ciz(l),l) = phi(i,j,k) .* w_g(l);
                    end
                end
            end
        end
    end

    phi(:,:,1) = phi(:,:,2);  
    phi(:,:,nz) = phi(:,:,nz-1); 
    phi(1,:,:) = phi(2,:,:); 
    phi(nx,:,:) = phi(nx-1,:,:); 
    phi(:,1,:) = phi(:,2,:); 
    phi(:,ny,:) = phi(:,ny-1,:); 

    if(mod(t,1) == 0)      
        if slicebool == 1
            x = 1:nx; y = 1:ny; z = 1:nz;
            h = slice(x, y, z, phi, [], ny/2, []); 
            shading interp; colorbar; 
            xlabel('X'); ylabel('Y'); zlabel('Z'); 
            title(['t = ', num2str(t)]);
            view(3); drawnow; 
        else
            hVol.Data = phi; 
            drawnow;
        end
    end

    disp(['Passo de tempo: ', num2str(t)]);

end




