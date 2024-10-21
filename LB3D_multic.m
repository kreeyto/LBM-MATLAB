%% D3Q19 
clc; clearvars; close all

%% Parâmetros Gerais

video = 0;

slicebool = 1;
nlinks = 19;
tau = 0.505;
cssq = 1/3;
omega = 1/tau;
sharp_c = 0.15*3;
sigma = 0.1;

[nx, ny, nz] = deal(15);
nsteps = 10000; 

f = zeros(nx,ny,nz,19); 
g = zeros(nx,ny,nz,15); 

%% Matrizes e Variáveis

[rho, u, v, w, ...
 fi, normx, normy, normz, ...
 curvature, indicator, ...
 bool_ind, ffx, ffy, ffz, ...
 grad_rho_x, grad_rho_y, grad_rho_z, ...
 mod_grad, isfluid] = deal(zeros(nx,ny,nz));

[pxx, pyy, pzz, ...
 pxy, pxz, pyz] = deal(ones(nx,ny,nz));

p = zeros(1,19);
p_g = zeros(1, 15);

fneq = zeros(19,1); 
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

cix = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0];
ciy = [0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1];
ciz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1];

%% Cálculo da Função de Distribuição em Função da Distância Radial

for i = 2:nx-1
    for j = 2:ny-1
        for k = 2:nz-1
            Ri = sqrt( ...
                    (i-(nx/2))^2/2.^2 + ...
                    (j-(ny/2))^2 + ...
                    (k-(nz/2))^2 ...
                );
            fi(i,j,k) = 1/2 + 1/2 * tanh(2*(20-Ri)/3);
        end
    end
end

%% Inicialização de Funções de Distribuição

for i = 1:19
    f(:,:,:,i) = p(i) * rho(:,:,:);
end

for i = 1:15
    g(:,:,:,i) = p_g(i) * fi(:,:,:);
end

%% Inicialização da Visualização

if video == true
    videoFilename = 'LBM.mp4'; 
    vidObj = VideoWriter(videoFilename, 'MPEG-4');
    vidObj.FrameRate = 30; 
    open(vidObj);
end

if slicebool == false
    hVol = volshow(fi, 'RenderingStyle', 'MaximumIntensityProjection');
    viewer = hVol.Parent;
    hFig = viewer.Parent;
end

%% Loop de Simulação

for t = 1:nsteps

    % Campo de fase
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                fi(i,j,k) = sum(g(i,j,k,:),4);
            end
        end
    end

    % Normal e arrays
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                if(isfluid(i,j,k) == 1)
                    for l = 1:19
                        grad_fix = 3 * p(l) .* cix(l) .* ((fi(i+cix(l),j+ciy(l),k+ciz(l))));
                        grad_fiy = 3 * p(l) .* ciy(l) .* ((fi(i+cix(l),j+ciy(l),k+ciz(l))));
                        grad_fiz = 3 * p(l) .* ciz(l) .* ((fi(i+cix(l),j+ciy(l),k+ciz(l))));
                    end
                    [mod_grad(i,j,k), ...
                     indicator(i,j,k)] = deal(sqrt(grad_fix.^2 + grad_fiy.^2 + grad_fiz.^2));
                    normx(i,j,k) = grad_fix ./ (mod_grad(i,j,k) + 10^-9);
                    normy(i,j,k) = grad_fiy ./ (mod_grad(i,j,k) + 10^-9);
                    normz(i,j,k) = grad_fiz ./ (mod_grad(i,j,k) + 10^-9);
                end
            end
        end
    end

    % Curvatura
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                if(isfluid(i,j,k) == 1)
                    curvature(i,j,k) = 0;
                    for l = 1:19
                        curvature(i,j,k) = curvature(i,j,k) - 3 .* p(l) .* (cix(l) .* (normx(i+cix(l),j+ciy(l),k+ciz(l))) + ...
                                                                            ciy(l) .* (normy(i+cix(l),j+ciy(l),k+ciz(l))) + ...
                                                                            ciz(l) .* (normz(i+cix(l),j+ciy(l),k+ciz(l))));
                    end
                    ffx(i,j,k) = sigma .* curvature(i,j,k) .* normx(i,j,k) .* indicator(i,j,k);
                    ffy(i,j,k) = sigma .* curvature(i,j,k) .* normy(i,j,k) .* indicator(i,j,k);
                    ffz(i,j,k) = sigma .* curvature(i,j,k) .* normz(i,j,k) .* indicator(i,j,k);
                end
            end
        end
    end

    % Momentos
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                if(isfluid(i,j,k) == 1)
                    
                    % LBM
                    u(i, j, k) = ( (f(i,j,k,7) + f(i,j,k,16) + f(i,j,k,11) + f(i,j,k,18) + f(i,j,k,13)) - (f(i,j,k,6) + f(i,j,k,10) + f(i,j,k,17) + f(i,j,k,12) + f(i,j,k,19)) ) ./ rho(i, j, k) + ffx(i, j, k) * 0.5 ./ rho(i, j, k);
                    v(i, j, k) = ( (f(i,j,k,4) + f(i,j,k,8) + f(i,j,k,15) + f(i,j,k,18) + f(i,j,k,12)) - (f(i,j,k,5) + f(i,j,k,14) + f(i,j,k,9) + f(i,j,k,13) + f(i,j,k,19)) ) ./ rho(i, j, k) + ffy(i, j, k) * 0.5 ./ rho(i, j, k);
                    w(i, j, k) = ( (f(i,j,k,2) + f(i,j,k,16) + f(i,j,k,10) + f(i,j,k,8) + f(i,j,k,14)) - (f(i,j,k,3) + f(i,j,k,11) + f(i,j,k,17) + f(i,j,k,15) + f(i,j,k,9)) ) ./ rho(i, j, k) + ffz(i, j, k) * 0.5 ./ rho(i, j, k);

                    % SUKOP
                    % u(i, j, k) = ( (f(i,j,k,2) + f(i,j,k,16) + f(i,j,k,10) + f(i,j,k,8) + f(i,j,k,14)) - (f(i,j,k,3) + f(i,j,k,11) + f(i,j,k,17) + f(i,j,k,15) + f(i,j,k,9)) ) ./ rho(i, j, k) + ffx(i, j, k) * 0.5 ./ rho(i, j, k);
                    % v(i, j, k) = ( (f(i,j,k,7) + f(i,j,k,16) + f(i,j,k,11) + f(i,j,k,18) + f(i,j,k,13)) - (f(i,j,k,6) + f(i,j,k,10) + f(i,j,k,17) + f(i,j,k,12) + f(i,j,k,19)) ) ./ rho(i, j, k) + ffy(i, j, k) * 0.5 ./ rho(i, j, k);
                    % w(i, j, k) = ( (f(i,j,k,4) + f(i,j,k,8) + f(i,j,k,15) + f(i,j,k,18) + f(i,j,k,12)) - (f(i,j,k,5) + f(i,j,k,14) + f(i,j,k,9) + f(i,j,k,13) + f(i,j,k,19)) ) ./ rho(i, j, k) + ffz(i, j, k) * 0.5 ./ rho(i, j, k);
                    
                    uu = 0.5 * (u(i,j,k).^ 2 + v(i,j,k).^ 2 + w(i,j,k).^2) / cssq;
                    rho(i,j,k) = sum(f(i,j,k,:),4);
                    for l = 1:19
                        udotc = (u(i,j,k) * cix(l) + v(i,j,k) * ciy(l) + w(i,j,k) * ciz(l));
                        HeF = (p(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5 .* udotc.^2 - uu))) ...
                            .* ((cix(l) - u(i,j,k)) .* ffx(i,j,k) + ...
                                (ciy(l) - v(i,j,k)) .* ffy(i,j,k) + ...
                                (ciz(l) - w(i,j,k)) .* ffz(i,j,k) ...
                               ) ./ (rho(i,j,k) .* cssq);
                        feq = p(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5 .* udotc.^2 - uu)) - 0.5 .* (HeF);
                        fneq(l) = f(i,j,k,l) - feq;
                    end

                    % CHECAR    
                    %{
                    pxx(i,j) = fneq(2) + fneq(4) + fneq(6) + fneq(7) + fneq(8) + fneq(9);
                    pyy(i,j) = fneq(3) + fneq(5) + fneq(6) + fneq(7) + fneq(8) + fneq(9);
                    pxy(i,j) = fneq(6) - fneq(7) + fneq(8) - fneq(9); 
                    %}

                end
            end
        end
    end

    % Colisão
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                if(isfluid(i,j,k) == 1)
                    uu = 0.5 * (u(i,j,k).^2 + v(i,j,k).^2 + w(i,j,k).^2) / cssq;
                    for l = 1:19
                        udotc = (u(i,j,k) * cix(l) + v(i,j,k) * ciy(l) + w(i,j,k) * ciz(l)) / cssq;
                        feq = p(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5 .* udotc.^2 - uu));
                        HeF = 0.5 .* (p(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5 .* udotc .^2 - uu))) .* ...
                            ((cix(l) - u(i,j,k)) .* ffx(i,j,k) + ...
                             (ciy(l) - v(i,j,k)) .* ffy(i,j,k) + ...
                             (ciz(l) - w(i,j,k)) .* ffz(i,j,k) ...
                            ) ./ (rho(i,j,k) .* cssq);
                        fneq = (cix(l) .* cix(l) - cssq) * pxx(i,j,k) + ...
                               (ciy(l) .* ciy(l) - cssq) * pyy(i,j,k) + ...
                               (ciz(l) .* ciz(l) - cssq) * pzz(i,j,k) + ...
                                2 * cix(l) .* ciy(l) .* pxy(i,j,k) + ...
                                2 * cix(l) .* ciz(l) .* pxz(i,j,k) + ...
                                2 * ciy(l) .* ciz(l) .* pyz(i,j,k);
                        f(i + cix(l), j + ciy(l), k + ciz(l), l) = feq + (1-omega) * (p(l) / (2*cssq^2)) * fneq + HeF;
                    end
                    for l = 1:15
                        udotc = (u(i,j,k) * cix(l) + v(i,j,k) * ciy(l) + w(i,j,k) * ciz(l)) / cssq;
                        feq = p_g(l) .* fi(i,j,k) .* (1 + udotc);
                        Hi = sharp_c .* fi(i,j,k) .* (1 - fi(i,j,k)) .* (cix(l) .* normx(i,j,k) + ciy(l) .* normy(i,j,k) + ciz(l) .* normz(i,j,k)); 
                        g(i,j,k,l) = feq + p_g(l) .* Hi;
                    end
                end
            end
        end
    end

    for l = 1:15
        g(:,:,:,l) = circshift(g(:,:,:,l),[cix(l),ciy(l),ciz(l),0]);
    end

    % Condições de contorno
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                if(isfluid(i,j,k) == 1)
                    for l = 1:19
                        if(i + cix(l) > 0 && j + ciy(l) > 0 && k + ciz(l) > 0)
                            f(i + cix(l), j + ciy(l), k + ciz(l), l) = rho(i,j,k) .* p(l); 
                        end
                    end
                    for l = 1:15
                        if(i + cix(l) > 0 && j + ciy(l) > 0 && k + ciz(l) > 0)
                            g(i + cix(l), j + ciy(l), k + ciz(l), l) = fi(i,j,k) .* p_g(l);
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

    if(mod(t,100) == 0)      
        if slicebool == 1
            hFig = figure(1); clf;
            x = 1:nx; y = 1:ny; z = 1:nz;
            h = slice(x, y, z, fi, [], ny/2, []); 
            set(h, 'EdgeColor', 'none'); 
            shading interp; colorbar; 
            xlabel('X'); ylabel('Y'); zlabel('Z'); 
            title(['Visualização 3D do campo de fase no tempo t = ', num2str(t)]);
            view(3); drawnow; 
        elseif slicebool == 2
            hFig = figure(1); clf;
            x = 1:nx; y = 1:ny; z = 1:nz;
            surfpatch = patch(isosurface(x, y, z, fi));
            set(surfpatch, 'FaceColor', 'red', 'EdgeColor', 'none'); 
            xlabel('X'); ylabel('Y'); zlabel('Z');
            axis equal;
            camlight; lighting phong; 
            title(['Isosurface do campo de fase no tempo t = ', num2str(t)]);
            view(3); drawnow;
        else
            hVol.Data = fi;
            drawnow;
        end
        if video == true
            frame = getframe(hFig); 
            writeVideo(vidObj, frame);
        end
    end

    disp(t)

end


