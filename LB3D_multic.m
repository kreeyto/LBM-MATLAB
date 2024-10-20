%% D3Q19 
clc; clearvars; close all

%% Parâmetros Gerais

nlinks = 19;
tau = 0.505;
cssq = 1/3;
omega = 1/tau;
sharp_c = 0.15*3;
sigma = 0.1;

[nx, ny, nz] = deal(150);
nsteps = 10000; 

f = zeros(nx,ny,nz,19); 
g = zeros(nx, ny, nz, 15); 

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
                    (i-(nx/2))^2/4 + ...
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

%% Loop de Simulação

for t = 1:nsteps

    % Cálculo do campo de fase
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                fi(i,j,k) = sum(g(i,j,k,:),4);
            end
        end
    end

    % Cálculo das normais e das matrizes
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                if(isfluid(i,j,k) == 1)
                    for l = 1:19
                        % Gradientes da função de campo de fase
                        grad_fix = 3 * p(l) .* cix(l) .* ((fi(i+cix(l),j+ciy(l),k+ciz(l))));
                        grad_fiy = 3 * p(l) .* ciy(l) .* ((fi(i+cix(l),j+ciy(l),k+ciz(l))));
                        grad_fiz = 3 * p(l) .* ciz(l) .* ((fi(i+cix(l),j+ciy(l),k+ciz(l))));
                    end
                    % Magnitude do gradiente e indicador
                    [mod_grad(i,j,k), ...
                     indicator(i,j,k)] = deal(sqrt(grad_fix.^2 + grad_fiy.^2 + grad_fiz.^2));
                    % Normalização dos gradientes
                    normx(i,j,k) = grad_fix ./ (mod_grad(i,j,k) + 10^-9);
                    normy(i,j,k) = grad_fiy ./ (mod_grad(i,j,k) + 10^-9);
                    normz(i,j,k) = grad_fiz ./ (mod_grad(i,j,k) + 10^-9);
                end
            end
        end
    end

    % Cálculo da curvatura e das forças de tensão superficial
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                if(isfluid(i,j,k) == 1)
                    % Curvatura
                    curvature(i,j,k) = 0;
                    for l = 1:19
                        curvature(i,j,k) = curvature(i,j,k) - 3 .* p(l) .* (cix(l) .* (normx(i+cix(l),j+ciy(l),k+ciz(l))) + ...
                                                                            ciy(l) .* (normy(i+cix(l),j+ciy(l),k+ciz(l))) + ...
                                                                            ciz(l) .* (normz(i+cix(l),j+ciy(l),k+ciz(l))));
                    end
                    % Forças de tensão superficial
                    ffx(i,j,k) = sigma .* curvature(i,j,k) .* normx(i,j,k) .* indicator(i,j,k);
                    ffy(i,j,k) = sigma .* curvature(i,j,k) .* normy(i,j,k) .* indicator(i,j,k);
                    ffz(i,j,k) = sigma .* curvature(i,j,k) .* normz(i,j,k) .* indicator(i,j,k);
                end
            end
        end
    end

    % Cálculo dos momentos
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                if(isfluid(i,j,k) == 1)
                    % AJUSTAR PARA 3D
                    % Cálculo das velocidades
                        % u(ii,jj)=(f(ii,jj,2)-f(ii,jj,4)+f(ii,jj,6)-f(ii,jj,7)-f(ii,jj,8)+f(ii,jj,9))./rho(ii,jj) + ffx(ii,jj)*0.5./rho(ii,jj) ;
                        % v(ii,jj)=(f(ii,jj,3)-f(ii,jj,5)+f(ii,jj,6)+f(ii,jj,7)-f(ii,jj,8)-f(ii,jj,9))./rho(ii,jj) + ffy(ii,jj)*0.5./rho(ii,jj) ;   
                    uu = 0.5 * (u(i,j,k) .^ 2 + v(i,j,k) .^ 2) / cssq;
                    % Cálculo da densidade 
                    rho(i,j,k) = sum(f(i,j,k,:),4);
                    % Cálculo dos momentos
                    for l = 1:19
                        udotc = (u(i,j,k) * cix(l) + v(i,j,k) * ciy(l) + w(i,j,k) * ciz(l));
                        % Cálculo de HeF (inseguro sobre primeira linha)
                        HeF = (p(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5 .* udotc.^2 - uu))) ...
                            .* ((cix(l) - u(i,j,k)) .* ffx(i,j,k) + ...
                                (ciy(l) - v(i,j,k)) .* ffy(i,j,k) + ...
                                (ciz(l) - w(i,j,k)) .* ffz(i,j,k) ...
                               ) ./ (rho(i,j,k) .* cssq);
                        % Distribuição de equilíbrio
                        feq = p(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5 .* udotc.^2 - uu)) - 0.5 .* (HeF);
                        fneq(l) = f(i,j,k,l) - feq;
                    end
                    % AJUSTAR PARA 3D
                    % Cálculo de momentos
                        % pxx(ii,jj)= fneq(2) + fneq(4) + fneq(6) + fneq(7) + fneq(8) + fneq(9);
                        % pyy(ii,jj)= fneq(3) + fneq(5) + fneq(6) + fneq(7) + fneq(8) + fneq(9);
                        % pxy(ii,jj)= fneq(6) - fneq(7) + fneq(8) - fneq(9); 
                end
            end
        end
    end

    % Cálculo da colisão
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                if(isfluid(i,j,k) == 1)
                    % Cálculo da energia cinética média
                    uu = 0.5 * (u(i,j,k).^2 + v(i,j,k).^2 + w(i,j,k).^2)/cssq;
                    for l = 1:19
                        % Produto escalar da velocidade do fluido e da direção da partícula
                        udotc = (u(i,j,k) * cix(l) + v(i,j,k) * ciy(l) + w(i,j,k) * ciz(l)) / cssq;
                        % Cálculo da função de distribuição de equilíbrio
                        feq = p(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5 .* udotc.^2 - uu));
                        % Cálculo de HeF (forças externas ou tensão superficial)
                        HeF = 0.5 .* (p(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5 .* udotc .^2 - uu))) .* ...
                            ((cix(l) - u(i,j,k)) .* ffx(i,j,k) + ...
                             (ciy(l) - v(i,j,k)) .* ffy(i,j,k) + ...
                             (ciz(l) - w(i,j,k)) .* ffz(i,j,k) ...
                            ) ./ (rho(i,j,k) .* cssq);
                        % AJUSTAR PARA 3D
                        % Cálculo da diferença em relação à função de equilíbrio
                            % fneq=(ex(kk).*ex(kk)-cssq)*pxx(ii,jj)+(ey(kk).*ey(kk)-cssq)*pyy(ii,jj) ...
                            %    + 2*ex(kk).*ey(kk).*pxy(ii,jj);
                        % Cálculo da diferença em relação à função de equilíbrio
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

    % Realiza a mudança de posição das funções de distribuição para simular movimento
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
                            % Atualiza a função de distribuição
                            f(i + cix(l), j + ciy(l), k + ciz(l), l) = rho(i,j,k) .* p(l); 
                        end
                    end
                    for l = 1:15
                        if(i + cix(l) > 0 && j + ciy(l) > 0 && k + ciz(l) > 0)
                            % Atualiza a função de distribuição para o segundo componente
                            g(i + cix(l), j + ciy(l), k + ciz(l), l) = fi(i,j,k) .* p_g(l);
                        end
                    end
                end
            end
        end
    end

    % AJUSTAR PARA 3D
    % Condições de contorno nos limites da malha
        % fi(:,1)=fi(:,2);
        % fi(:,ny)=fi(:,ny-1);
        % fi(1,:)=fi(2,:);
        % fi(nx,:)=fi(nx-1,:);

    if(mod(tt,100) == 0)      
        imagesc(fi')
        axis xy
        axis equal
        colorbar
        drawnow;
    end

    disp(tt)

end


