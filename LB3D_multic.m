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

% função de distribuição de equilíbrio local?
f = zeros(nx,ny,nz,19); 
% função de distribuição de equilíbrio?
% g = zeros(nx, ny, nz, ?); 

%% Matrizes e Variáveis

[rho, u, v, w, ...
 fi, normx, normy, ...
 curvature, indicator, ...
 bool_ind, ffx, ffy, ...
 grad_rho_x, grad_rho_y, ...
 mod_grad, isfluid] = deal(zeros(nx,ny,nz));

[pxx, pyy, pzz, ...
 pxy, pxz, pyz] = deal(ones(nx,ny,nz));

p = zeros(1,19);
% p_g = zeros(1, ?);

fneq = zeros(19,1); 
isfluid(2:nx-1,2:ny-1,2:nz-1) = 1;
rho(:,:,:) = 1;

%% Propriedades do Modelo

p(1) = 1/3;
p(2:6) = 1/18;
p(7:19) = 1/36;

% PESOS PARA O CAMPO DE FASE?
% p_g = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% for i = 1:?
%     g(:,:,:,i) = p_g(i) * fi(:,:,:);
% end

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
                     indicator(i,j,k)] = deal(norm(grad_fix,grad_fiy,grad_fiz));
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
                    % Cálculo das velocidades
                        % fazer
                    % Cálculo da densidade
                    rho(i,j,k) = sum(f(i,j,k,:),4);
                    % Cálculo dos momentos
                    for l = 1:19
                        udotc = (u(i,j,k) * cix(l) + v(i,j,k) * ciy(l) + w(i,j,k) * ciz(l));
                        % Cálculo de HeF
                            % fazer
                        % Distribuição de equilíbrio
                            % fazer
                        fneq(l) = f(i,j,k,l) - feq;
                    end
                    % Cálculo de momentos
                        % fazer
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
                        udotc = (u(i,j,k) * cix(l) + v(i,j,k) * ciy(l) + w(i,j,k) * ciz(l))/cssq;
                        % Cálculo da função de distribuição de equilíbrio
                            % fazer
                        % Cálculo de HeF (forças externas ou tensão superficial)
                            % fazer
                        % Cálculo da diferença em relação à função de equilíbrio
                            % fazer
                        % Cálculo da diferença em relação à função de equilíbrio
                            % fazer
                    end
                    % for l = 1:?
                        % fazer
                    % end
                end
            end
        end
    end
    % Realiza a mudança de posição das funções de distribuição para simular movimento
    % for l = 1:?
        % FAZER
    % end

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
                    % for l = 1:?
                    % Atualiza a função de distribuição para o segundo componente
                        % FAZER
                    % end
                end
            end
        end
    end

    % Condições de contorno nos limites da malha
        % FAZER

    if(mod(tt,100) == 0)      
        imagesc(fi')
        axis xy
        axis equal
        colorbar
        drawnow;
    end

    disp(tt)

end


