%% Lattice Boltzmann Multicomponent Single Distribution
clc; clearvars; close all

%% Thread-Safe LB
nlinks = 9; % Número de direções no modelo D2Q9 (2 dimensões, 9 direções).
tau = 0.505; % Tempo de relaxamento.
cssq = 1.0 / 3.0; % Valor constante relacionado à velocidade do lattice Boltzmann.
visc_LB = cssq * (tau - 0.5); % Viscosidade no modelo lattice Boltzmann.
visc1 = visc_LB; % Define viscosidade1.
tau = 0.505; % Tempo de relaxamento novamente (redundante).
visc_LB = cssq * (tau - 0.5); % Recalcula viscosidade LB.
visc2 = visc_LB; % Define viscosidade2.
taud = 1; % Tempo de relaxamento para a fase de cálculo (não utilizado).
omega = 1 / tau; % Parâmetro de colisão.

one_ov_nu = 1.0 / visc_LB; % Inverso da viscosidade lattice Boltzmann.
sharp_c = 0.15 * 3; % Coeficiente de afiação da interface, termo antidifusão.
sigma = 0.1; % Coeficiente para a força de tensão superficial.

nx = 150; % Número de pontos na direção x.
ny = 150; % Número de pontos na direção y.
nsteps = 10000; % Número de passos de simulação.
stamp = 10; % Intervalo de exibição (não utilizado).

ggx = 0 * 10^-6; % Força externa na direção x (zero neste caso).
ggy = 0; % Força externa na direção y (zero neste caso).

% Inicializa as variáveis de distribuição e campos
f = zeros(nx, ny, 9); % Campos de distribuição de partículas (9 direções).
g = zeros(nx, ny, 5); % Campos de distribuição de fase (5 direções).
rho = 0 .* ones(nx, ny); % Densidade (inicialmente zero).
u = zeros(nx, ny); % Velocidade na direção x.
v = zeros(nx, ny); % Velocidade na direção y.
pxx = ones(nx, ny); % Tensor de pressão na direção xx.
pyy = ones(nx, ny); % Tensor de pressão na direção yy.
pxy = ones(nx, ny); % Tensor de pressão na direção xy.

% Variáveis auxiliares
fi = ones(nx, ny) .* 0; % Solução atual do campo de fase.
normx = zeros(nx, ny); % Gradiente da fase na direção x.
normy = zeros(nx, ny); % Gradiente da fase na direção y.
curvature = zeros(nx, ny); % Curvatura da interface.
indicator = zeros(nx, ny); % Indicador para a curvatura.
bool_ind = zeros(nx, ny); % Indicador booleano (não utilizado).
ffx = zeros(nx, ny); % Força na direção x.
ffy = zeros(nx, ny); % Força na direção y.
grad_rho_x = zeros(nx, ny); % Gradiente da densidade na direção x.
grad_rho_y = zeros(nx, ny); % Gradiente da densidade na direção y.
mod_grad = zeros(nx, ny); % Módulo do gradiente da fase.

% Conjunto de velocidades D2Q9
p = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]; % Pesos para cada direção.
p_g = [2/6, 1/6, 1/6, 1/6, 1/6]; % Pesos para as direções de fase.

% v = [1, 2, 3, 4, 5, 6, 7, 8, 9]; % Direções normais.
opp = [1, 4, 5, 2, 3, 8, 9, 6, 7]; % Direções opostas.

ex = [0, 1, 0, -1, 0, 1, -1, -1, 1]; % Componentes x das direções.
ey = [0, 0, 1, 0, -1, 1, 1, -1, -1]; % Componentes y das direções.
fneq = zeros(9, 1); % Termos fora do equilíbrio.
isfluid = zeros(nx, ny); % Matriz de fluidez.

isfluid(2:nx-1, 2:ny-1) = 1; % Marca o domínio fluido.

% Inicializa os campos de densidade e fase
rho(:,:) = 1; % Define densidade inicial como 1.
% Preenchimento da fase com um perfil de interface
for ii = 2:(nx-1)
    for jj = 2:ny-1
        Ri = sqrt((ii - (nx / 2))^2 / 2.^2 + (jj - (ny / 2))^2); % Distância radial.
        fi(ii, jj) = 0.5 + 0.5 * tanh(2 * (20 - Ri) / 3); % Perfil da fase.
    end
end

% Inicializa as distribuições f e g com base nas densidades
f(:,:,1) = p(1) * rho(:,:);
f(:,:,2) = p(2) * rho(:,:);
f(:,:,3) = p(2) * rho(:,:);
f(:,:,4) = p(2) * rho(:,:);
f(:,:,5) = p(2) * rho(:,:);
f(:,:,6) = p(6) * rho(:,:);
f(:,:,7) = p(6) * rho(:,:);
f(:,:,8) = p(6) * rho(:,:);
f(:,:,9) = p(6) * rho(:,:);
g(:,:,1) = p_g(1) * fi(:,:);
g(:,:,2) = p_g(2) * fi(:,:);
g(:,:,3) = p_g(2) * fi(:,:);
g(:,:,4) = p_g(2) * fi(:,:);
g(:,:,5) = p_g(2) * fi(:,:);

% Loop principal da simulação
for tt = 1:nsteps
    %% Cálculo do campo de fase
    for ii = 1:nx
        for jj = 1:ny
            if isfluid(ii, jj) == 1
                fi(ii, jj) = sum(g(ii, jj, :), 3); % Atualiza o campo de fase.
            end
        end 
    end

    for ii = 1:nx
        for jj = 1:ny
            if isfluid(ii, jj) == 1
                %% Cálculo dos gradientes e arrays
                grad_fix = 0; % Gradiente na direção x.
                grad_fiy = 0; % Gradiente na direção y.
                for kk = 1:9
                    grad_fix = grad_fix + 3 * p(kk) * ex(kk) * fi(ii + ex(kk), jj + ey(kk));
                    grad_fiy = grad_fiy + 3 * p(kk) * ey(kk) * fi(ii + ex(kk), jj + ey(kk));
                end
                mod_grad(ii, jj) = sqrt(grad_fix.^2 + grad_fiy.^2); % Módulo do gradiente.
                normx(ii, jj) = grad_fix / (mod_grad(ii, jj) + 10^-9); % Normalizado na direção x.
                normy(ii, jj) = grad_fiy / (mod_grad(ii, jj) + 10^-9); % Normalizado na direção y.
                indicator(ii, jj) = sqrt(grad_fix.^2 + grad_fiy.^2); % Indicador do gradiente.
            end
        end
    end

    %% Curvatura
    for ii = 1:nx
        for jj = 1:ny
            if isfluid(ii, jj) == 1
                %% Força de tensão superficial
                curvature(ii, jj) = 0; % Inicializa curvatura.
                for kk = 1:9
                    curvature(ii, jj) = curvature(ii, jj) - 3 * p(kk) * (ex(kk) * (normx(ii + ex(kk), jj + ey(kk))) + ...
                        ey(kk) * (normy(ii + ex(kk), jj + ey(kk))));
                end
                ffx(ii, jj) = sigma * curvature(ii, jj) * normx(ii, jj) * indicator(ii, jj);
                ffy(ii, jj) = sigma * curvature(ii, jj) * normy(ii, jj) * indicator(ii, jj);
            end
            end5
    end

    %% Momentos
    for ii = 1:nx
        for jj = 1:ny
            if isfluid(ii, jj) == 1
                u(ii, jj) = (f(ii, jj, 2) - f(ii, jj, 4) + f(ii, jj, 6) - f(ii, jj, 7) - f(ii, jj, 8) + f(ii, jj, 9)) / rho(ii, jj) + ffx(ii, jj) * 0.5 / rho(ii, jj);
                v(ii, jj) = (f(ii, jj, 3) - f(ii, jj, 5) + f(ii, jj, 6) + f(ii, jj, 7) - f(ii, jj, 8) - f(ii, jj, 9)) / rho(ii, jj) + ffy(ii, jj) * 0.5 / rho(ii, jj);
                uu = 0.5 * (u(ii, jj).^2 + v(ii, jj).^2) / cssq; % Velocidade quadrática.
                rho(ii, jj) = sum(f(ii, jj, :), 3); % Atualiza a densidade.
                for kk = 1:9
                    udotc = (u(ii, jj) * ex(kk) + v(ii, jj) * ey(kk)) / cssq;
                    HeF = (p(kk) * (rho(ii, jj) + rho(ii, jj) * (udotc + 0.5 * udotc.^2 - uu))) * ((ex(kk) - u(ii, jj)) * ffx(ii, jj) + (ey(kk) - v(ii, jj)) * ffy(ii, jj)) / (rho(ii, jj) * cssq);
                    feq = p(kk) * (rho(ii, jj) + rho(ii, jj) * (udotc + 0.5 * udotc.^2 - uu)) - 0.5 * (HeF);
                    fneq(kk) = f(ii, jj, kk) - feq; % Termos fora do equilíbrio.
                end
                pxx(ii, jj) = fneq(2) + fneq(4) + fneq(6) + fneq(7) + fneq(8) + fneq(9);
                pyy(ii, jj) = fneq(3) + fneq(5) + fneq(6) + fneq(7) + fneq(8) + fneq(9);
                pxy(ii, jj) = fneq(6) - fneq(7) + fneq(8) - fneq(9);
            end
        end
    end

    %% Colisão
    for ii = 1:nx
        for jj = 1:ny
            if isfluid(ii, jj) == 1
                uu = 0.5 * (u(ii, jj).^2 + v(ii, jj).^2) / cssq;
                for kk = 1:9
                    udotc = (u(ii, jj) * ex(kk) + v(ii, jj) * ey(kk)) / cssq;
                    feq = p(kk) * (rho(ii, jj) + rho(ii, jj) * (udotc + 0.5 * udotc.^2 - uu));
                    HeF = 0.5 * (p(kk) * (rho(ii, jj) + rho(ii, jj) * (udotc + 0.5 * udotc.^2 - uu))) * ((ex(kk) - u(ii, jj)) * ffx(ii, jj) + (ey(kk) - v(ii, jj)) * ffy(ii, jj)) / (rho(ii, jj) * cssq);
                    fneq = (ex(kk) * ex(kk) - cssq) * pxx(ii, jj) + (ey(kk) * ey(kk) - cssq) * pyy(ii, jj) + 2 * ex(kk) * ey(kk) * pxy(ii, jj);
                    f(ii + ex(kk), jj + ey(kk), kk) = feq + (1.0 - omega) * (p(kk) / (2 * cssq^2)) * fneq + HeF; % Atualiza a distribuição.
                end

                for kk = 1:5
                    udotc = (u(ii, jj) * ex(kk) + v(ii, jj) * ey(kk)) / cssq;
                    feq = p_g(kk) * fi(ii, jj) * (1 + udotc);
                    Hi = sharp_c * fi(ii, jj) * (1 - fi(ii, jj)) * (ex(kk) * normx(ii, jj) + ey(kk) * normy(ii, jj));
                    g(ii, jj, kk) = feq + p_g(kk) * Hi; % Atualiza a distribuição de fase.
                end
            end
        end
    end

    % Desloca as distribuições em torno do domínio
    for kk = 1:5
        g(:,:,kk) = circshift(g(:,:,kk), [ex(kk), ey(kk), 0]); % Deslocamento periódico.
    end

    % Condições de contorno
    for ii = 1:nx
        for jj = 1:ny
            if isfluid(ii, jj) == 0
                for kk = 1:9
                    if ii + ex(kk) > 0 && jj + ey(kk) > 0
                        f(ii + ex(kk), jj + ey(kk), kk) = rho(ii, jj) * p(kk); % Condições de contorno de densidade.
                    end
                end
                for kk = 1:5
                    if ii + ex(kk) > 0 && jj + ey(kk) > 0
                        g(ii + ex(kk), jj + ey(kk), kk) = fi(ii, jj) * p_g(kk); % Condições de contorno de fase.
                    end
                end
            end
        end
    end

    % Condições de contorno periódicas
    fi(:,1) = fi(:,2);
    fi(:,ny) = fi(:,ny-1);
    fi(1,:) = fi(2,:);
    fi(nx,:) = fi(nx-1,:);

    % Exibição periódica dos resultados
    if mod(tt, 100) == 0
        imagesc(fi') % Exibe o campo de fase.
        axis xy
        axis equal
        colorbar
        drawnow % Atualiza a figura.
    end

    disp(tt) % Exibe o número do passo atual.
end
