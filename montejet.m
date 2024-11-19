%% D2Q9 Lattice Boltzmann - Jato Turbulento (Uma Fase)
clc; clearvars; close all

% Parâmetros do Lattice Boltzmann
nlinks = 9;
tau = 0.505;  % Tempo de relaxação próximo de 0.5 para minimizar difusão
cssq = 1.0 / 3.0; % Velocidade do som ao quadrado
nu = cssq * (tau - 0.5); % Viscosidade cinemática
omega = 1 / tau;

% Número de Reynolds
L = 20;        % Largura do jato (em unidades de lattice)
U = 0.2;       % Velocidade característica do jato
Re = (U * L) / nu; % Número de Reynolds
disp(['Número de Reynolds: ', num2str(Re)]);

% Dimensões do domínio
nx = 300; % Comprimento do domínio
ny = 150; % Altura do domínio

% Passos de tempo
nsteps = 30000;
stamp = 100;

% Inicialização de variáveis
f = zeros(nx, ny, nlinks);
rho = ones(nx, ny);  % Densidade inicial
u = zeros(nx, ny);   % Velocidade em x
v = zeros(nx, ny);   % Velocidade em y

% Padrões do LB
p = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
ex = [0, 1, 0, -1, 0, 1, -1, -1, 1];
ey = [0, 0, 1, 0, -1, 1, 1, -1, -1];

% Condição inicial do jato - centro do domínio
jet_center = round(ny / 2); % Centro do domínio
jet_width = round(L / 2);   % Metade da largura do jato

% Velocidade inicial do jato
u(:, jet_center - jet_width:jet_center + jet_width) = U;  % Velocidade do jato
v(:, jet_center - jet_width:jet_center + jet_width) = 0; % Perturbações mínimas

% Inicializar distribuições
for kk = 1:nlinks
    udotc = ex(kk) * u + ey(kk) * v;
    feq = p(kk) * (rho + 3 * rho .* udotc + 4.5 * rho .* udotc.^2 - 1.5 * rho .* (u.^2 + v.^2));
    f(:, :, kk) = feq;
end

% Loop no tempo
for tt = 1:nsteps
    % Propagação
    for kk = 1:nlinks
        f(:, :, kk) = circshift(f(:, :, kk), [ex(kk), ey(kk)]);
    end

    % Condições de contorno - entrada do jato
    for jj = jet_center - jet_width:jet_center + jet_width
        rho(1, jj) = 1; % Densidade constante na entrada
       
        for kk = 1:nlinks
            udotc = ex(kk) * u(1, jj) + ey(kk) * v(1, jj);
            feq = p(kk) * (rho(1, jj) + 3 * rho(1, jj) * udotc ...
                + 4.5 * rho(1, jj) * udotc^2 - 1.5 * rho(1, jj) * (u(1, jj)^2 + v(1, jj)^2));
            f(1, jj, kk) = feq; % Atualiza distribuição de entrada
        end
    end

    % Saída aberta
    f(nx, :, :) = f(nx-1, :, :);

    % Colisão
    rho = sum(f, 3);
    u = (sum(f(:, :, [2, 6, 9]), 3) - sum(f(:, :, [4, 7, 8]), 3)) ./ rho;
    v = (sum(f(:, :, [3, 6, 7]), 3) - sum(f(:, :, [5, 8, 9]), 3)) ./ rho;

    for kk = 1:nlinks
        udotc = ex(kk) * u + ey(kk) * v;
        feq = p(kk) * (rho + 3 * rho .* udotc + 4.5 * rho .* udotc.^2 - 1.5 * rho .* (u.^2 + v.^2));
        f(:, :, kk) = f(:, :, kk) - omega * (f(:, :, kk) - feq);
    end

    % Visualização
    if mod(tt, stamp) == 0
        imagesc(sqrt(u.^2 + v.^2)');    
        axis equal off;
        colorbar;
        title(['Tempo: ', num2str(tt)]);
        drawnow;
    end
end
