%% Lattice Boltzmann Multicomponent Single Distribution
clc; clearvars; close all

%% Thread-Safe LB
nlinks = 19; % Número de direções no modelo D3Q19.
tau = 0.505; % Tempo de relaxamento. % original 0.505
cssq = 1/3; % Valor constante relacionado à velocidade do lattice Boltzmann.
visc_LB = cssq * (tau - 0.5); % Viscosidade no modelo lattice Boltzmann.
visc1 = visc_LB; % Define viscosidade 1.
visc2 = visc_LB; % Define viscosidade 2.
omega = 1/tau; % Parâmetro de colisão.

one_ov_nu = 1/visc_LB; % Inverso da viscosidade lattice Boltzmann.
sharp_c = 0.15 * 3; % Coeficiente de afiação da interface, termo antidifusão. % original 0.15 * 3
sigma = 0.1; % Coeficiente para a força de tensão superficial. % original 0.1

npoints = 50;
nx = npoints; % Número de pontos na direção x.
ny = npoints; % Número de pontos na direção y.
nz = npoints; % Número de pontos na direção z.
nsteps = 1000; % Número de passos de simulação.

ggx = 0; % Força externa na direção x (zero neste caso).
ggy = 0; % Força externa na direção y (zero neste caso).
ggz = 0; % Força externa na direção z (zero neste caso).

%% Inicializa as variáveis de distribuição e campos
f = zeros(nx, ny, nz, 19); % Campos de distribuição de partículas (19 direções).
g = zeros(nx, ny, nz, 19); % Campos de distribuição de fase (5 direções).
rho = ones(nx, ny, nz) .* 0; % Densidade (inicialmente zero).

u = zeros(nx, ny, nz); % Velocidade na direção x.
v = zeros(nx, ny, nz); % Velocidade na direção y.
w = zeros(nx, ny, nz); % Velocidade na direção z.

pxx = ones(nx, ny, nz); % Tensor de pressão na direção xx.
pyy = ones(nx, ny, nz); % Tensor de pressão na direção yy.
pzz = ones(nx, ny, nz); % Tensor de pressão na direção zz.
pxy = ones(nx, ny, nz); % Tensor de pressão na direção xy.
pxz = ones(nx, ny, nz); % Tensor de pressão na direção xz
pyz = ones(nx, ny, nz); % Tensor de pressão na direção yz.

%% Variáveis auxiliares
fi = ones(nx, ny, nz) .* 0; % Solução atual do campo de fase.

normx = zeros(nx, ny, nz); % Gradiente da fase na direção x.
normy = zeros(nx, ny, nz); % Gradiente da fase na direção y.
normz = zeros(nx, ny, nz); % Gradiente da fase na direção z.

curvature = zeros(nx, ny, nz); % Curvatura da interface.
indicator = zeros(nx, ny, nz); % Indicador para a curvatura.

ffx = zeros(nx, ny, nz); % Força na direção x.
ffy = zeros(nx, ny, nz); % Força na direção y.
ffz = zeros(nx, ny, nz); % Força na direção z.

grad_rho_x = zeros(nx, ny, nz); % Gradiente da densidade na direção x.
grad_rho_y = zeros(nx, ny, nz); % Gradiente da densidade na direção y.
grad_rho_z = zeros(nx, ny, nz); % Gradiente da densidade na direção z.

mod_grad = zeros(nx, ny, nz); % Módulo do gradiente da fase.

%% Conjunto de velocidades D3Q19
% Pesos
p = [1/3, ...
    1/18, 1/18, 1/18, 1/18, 1/18, 1/18, ...
    1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36];

% Pesos para as direções de fase
p_g = p;

% Direções opostas no D3Q19
opp = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19];

ex = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0]; % Componentes x das direções.
ey = [0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1]; % Componentes y das direções.
ez = [0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1]; % Componentes z das direções.

fneq = zeros(nx, ny, nz, nlinks); % Termos fora do equilíbrio.
isfluid = zeros(nx, ny, nz); % Matriz de fluidez.

isfluid(2:nx-1, 2:ny-1, 2:nz-1) = 1; % Marca o domínio fluido.

%% Inicializa os campos de densidade e fase
rho(:,:,:) = 1; % Define densidade inicial como 1.

%% Preenchimento da fase com um perfil de interface
for ii = 2:nx-1
    for jj = 2:ny-1
        for kk = 2:nz-1
            Ri = sqrt((ii - (nx / 2))^2 / 2.^2 + (jj - (ny / 2))^2 + (kk - (nz / 2))^2); % Distância radial.
            fi(ii, jj, kk) = 0.5 + 0.5 * tanh(2 * (20 - Ri) / 3); % Perfil da fase.
        end
    end
end

%% Inicializa as distribuições f e g com base nas densidades
f(:,:,:,1) = p(1) * rho(:,:,:); % Direção 0
f(:,:,:,2) = p(2) * rho(:,:,:); % Direção [1, 0, 0]
f(:,:,:,3) = p(3) * rho(:,:,:); % Direção [-1, 0, 0]
f(:,:,:,4) = p(4) * rho(:,:,:); % Direção [0, 1, 0]
f(:,:,:,5) = p(5) * rho(:,:,:); % Direção [0, -1, 0]
f(:,:,:,6) = p(6) * rho(:,:,:); % Direção [0, 0, 1]
f(:,:,:,7) = p(7) * rho(:,:,:); % Direção [0, 0, -1]
f(:,:,:,8) = p(8) * rho(:,:,:); % Direção [1, 1, 0]
f(:,:,:,9) = p(9) * rho(:,:,:); % Direção [-1, 1, 0]
f(:,:,:,10) = p(10) * rho(:,:,:); % Direção [-1, -1, 0]
f(:,:,:,11) = p(11) * rho(:,:,:); % Direção [1, -1, 0]
f(:,:,:,12) = p(12) * rho(:,:,:); % Direção [1, 0, 1]
f(:,:,:,13) = p(13) * rho(:,:,:); % Direção [-1, 0, 1]
f(:,:,:,14) = p(14) * rho(:,:,:); % Direção [-1, 0, -1]
f(:,:,:,15) = p(15) * rho(:,:,:); % Direção [1, 0, -1]
f(:,:,:,16) = p(16) * rho(:,:,:); % Direção [0, 1, 1]
f(:,:,:,17) = p(17) * rho(:,:,:); % Direção [0, -1, 1]
f(:,:,:,18) = p(18) * rho(:,:,:); % Direção [0, -1, -1]
f(:,:,:,19) = p(19) * rho(:,:,:); % Direção [0, 1, -1]

numg = 19;
for i = 1:numg
    g(:,:,:,i) = p_g(i) * fi(:,:,:); % Direção 0 
end

videoFilename = 'LBM.mp4'; 
vidObj = VideoWriter(videoFilename, 'MPEG-4');
vidObj.FrameRate = 5; 
open(vidObj);

hVol = volshow(fi, 'RenderingStyle', 'MaximumIntensityProjection');
viewer = hVol.Parent;
hFig = viewer.Parent;

%% Loop principal da simulação
for tt = 1:nsteps
    %% Cálculo do campo de fase
    for ii = 1:nx
        for jj = 1:ny
            for kk = 1:nz
                if isfluid(ii, jj, kk) == 1
                    fi(ii, jj, kk) = sum(g(ii, jj, kk, :), 4); % Atualiza o campo de fase.
                end
            end
        end 
    end

    for ii = 1:nx
        for jj = 1:ny
            for kk = 1:nz
                if isfluid(ii, jj, kk) == 1
                    %% Cálculo dos gradientes e arrays
                    grad_fix = 0; % Gradiente na direção x.
                    grad_fiy = 0; % Gradiente na direção y.
                    grad_fiz = 0; % Gradiente na direção z.
                    
                    for ll = 1:19
                        grad_fix = grad_fix + 3 * p(ll) * ex(ll) * fi(ii + ex(ll), jj + ey(ll), kk + ez(ll));
                        grad_fiy = grad_fiy + 3 * p(ll) * ey(ll) * fi(ii + ex(ll), jj + ey(ll), kk + ez(ll));
                        grad_fiz = grad_fiz + 3 * p(ll) * ez(ll) * fi(ii + ex(ll), jj + ey(ll), kk + ez(ll));
                    end
                    
                    mod_grad(ii, jj, kk) = sqrt(grad_fix.^2 + grad_fiy.^2 + grad_fiz.^2); % Módulo do gradiente.
                    normx(ii, jj, kk) = grad_fix / (mod_grad(ii, jj, kk) + 10^-9); % Normalizado na direção x.
                    normy(ii, jj, kk) = grad_fiy / (mod_grad(ii, jj, kk) + 10^-9); % Normalizado na direção y.
                    normz(ii, jj, kk) = grad_fiz / (mod_grad(ii, jj, kk) + 10^-9); % Normalizado na direção z.
                    indicator(ii, jj, kk) = sqrt(grad_fix.^2 + grad_fiy.^2 + grad_fiz.^2); % Indicador do gradiente.
                end
            end
        end
    end

    %% Curvatura
    for ii = 1:nx
        for jj = 1:ny
            for kk = 1:nz
                if isfluid(ii, jj, kk) == 1
                    %% Inicializa curvatura e forças de tensão superficial
                    curvature(ii, jj, kk) = 0; % Curvatura
                    ffx(ii, jj, kk) = 0; % Força de tensão superficial na direção x
                    ffy(ii, jj, kk) = 0; % Força de tensão superficial na direção y
                    ffz(ii, jj, kk) = 0; % Força de tensão superficial na direção z
    
                    %% Cálculo da curvatura
                    for ll = 1:19
                        curvature(ii, jj, kk) = curvature(ii, jj, kk) - 3 * p(ll) * ( ...
                            ex(ll) * (normx(ii + ex(ll), jj + ey(ll), kk + ez(ll))) + ...
                            ey(ll) * (normy(ii + ex(ll), jj + ey(ll), kk + ez(ll))) + ...
                            ez(ll) * (normz(ii + ex(ll), jj + ey(ll), kk + ez(ll))) );
                    end
    
                    %% Cálculo das forças de tensão superficial
                    ffx(ii, jj, kk) = sigma * curvature(ii, jj, kk) * normx(ii, jj, kk) * indicator(ii, jj, kk);
                    ffy(ii, jj, kk) = sigma * curvature(ii, jj, kk) * normy(ii, jj, kk) * indicator(ii, jj, kk);
                    ffz(ii, jj, kk) = sigma * curvature(ii, jj, kk) * normz(ii, jj, kk) * indicator(ii, jj, kk);
                end
            end
        end
    end

    %% Momentos
    for ii = 1:nx
        for jj = 1:ny
            for kk = 1:nz
                if isfluid(ii, jj, kk) == 1
                    % Calcular velocidades u, v, w
                    u(ii, jj, kk) = (f(ii, jj, kk, 2) - f(ii, jj, kk, 5) + f(ii, jj, kk, 8) - f(ii, jj, kk, 9) + f(ii, jj, kk, 14) - f(ii, jj, kk, 17)) / rho(ii, jj, kk) + ffx(ii, jj, kk) * 0.5 / rho(ii, jj, kk);
                    v(ii, jj, kk) = (f(ii, jj, kk, 3) - f(ii, jj, kk, 6) + f(ii, jj, kk, 8) + f(ii, jj, kk, 9) - f(ii, jj, kk, 14) - f(ii, jj, kk, 17)) / rho(ii, jj, kk) + ffy(ii, jj, kk) * 0.5 / rho(ii, jj, kk);
                    w(ii, jj, kk) = (f(ii, jj, kk, 4) - f(ii, jj, kk, 7) + f(ii, jj, kk, 10) - f(ii, jj, kk, 11) - f(ii, jj, kk, 15) + f(ii, jj, kk, 18)) / rho(ii, jj, kk) + ffz(ii, jj, kk) * 0.5 / rho(ii, jj, kk);
                    
                    % Velocidade quadrática
                    uu = 0.5 * (u(ii, jj, kk).^2 + v(ii, jj, kk).^2 + w(ii, jj, kk).^2) / cssq;
                    
                    % Atualiza a densidade
                    rho(ii, jj, kk) = sum(f(ii, jj, kk, :), 4);
                    
                    % Calcular momentos fora do equilíbrio

                    for ll = 1:19
                        udotc = (u(ii, jj, kk) * ex(ll) + v(ii, jj, kk) * ey(ll) + w(ii, jj, kk) * ez(ll)) / cssq;
                        HeF = (p(ll) * (rho(ii, jj, kk) + rho(ii, jj, kk) * (udotc + 0.5 * udotc.^2 - uu))) * ...
                            ((ex(ll) - u(ii, jj, kk)) * ffx(ii, jj, kk) + (ey(ll) - v(ii, jj, kk)) * ffy(ii, jj, kk) + (ez(ll) - w(ii, jj, kk)) * ffz(ii, jj, kk)) / (rho(ii, jj, kk) * cssq);
                        feq = p(ll) * (rho(ii, jj, kk) + rho(ii, jj, kk) * (udotc + 0.5 * udotc.^2 - uu)) - 0.5 * (HeF);
                        fneq(ll) = f(ii, jj, kk, ll) - feq; % Termos fora do equilíbrio
                    end
                    
                    % Calcular tensores de pressão
                    pxx(ii, jj, kk) = fneq(2) + fneq(5) + fneq(8) + fneq(9) + fneq(14) + fneq(17);
                    pyy(ii, jj, kk) = fneq(3) + fneq(6) + fneq(8) + fneq(9) + fneq(14) + fneq(17);
                    pzz(ii, jj, kk) = fneq(4) + fneq(7) + fneq(10) + fneq(11) + fneq(15) + fneq(18);
                    pxy(ii, jj, kk) = fneq(8) - fneq(9) + fneq(14) - fneq(17);
                    pxz(ii, jj, kk) = fneq(10) - fneq(11) + fneq(15) - fneq(18);
                    pyz(ii, jj, kk) = fneq(12) - fneq(13) + fneq(16) - fneq(19);
                end
            end
        end
    end

    %% Colisão
    for ii = 1:nx
        for jj = 1:ny
            for kk = 1:nz
                if isfluid(ii, jj, kk) == 1
                    % Calcular velocidades
                    uu = 0.5 * (u(ii, jj, kk).^2 + v(ii, jj, kk).^2 + w(ii, jj, kk).^2) / cssq;
                    
                    % Atualizar as distribuições de função
                    for ll = 1:19
                        udotc = (u(ii, jj, kk) * ex(ll) + v(ii, jj, kk) * ey(ll) + w(ii, jj, kk) * ez(ll)) / cssq;
                        feq = p(ll) * (rho(ii, jj, kk) + rho(ii, jj, kk) * (udotc + 0.5 * udotc.^2 - uu));
                        HeF = 0.5 * (p(ll) * (rho(ii, jj, kk) + rho(ii, jj, kk) * (udotc + 0.5 * udotc.^2 - uu))) * ...
                              ((ex(ll) - u(ii, jj, kk)) * ffx(ii, jj, kk) + (ey(ll) - v(ii, jj, kk)) * ffy(ii, jj, kk) + (ez(ll) - w(ii, jj, kk)) * ffz(ii, jj, kk)) / (rho(ii, jj, kk) * cssq);
                        fneq = (ex(ll)^2 - cssq) * pxx(ii, jj, kk) + (ey(ll)^2 - cssq) * pyy(ii, jj, kk) + (ez(ll)^2 - cssq) * pzz(ii, jj, kk) + ...
                               2 * ex(ll) * ey(ll) * pxy(ii, jj, kk) + 2 * ex(ll) * ez(ll) * pxz(ii, jj, kk) + 2 * ey(ll) * ez(ll) * pyz(ii, jj, kk);
                        f(ii + ex(ll), jj + ey(ll), kk + ez(ll), ll) = feq + (1.0 - omega) * (p(ll) / (2 * cssq^2)) * fneq + HeF; % Atualiza a distribuição
                    end
    
                    % Atualizar distribuições de fase
                    for ll = 1:numg
                        udotc = (u(ii, jj, kk) * ex(ll) + v(ii, jj, kk) * ey(ll) + w(ii, jj, kk) * ez(ll)) / cssq;
                        feq = p_g(ll) * fi(ii, jj, kk) * (1 + udotc);
                        Hi = sharp_c * fi(ii, jj, kk) * (1 - fi(ii, jj, kk)) * (ex(ll) * normx(ii, jj, kk) + ey(ll) * normy(ii, jj, kk) + ez(ll) * normz(ii, jj, kk));
                        g(ii, jj, kk, ll) = feq + p_g(ll) * Hi; % Atualiza a distribuição de fase
                    end
                end
            end
        end
    end

    % Desloca as distribuições em torno do domínio
    for ll = 1:19
        g(:,:,:,ll) = circshift(g(:,:,:,ll), [ex(ll), ey(ll), ez(ll)]); % Deslocamento periódico em 3D.
    end


    % Condições de contorno
    for ii = 1:nx
        for jj = 1:ny
            for kk = 1:nz
                if isfluid(ii, jj, kk) == 0
                    for ll = 1:19
                        if ii + ex(ll) > 0 && jj + ey(ll) > 0 && kk + ez(ll) > 0
                            f(ii + ex(ll), jj + ey(ll), kk + ez(ll), ll) = rho(ii, jj, kk) * p(ll); % Condições de contorno de densidade.
                        end
                    end
                    for ll = 1:numg
                        if ii + ex(ll) > 0 && jj + ey(ll) > 0 && kk + ez(ll) > 0
                            g(ii + ex(ll), jj + ey(ll), kk + ez(ll), ll) = fi(ii, jj, kk) * p_g(ll); % Condições de contorno de fase.
                        end
                    end
                end
            end
        end
    end

    % Condições de contorno periódicas (ou bordas)
    fi(:, :, 1) = fi(:, :, 2); % Fronteira em z=1
    fi(:, :, nz) = fi(:, :, nz-1); % Fronteira em z=nz
    fi(1, :, :) = fi(2, :, :); % Fronteira em x=1
    fi(nx, :, :) = fi(nx-1, :, :); % Fronteira em x=nx
    fi(:, 1, :) = fi(:, 2, :); % Fronteira em y=1
    fi(:, ny, :) = fi(:, ny-1, :); % Fronteira em y=ny

    % Exibição periódica dos resultados
    if mod(tt, 100) == 0
        hVol.Data = fi;
        drawnow;
        frame = getframe(hFig); 
        writeVideo(vidObj, frame);
    end

    disp(tt)
    
end

close(vidObj);


