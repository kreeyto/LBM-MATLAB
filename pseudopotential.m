clc; clearvars; close all

% parametros
nx = 64; ny = 64;
npop = 9;

tau_water = 1.0; tau_air = 1.0;
omega_water = 1 / tau_water;
omega_air   = 1 / tau_air;

rho_water0 = 1.0;
rho_air0 = 1.0;

cs2 = 1/3;
G = -1;
g = 0.01;

radius = 10;
nsteps = 1000;
output_interval = 1;

% D2Q9
w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
cx = [0, 1, 0, -1, 0, 1, -1, -1, 1];
cy = [0, 0, 1, 0, -1, 1, 1, -1, -1];

% pseudopotencial
psi = @(rho,rho0) rho0 .* (1 - exp(-rho./rho));

% inicializaçao
[xg, yg] = meshgrid(1:nx,1:ny);
mask = (xg-nx/2).^2 + (yg-ny/2).^2 <= radius^2;

rho_water = ones(nx,ny) * rho_water0;
rho_air = ones(nx,ny) * 1e-5 * rho_air0;
rho_water(mask) = 1e-5 * rho_water0;
rho_air(mask) = rho_air0;

f_water = zeros(nx,ny,npop);
f_air = zeros(nx,ny,npop);
for i = 1:npop
    f_water(:,:,i) = w(i) * rho_water;
    f_air(:,:,i) = w(i) * rho_air;
end

% loop
for t = 1:nsteps

    % computar densidades
    rho_water = sum(f_water,3);
    rho_air = sum(f_air,3);

    % velocidades
    ux = zeros(nx,ny);
    uy = zeros(nx,ny);
    for i = 1:npop
        ux = ux + (f_water(:,:,i) + f_air(:,:,i)) * cx(i);
        uy = uy + (f_water(:,:,i) + f_air(:,:,i)) * cy(i);
    end
    rho_total = rho_water + rho_air + 1e-9;
    ux = ux ./ rho_total;
    uy = uy ./ rho_total;

    psi_w = psi(rho_water,rho_water0);
    psi_a = psi(rho_air,rho_air0);

    % shan-chen
    Fx_w = zeros(nx,ny); Fy_w = zeros(nx,ny);
    Fx_a = zeros(nx,ny); Fy_a = zeros(nx,ny);
    for i = 2:npop
        psi_a_shift = circshift(psi_a, [-cy(i), -cx(i)]);
        Fx_w = Fx_w + w(i) * cx(i) * psi_a_shift;
        Fy_w = Fy_w + w(i) * cy(i) * psi_a_shift;

        psi_w_shift = circshift(psi_w, [-cy(i), -cx(i)]);
        Fx_a = Fx_a + w(i) * cx(i) * psi_w_shift;
        Fy_a = Fy_a + w(i) * cy(i) * psi_w_shift;
    end

    Fx_w = -G * psi_w .* Fx_w;
    Fy_w = -G * psi_w .* Fy_w;
    Fx_a = -G * psi_a .* Fx_a;
    Fy_a = -G * psi_a .* Fy_a;

    % no caso de voces, zera aqui
    %Fx_w(:) = 0; Fy_w(:) = 0;
    %Fx_a(:) = 0; Fy_a(:) = 0;

    u1w = ux + 0.5 * Fx_w ./ (rho_water + 1e-9);
    u2w = uy + 0.5 * Fy_w ./ (rho_water + 1e-9);
    u1a = ux + 0.5 * Fx_a ./ (rho_air + 1e-9);
    u2a = uy + 0.5 * Fy_a ./ (rho_air + 1e-9);

    % colisao regularizada c guo
    for i = 1:npop
        cuw = cx(i)*u1w + cy(i)*u2w;
        cua = cx(i)*u1a + cy(i)*u2a;

        usq_w = u1w.^2 + u2w.^2;
        usq_a = u1a.^2 + u2a.^2;

        feq_w = w(i) * rho_water .* (1 + 3*cuw + 4.5*cuw.^2 - 1.5*usq_w);
        feq_a = w(i) * rho_air   .* (1 + 3*cua + 4.5*cua.^2 - 1.5*usq_a);

        FdotCw = cx(i)*Fx_w + cy(i)*Fy_w;
        FdotCa = cx(i)*Fx_a + cy(i)*Fy_a;

        fguo_w = w(i) * (1 - 0.5 * omega_water) .* ...
                 (3*FdotCw + 9*cuw.*FdotCw);
        fguo_a = w(i) * (1 - 0.5 * omega_air) .* ...
                 (3*FdotCa + 9*cua.*FdotCa);

        f_neq_w = f_water(:,:,i) - feq_w;
        f_neq_a = f_air(:,:,i)   - feq_a;

        Pi_xx_w = cx(i)^2 * f_neq_w;
        Pi_xy_w = cx(i)*cy(i) * f_neq_w;
        Pi_yy_w = cy(i)^2 * f_neq_w;

        Pi_xx_a = cx(i)^2 * f_neq_a;
        Pi_xy_a = cx(i)*cy(i) * f_neq_a;
        Pi_yy_a = cy(i)^2 * f_neq_a;

        Hab_Pi_w = (cx(i)^2 - cs2) .* Pi_xx_w ...
                 + 2*cx(i)*cy(i).*Pi_xy_w ...
                 + (cy(i)^2 - cs2).*Pi_yy_w;

        Hab_Pi_a = (cx(i)^2 - cs2) .* Pi_xx_a ...
                 + 2*cx(i)*cy(i).*Pi_xy_a ...
                 + (cy(i)^2 - cs2).*Pi_yy_a;

        f_water(:,:,i) = feq_w + ...
            (1 - 1/tau_water) * w(i) * Hab_Pi_w / (2*cs2^2) + fguo_w;

        f_air(:,:,i) = feq_a + ...
            (1 - 1/tau_air) * w(i) * Hab_Pi_a / (2*cs2^2) + fguo_a;
    end

    % streaming
    for i = 1:npop
        f_water(:,:,i) = circshift(f_water(:,:,i),[cy(i),cx(i)]);
        f_air(:,:,i) = circshift(f_air(:,:,i),[cy(i),cx(i)]);
    end

    % visualizar
    if mod(t,output_interval) == 0
        imagesc(rho_air'); axis equal off; colorbar;
        title(['step ', num2str(t)]);
        drawnow;
    end
end
