clc; clearvars; close all

nx = 128;
ny = 64;

u_in = 0.2; 
tau = 0.6;
omega = 1 / tau; 
cssq = 1/3;

nsteps = 1000;
stamp = 25;     

fpoints = 9;

f = zeros(nx,ny,fpoints);
[ux, uy, rho] = deal(zeros(nx,ny));

w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

cix = [0, 1, 0, -1, 0,  1, -1, -1,  1];
ciy = [0, 0, 1,  0, -1, 1,  1, -1, -1];

rho(:,:) = 1.0; 

ux(1,:) = u_in;  
uy(1,:) = 0.0;

for j = 1:ny
    rho(1,j) = 1.0; 
end

for l = 1:fpoints
    f(:,:,l) = w(l) .* rho(:,:);
end

nozzle_radius = ny / 10;  
j_center = ny / 2;
j_start = round(j_center - nozzle_radius);
j_end   = round(j_center + nozzle_radius);

for t = 1:nsteps

    for i = 2:nx-1
        for j = 2:ny-1
            rho(i,j) = sum(f(i,j,:),3);

            ux(i,j) = (f(i,j,2) - f(i,j,4) + f(i,j,6) - f(i,j,7) - f(i,j,8) + f(i,j,9)) ./ rho(i,j);
            uy(i,j) = (f(i,j,3) - f(i,j,5) + f(i,j,6) + f(i,j,7) - f(i,j,8) - f(i,j,9)) ./ rho(i,j);

            uu = 0.5 * (ux(i,j)^2 + uy(i,j)^2) / cssq;
            for k = 1:fpoints
                udotc = (ux(i,j) * cix(k) + uy(i,j) * ciy(k)) / cssq;
                feq = w(k) * (rho(i,j) + rho(i,j) * (udotc + 0.5*udotc^2 - uu));
                f(i+cix(k),j+ciy(k),k) = feq + (1 - omega) * (f(i,j,k) - feq);
            end
        end
    end

    for j = 1:ny
        rho(1,j) = 1.0;

        if j >= j_start && j <= j_end
            ux(1,j) = u_in;
        else
            ux(1,j) = 0.0;
        end
        uy(1,j) = 0.0;

        uu = 0.5 * (ux(1,j)^2 + uy(1,j)^2) / cssq;
        for k = 1:fpoints
            udotc = (ux(1,j) * cix(k) + uy(1,j) * ciy(k)) / cssq;
            feq = w(k) * (rho(1,j) + rho(1,j) * (udotc + 0.5*udotc^2 - uu));
            ix = 1 + cix(k);
            iy = j + ciy(k);
            if ix >= 1 && ix <= nx && iy >= 1 && iy <= ny
                f(ix, iy, k) = feq;
            end
        end
    end

    if mod(t, stamp) == 0  
        imagesc(ux'); 
        axis xy; axis equal tight; colorbar;
        colormap(jet);
        drawnow;
    end

    disp(['passo: ', num2str(t)]);
end

