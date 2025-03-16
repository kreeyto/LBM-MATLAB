clc; clearvars; close all

res = 1;

nx = 256*res;
ny = 128*res;

Re = 100;  

u_in = 0.2; 
diam = 10;     

nu = u_in * diam / Re; 
tau = (nu / (1/3)) + 0.5;
omega = 1 / tau; 
cssq = 1/3;

nsteps = 10000;
stamp = 100;     

fpoints = 9;

f = zeros(nx,ny,fpoints);
[ux, uy, rho] = deal(zeros(nx,ny));

w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

cix = [0, 1, 0, -1, 0,  1, -1, -1,  1];
ciy = [0, 0, 1,  0, -1, 1,  1, -1, -1];

ii = 2:nx-1; jj = 2:ny-1;

rho(:,:) = 1.0; 

nozzle_center = floor(ny/2);
nozzle_half = 5;

for j = (nozzle_center-nozzle_half):(nozzle_center+nozzle_half)
    rho(1,j) = 1.0;  
    ux(1,j) = u_in;
    uy(1,j) = 0.0;
end

for l = 1:fpoints
    f(:,:,l) = w(l) .* rho(:,:);
end

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

    for j = (nozzle_center-nozzle_half):(nozzle_center+nozzle_half)
        rho(1,j) = 1.0;     
        ux(1,j) = u_in;
        uy(1,j) = 0.0;
        
        uu = 0.5 * (u_in^2) / cssq;
        for k = 1:fpoints
            udotc = (ux(1,j) * cix(k) + uy(1,j) * ciy(k)) / cssq;
            feq = w(k) * (rho(1,j) + rho(1,j) * (udotc + 0.5*udotc^2 - uu));
            if (1+cix(k) > 1 && 1+cix(k) <= nx && j+ciy(k) > 0 && j+ciy(k) <= ny)
                f(1+cix(k),j+ciy(k),k) = feq;
            end
        end
    end  

    if (mod(t,stamp) == 0)      
        imagesc(ux'); 
        axis xy; axis equal; colorbar;
        title(['t = ',num2str(t)]);
        pause(0.01);
    end

    disp(t);
end
