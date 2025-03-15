clc; clearvars; close all

res = 1;

nx = 128*res;
ny = 64*res;

%Re = 5000;  
We = 1000;

u_in = 0.09; 
diam = 10;     

%nu = u_in * L / Re; 
%tau = (nu / (1/3)) + 0.5;
tau = 0.505;
omega = 1 / tau; 
cssq = 1/3;
sharp_c = 0.15*3;  
sigma = u_in * u_in * (diam+diam) / We;

nsteps = 10000;
stamp = 100;     

fpoints = 9;
gpoints = 9;

f = zeros(nx,ny,fpoints);
g = zeros(nx,ny,gpoints);

[ux, uy, phi, normx, normy, mod_grad, curvature, ffx, ffy] = deal(zeros(nx,ny));
[pxx, pyy, pxy, rho]                                       = deal(ones(nx,ny));

fneq = zeros(fpoints,1);

w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
w_g = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

cix = [0, 1, 0, -1, 0,  1, -1, -1,  1];
ciy = [0, 0, 1,  0, -1, 1,  1, -1, -1];

ii = 2:nx-1; jj = 2:ny-1;

phi(:,:) = 0; 

nozzle_center = floor(ny/2);
nozzle_half = 5;

for j = (nozzle_center-nozzle_half):(nozzle_center+nozzle_half)
    phi(1,j) = 1.0;  
end

for l = 1:fpoints
    f(:,:,l) = w(l) .* rho(:,:);
end
for l = 1:gpoints
    g(:,:,l) = w_g(l) .* phi(:,:);
end

for t = 1:nsteps
    
    phi(ii,jj) = sum(g(ii,jj,:),3);
    
    for i = 2:nx-1
        for j = 2:ny-1
            grad_fix = 0; grad_fiy = 0; 
            for k = 1:fpoints
                grad_fix = grad_fix + 3 * w(k) * cix(k) * phi(i+cix(k),j+ciy(k));
                grad_fiy = grad_fiy + 3 * w(k) * ciy(k) * phi(i+cix(k),j+ciy(k));
            end
            mod_grad(i,j) = sqrt(grad_fix.^2 + grad_fiy.^2);
            normx(i,j) = grad_fix ./ (mod_grad(i,j) + 1e-9);
            normy(i,j) = grad_fiy ./ (mod_grad(i,j) + 1e-9);
        end
    end

    for i = 2:nx-1
        for j = 2:ny-1
            curvature(i,j) = 0;
            for k = 1:fpoints
                curvature(i,j) = curvature(i,j) - 3 .* w(k) .* (cix(k) .* (normx(i+cix(k),j+ciy(k))) + ...
                                                                ciy(k) .* (normy(i+cix(k),j+ciy(k))));
            end 
            ffx(i,j) = sigma .* curvature(i,j) .* normx(i,j) .* mod_grad(i,j);
            ffy(i,j) = sigma .* curvature(i,j) .* normy(i,j) .* mod_grad(i,j);
        end
    end

    for i = 2:nx-1
        for j = 2:ny-1

            rho(i,j) = sum(f(i,j,:),3);

            ux(i,j) = (f(i,j,2) - f(i,j,4) + f(i,j,6) - f(i,j,7) - f(i,j,8) + f(i,j,9)) ./ rho(i,j) + ffx(i,j) * 0.5 ./ rho(i,j);
            uy(i,j) = (f(i,j,3) - f(i,j,5) + f(i,j,6) + f(i,j,7) - f(i,j,8) - f(i,j,9)) ./ rho(i,j) + ffy(i,j) * 0.5 ./ rho(i,j);

            uu = 0.5 * (ux(i,j)^2 + uy(i,j)^2) / cssq;
            for k = 1:fpoints
                udotc = (ux(i,j) * cix(k) + uy(i,j) * ciy(k)) / cssq;
                HeF = (w(k) * (rho(i,j) + rho(i,j) .* (udotc + 0.5.*udotc.^2 - uu))) .* ...
                              ((cix(k) - ux(i,j)) .* ffx(i,j) + (ciy(k) - uy(i,j)) .* ffy(i,j)) ./ (rho(i,j) .* cssq);
                feq = w(k) * (rho(i,j) + rho(i,j) * (udotc + 0.5*udotc^2 - uu)) - 0.5 .* HeF;
                fneq(k) = f(i,j,k) - feq;  
            end
            pxx(i,j) = fneq(2) + fneq(4) + fneq(6) + fneq(7) + fneq(8) + fneq(9);
            pyy(i,j) = fneq(3) + fneq(5) + fneq(6) + fneq(7) + fneq(8) + fneq(9);
            pxy(i,j) = fneq(6) - fneq(7) + fneq(8) - fneq(9);
        end
    end

    for i = 2:nx-1
        for j = 2:ny-1
             uu = 0.5 * (ux(i,j)^2 + uy(i,j)^2) / cssq;
             for k = 1:fpoints
                 udotc = (ux(i,j) * cix(k) + uy(i,j) * ciy(k)) / cssq;
                 feq = w(k) * (rho(i,j) + rho(i,j) * (udotc + 0.5*udotc^2 - uu));
                 HeF = 0.5 .* (w(k) * (rho(i,j) + rho(i,j) .* (udotc + 0.5.*udotc.^2 - uu))) .* ...
                              ((cix(k) - ux(i,j)) .* ffx(i,j) + (ciy(k) - uy(i,j)) .* ffy(i,j)) ./ (rho(i,j) .* cssq);
                 fneqc = (w(k) / (2*cssq^2)) * ((cix(k)^2 - cssq) * pxx(i,j) + ...
                                                (ciy(k)^2 - cssq) * pyy(i,j) + ...
                                                 2 * cix(k) * ciy(k) * pxy(i,j));
                 f(i+cix(k),j+ciy(k),k) = feq + (1 - omega) * fneqc + HeF;
             end

             for k = 1:gpoints
                 udotc = (ux(i,j) * cix(k) + uy(i,j) * ciy(k)) / cssq;
                 geq = w_g(k) * phi(i,j) * (1 + udotc);
                 Hi = sharp_c * phi(i,j) * (1 - phi(i,j)) * (cix(k) * normx(i,j) + ciy(k) * normy(i,j));
                 g(i+cix(k),j+ciy(k),k) = geq + w_g(k) * Hi; %+ (1-omega_d).*(g(ii,jj,kk)-feq) ;
             end
        end
    end

    for j = (nozzle_center-nozzle_half):(nozzle_center+nozzle_half)
        rho(1,j) = 1.0;     
        phi(1,j) = 1.0;     
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
        
        for k = 1:gpoints
            udotc = (ux(1,j) * cix(k) + uy(1,j) * ciy(k)) / cssq;
            feq_g = w_g(k) * phi(1,j) * (1 + udotc);
            Hi = sharp_c * phi(1,j) * (1 - phi(1,j)) * (cix(k) * normx(1,j) + ciy(k) * normy(1,j));
            g(1,j,k) = feq_g + w_g(k) * Hi;
            if (1+cix(k) > 1 && 1+cix(k) <= nx && j+ciy(k) > 0 && j+ciy(k) <= ny)
                g(1+cix(k),j+ciy(k),k) = g(1,j,k);
            end
        end
    end  

    if (mod(t,stamp) == 0)      
        %quiver(ux', uy');
        imagesc(phi'); 
        axis xy; axis equal; colorbar;
        title(['t = ',num2str(t)]);
        pause(0.01);
    end

    disp(t);
end
