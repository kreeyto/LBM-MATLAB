clc; clearvars; close all

res = 0.5;

nx = 256*res;
ny = 64*res;
diam = ny/10;  

u_in = 0.05; 
Re = 5000;  
We = 500;

u_in2 = 0.01;
Re2 = 1200;
diam2 = 10;

visc1 = (u_in * diam) / Re;
visc2 = (u_in2 * diam2) / Re2;
tau = 0.5 + 3 * visc1;
omega = 1 / tau; 
cssq = 1/3;
sharp_c = 0.15*7;  
sigma = (u_in * u_in * diam) / We;

nsteps = 10000;
stamp = 100;     

fpoints = 9;
gpoints = 5;

f = zeros(nx,ny,fpoints);
g = zeros(nx,ny,gpoints);

[ux, uy, phi, normx, normy, mod_grad, ffx, ffy, rho] = deal(zeros(nx,ny));
[pxx, pyy, pxy]                                      = deal(ones(nx,ny));

fneq = zeros(fpoints,1);

w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
w_g = [2/6, 1/6, 1/6, 1/6, 1/6];

cix = [0, 1, 0, -1, 0,  1, -1, -1,  1];
ciy = [0, 0, 1,  0, -1, 1,  1, -1, -1];

ii = 2:nx-1; jj = 2:ny-1;

isfluid = zeros(nx,ny);
isfluid(ii,jj) = 1;

rho(:,:) = 1;
phi(:,:) = 0;

rho1 = 1.0;
rho2 = 1.0;

nozzle_center = floor(ny/2);
nozzle_half = round(diam*res);

for l = 1:fpoints
    f(:,:,l) = w(l) .* rho(:,:);
end
for l = 1:gpoints
    g(:,:,l) = w_g(l) .* phi(:,:);
end

for t = 1:nsteps

    for j = (nozzle_center-nozzle_half):(nozzle_center+nozzle_half)
        rho(1,j) = rho1;     
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
            if (1+cix(k) > 1 && 1+cix(k) <= nx && j+ciy(k) > 0 && j+ciy(k) <= ny)
                g(1+cix(k),j+ciy(k),k) = feq_g;
            end
        end
    end  
    
    phi(ii,jj) = sum(g(ii,jj,:),3);
    
    for i = 2:nx-1
        for j = 2:ny-1
            grad_fix = 0; grad_fiy = 0; 
            for k = 1:gpoints
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
            curvature = 0;
            for k = 1:gpoints
                curvature = curvature - 3 .* w(k) .* (cix(k) .* (normx(i+cix(k),j+ciy(k))) + ...
                                                      ciy(k) .* (normy(i+cix(k),j+ciy(k))));
            end 
            ffx(i,j) = sigma .* curvature .* normx(i,j) .* mod_grad(i,j);
            ffy(i,j) = sigma .* curvature .* normy(i,j) .* mod_grad(i,j);
        end
    end

    for i = 2:nx-1
        for j = 2:ny-1

            rho(i,j) = sum(f(i,j,:),3);
            %rho(i,j) = phi(i,j) * rho1 + (1 - phi(i,j)) * rho2;

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
             visc_loc = phi(i,j) .* visc1 + (1-phi(i,j)) .* visc2;
             omega = 1 / (visc_loc*3 + 0.5);
             for k = 1:fpoints
                 udotc = (ux(i,j) * cix(k) + uy(i,j) * ciy(k)) / cssq;
                 feq = w(k) * (rho(i,j) + rho(i,j) * (udotc + 0.5*udotc^2 - uu));
                 HeF = 0.5 .* (w(k) * (rho(i,j) + rho(i,j) .* (udotc + 0.5.*udotc.^2 - uu))) .* ...
                              ((cix(k) - ux(i,j)) .* ffx(i,j) + (ciy(k) - uy(i,j)) .* ffy(i,j)) ./ (rho(i,j) .* cssq);
                 fneq_reg = (w(k) / (2*cssq^2)) * ((cix(k)^2 - cssq) * pxx(i,j) + ...
                                                (ciy(k)^2 - cssq) * pyy(i,j) + ...
                                                 2 * cix(k) * ciy(k) * pxy(i,j));
                 f(i+cix(k),j+ciy(k),k) = feq + (1 - omega) * fneq_reg + HeF;
             end

             for k = 1:gpoints
                 udotc = (ux(i,j) * cix(k) + uy(i,j) * ciy(k)) / cssq;
                 geq = w_g(k) * phi(i,j) * (1 + udotc);
                 Hi = sharp_c * phi(i,j) * (1 - phi(i,j)) * (cix(k) * normx(i,j) + ciy(k) * normy(i,j));
                 g(i+cix(k),j+ciy(k),k) = geq + w_g(k) * Hi; %+ (1-omega_d).*(g(ii,jj,kk)-feq) ;
             end
        end
    end

    for i = 1:nx
        for j = 1:ny
            if (isfluid(i,j) == 0)
                uu = 0.5 * (ux(i,j)^2 + uy(i,j)^2) / cssq;
                for k = 1:fpoints
                    udotc = (ux(i,j) * cix(k) + uy(i,j) * ciy(k)) / cssq;
                    feq = w(k) * (rho(i,j) + rho(i,j) * (udotc + 0.5*udotc^2 - uu));
                    if (i+cix(k) >= 1 && i+cix(k) <= nx && ...
                        j+ciy(k) >= 1 && j+ciy(k) <= ny)
                        fneq_reg = (w(k) / (2*cssq^2)) * ((cix(k)^2 - cssq) * pxx(i+cix(k),j+ciy(k)) + ...
                                                          (ciy(k)^2 - cssq) * pyy(i+cix(k),j+ciy(k)) + ...
                                                        2 * cix(k) * ciy(k) * pxy(i+cix(k),j+ciy(k)));
                        f(i,j,k) = feq + (1-omega) * fneq_reg;
                    end
                end
                for k = 1:gpoints 
                    udotc = (ux(i,j) * cix(k) + uy(i,j) * ciy(k)) / cssq;
                    geq = w_g(k) * phi(i,j) * (1 + udotc);
                    if (i+cix(k)>0 && j+ciy(k)>0)
                        g(i+cix(k),j+ciy(k),k) = geq;
                    end
                end
            end
        end
    end
    
    phi(:,1) = phi(:,ny-1); 
    phi(:,ny) = phi(:,2);
    phi(nx,:) = phi(nx-1,:);

    if (mod(t,stamp) == 0)      
        imagesc(phi); 
        axis xy; axis equal; colorbar;
        title(['t = ',num2str(t)]);
        pause(0.01);
    end

    disp(t);
end
