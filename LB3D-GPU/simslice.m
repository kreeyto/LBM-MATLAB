% D3Q19 
clc; clearvars; close all

%% parameters

% phase field velocity set
pfvs = "D3Q15";

tau = 0.505;
cssq = 1/3;
omega = 1/tau;
sharp_c = 0.15*3;
sigma = 0.1;
stamp = 1;

radius = 20;
res = 1;

[nx, ny, nz] = deal(150*res);
nsteps = 10000; 

fpoints = 19; 
if pfvs == "D3Q19"
    gpoints = 19;
elseif pfvs == "D3Q15"
    gpoints = 15;
elseif pfvs == "D3Q7"
    gpoints = 7;
end
f = zeros(nx,ny,nz,fpoints); 
g = zeros(nx,ny,nz,gpoints); 

%% index pre-alloc

ix = round(nx/2); iy = 2:ny-1; iz = 2:nz-1;
phi_slice = zeros(length(iy), length(iz));

%% arrays and variables

[rho, ux, uy, uz, ...
 phi, normx, normy, normz, ...
 curvature, indicator, ...
 ffx, ffy, ffz, ...
 mod_grad, isfluid] = deal(zeros(nx,ny,nz));

[pxx, pyy, pzz, ...
 pxy, pxz, pyz] = deal(ones(nx,ny,nz));

w = zeros(1,fpoints);
w_g = zeros(1,gpoints);

fneq = zeros(fpoints,1,1); 

isfluid(ix,iy,iz) = 1;
rho(:,:,:) = 1;

%% velocity set properties

w(1) = 1/3;
w(2:6) = 1/18;
w(7:19) = 1/36;

if pfvs == "D3Q19"
    w_g = w;
elseif pfvs == "D3Q15"
    w_g(1) = 2/9;
    w_g(2:7) = 1/9;
    w_g(8:15) = 1/72;
elseif pfvs == "D3Q7"
    w_g(1) = 1/4;
    w_g(2:7) = 1/8;
end

% opp = [1, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16, 19, 18];
cix = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0];
ciy = [0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1];
ciz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1];

%% phase field init

nx2 = nx/2; ny2 = ny/2; nz2 = nz/2; 
for j = iy
    for k = iz
        Ri = sqrt((ix-nx2)^2 + (j-ny2)^2/2.^2 + (k-nz2)^2);
        phi(ix,j,k) = 0.5 + 0.5 * tanh(2*(radius*res-Ri)/(3*res));
    end
end

%% distribution function init

for i = 1:fpoints
    f(:,:,:,i) = w(i) * rho(:,:,:);
end

for i = 1:gpoints
    g(:,:,:,i) = w_g(i) * phi(:,:,:);
end

%% simulation loop

for t = 1:nsteps

    % phase field calc
    phi(ix,iy,iz) = sum(g(ix,iy,iz,:),4);

    % normal and arrays
    [grad_fix, grad_fiy, grad_fiz] = deal(0);
    for l = 1:fpoints
        grad_fix = grad_fix + 3 * w(l) .* cix(l) .* ((phi((ix)+cix(l),(iy)+ciy(l),(iz)+ciz(l))));
        grad_fiy = grad_fiy + 3 * w(l) .* ciy(l) .* ((phi((ix)+cix(l),(iy)+ciy(l),(iz)+ciz(l))));
        grad_fiz = grad_fiz + 3 * w(l) .* ciz(l) .* ((phi((ix)+cix(l),(iy)+ciy(l),(iz)+ciz(l))));
    end
    mod_grad(ix,iy,iz) = sqrt(grad_fix.^2 + grad_fiy.^2 + grad_fiz.^2);
    normx(ix,iy,iz) = grad_fix ./ (mod_grad(ix,iy,iz) + 1e-9);
    normy(ix,iy,iz) = grad_fiy ./ (mod_grad(ix,iy,iz) + 1e-9);
    normz(ix,iy,iz) = grad_fiz ./ (mod_grad(ix,iy,iz) + 1e-9);
    indicator(ix,iy,iz) = sqrt(grad_fix.^2 + grad_fiy.^2 + grad_fiz.^2);

    % curvature
    curvature(ix,iy,iz) = 0;
    for l = 1:fpoints
        curvature(ix,iy,iz) = curvature(ix,iy,iz) - 3 .* w(l) .* ...
        (cix(l) .* (normx((ix)+cix(l),(iy)+ciy(l),(iz)+ciz(l))) + ...
         ciy(l) .* (normy((ix)+cix(l),(iy)+ciy(l),(iz)+ciz(l))) + ...
         ciz(l) .* (normz((ix)+cix(l),(iy)+ciy(l),(iz)+ciz(l))) ...
        );
    end
    ffx(ix,iy,iz) = sigma .* curvature(ix,iy,iz) .* normx(ix,iy,iz) .* indicator(ix,iy,iz);
    ffy(ix,iy,iz) = sigma .* curvature(ix,iy,iz) .* normy(ix,iy,iz) .* indicator(ix,iy,iz);
    ffz(ix,iy,iz) = sigma .* curvature(ix,iy,iz) .* normz(ix,iy,iz) .* indicator(ix,iy,iz);

    % moments
    for j = iy
        for k = iz
            ux(ix,j,k) = ( sum(f(ix,j,k,[2,16,10,8,14])) - sum(f(ix,j,k,[3,11,17,15,8])) ) ./ rho(ix,j,k) + ffx(ix,j,k) * 0.5 ./ rho(ix,j,k);
            uy(ix,j,k) = ( sum(f(ix,j,k,[4,8,15,18,12])) - sum(f(ix,j,k,[5,14,9,13,19])) ) ./ rho(ix,j,k) + ffx(ix,j,k) * 0.5 ./ rho(ix,j,k);
            uz(ix,j,k) = ( sum(f(ix,j,k,[7,16,11,18,13])) - sum(f(ix,j,k,[6,10,17,12,19])) ) ./ rho(ix,j,k) + ffx(ix,j,k) * 0.5 ./ rho(ix,j,k);
            uu = 0.5 * (ux(ix,j,k).^2 + uy(ix,j,k).^2 + uz(ix,j,k).^2) / cssq;
            rho(ix,j,k) = sum(f(ix,j,k,:),4);
            for l = 1:fpoints
                udotc = (ux(ix,j,k) * cix(l) + uy(ix,j,k) * ciy(l) + uz(ix,j,k) * ciz(l)) / cssq;
                HeF = (w(l) * (rho(ix,j,k) + rho(ix,j,k) .* (udotc + 0.5.*udotc.^2 - uu))) ...
                    .* ((cix(l) - ux(ix,j,k)) .* ffx(ix,j,k) + ...
                        (ciy(l) - uy(ix,j,k)) .* ffy(ix,j,k) + ...
                        (ciz(l) - uz(ix,j,k)) .* ffz(ix,j,k) ...
                       ) ./ (rho(ix,j,k) .* cssq);
                feq = w(l) * (rho(ix,j,k) + rho(ix,j,k) .* (udotc + 0.5.*udotc.^2 - uu)) - 0.5 .* HeF;
                fneq(l) = f(ix,j,k,l) - feq;
            end
            pxx(ix,j,k) = sum(fneq([2,3,8,9,10,11,14,15,16,17]));
            pyy(ix,j,k) = sum(fneq([4,5,8,9,12,13,14,15,18,19]));
            pzz(ix,j,k) = sum(fneq([6,7,10,11,12,13,16,17,18,19]));
            pxy(ix,j,k) = fneq(8) + fneq(9) - fneq(14) - fneq(15);
            pxz(ix,j,k) = fneq(10) + fneq(11) - fneq(16) - fneq(17);
            pyz(ix,j,k) = fneq(12) + fneq(13) - fneq(18) - fneq(19);
        end
    end

    % collision
    for j = iy
        for k = iz
            uu = 0.5 * (ux(ix,j,k).^2 + uy(ix,j,k).^2 + uz(ix,j,k).^2) / cssq;
            for l = 1:fpoints
                udotc = (ux(ix,j,k) * cix(l) + uy(ix,j,k) * ciy(l) + uz(ix,j,k) * ciz(l)) / cssq;
                feq = w(l) * (rho(ix,j,k) + rho(ix,j,k) .* (udotc + 0.5.*udotc.^2 - uu));
                HeF = 0.5 .* (w(l) * (rho(ix,j,k) + rho(ix,j,k) .* (udotc + 0.5.*udotc.^2 - uu))) .* ...
                    ((cix(l) - ux(ix,j,k)) .* ffx(ix,j,k) + ...
                     (ciy(l) - uy(ix,j,k)) .* ffy(ix,j,k) + ...
                     (ciz(l) - uz(ix,j,k)) .* ffz(ix,j,k) ...
                    ) ./ (rho(ix,j,k) .* cssq);
                fneq = (cix(l) .* cix(l) - cssq) * pxx(ix,j,k) + ...
                       (ciy(l) .* ciy(l) - cssq) * pyy(ix,j,k) + ...
                       (ciz(l) .* ciz(l) - cssq) * pzz(ix,j,k) + ...
                       2 * cix(l) .* ciy(l) .* pxy(ix,j,k) + ...
                       2 * cix(l) .* ciz(l) .* pxz(ix,j,k) + ...
                       2 * ciy(l) .* ciz(l) .* pyz(ix,j,k);
                f(ix+cix(l),j+ciy(l),k+ciz(l),l) = feq + (1-omega) * (w(l) / (2*cssq^2)) * fneq + HeF;
            end
            for l = 1:gpoints
                udotc = (ux(ix,j,k) * cix(l) + uy(ix,j,k) * ciy(l) + uz(ix,j,k) * ciz(l)) / cssq;
                feq = w_g(l) .* phi(ix,j,k) .* (1 + udotc);
                Hi = sharp_c .* phi(ix,j,k) .* (1 - phi(ix,j,k)) .* (cix(l) .* normx(ix,j,k) + ciy(l) .* normy(ix,j,k) + ciz(l) .* normz(ix,j,k)); 
                g(ix,j,k,l) = feq + w_g(l) .* Hi;
            end
        end
    end

    for l = 1:gpoints
        g(:,:,:,l) = circshift(g(:,:,:,l),[cix(l),ciy(l),ciz(l)]);
    end

    % boundary conditions
    for i = [1,nx]
        for j = [1,ny]
            for k = [1,nz]
                for l = 1:fpoints
                    if (i+cix(l)>0 && j+ciy(l)>0 && k+ciz(l)>0)
                        f(i+cix(l),j+ciy(l),k+ciz(l),l) = rho(i,j,k) .* w(l); 
                    end
                end
                for l = 1:gpoints
                    if (i+cix(l)>0 && j+ciy(l)>0 && k+ciz(l)>0)
                        g(i+cix(l),j+ciy(l),k+ciz(l),l) = phi(i,j,k) .* w_g(l);
                    end
                end
            end
        end
    end

    phi(:,:,1) = phi(:,:,2);  
    phi(:,:,nz) = phi(:,:,nz-1); 
    phi(1,:,:) = phi(2,:,:); 
    phi(nx,:,:) = phi(nx-1,:,:); 
    phi(:,1,:) = phi(:,2,:); 
    phi(:,ny,:) = phi(:,ny-1,:); 

    if(mod(t,stamp) == 0)      
        phi_slice(:,:) = phi(ix, iy, iz);
        imagesc(iz, iy, phi_slice'); colorbar;
        title(['t = ', num2str(t)]);
        xlabel('$z$'); ylabel('$y$');
        axis equal tight; 
        drawnow;
    end

    disp(['tstep = ', num2str(t)]);

end




