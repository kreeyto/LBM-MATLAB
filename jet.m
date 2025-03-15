clc; clearvars; close all

res = 1;

tau = 0.505;
omega = 1/tau;
cssq = 1/3;
tau = 0.505;
sharp_c = 0.15*3;
sigma = 0.1;

nx = 128*res;
ny = 64*res;
radius = 20*res;

nsteps = 10000;
stamp = 10;

fpoints = 9;
gpoints = 9;

f = zeros(nx,ny,fpoints);
g = zeros(nx,ny,gpoints);

[ux, uy, ...
 phi, normx, normy, ...
 mod_grad] = deal(zeros(nx,ny));

[pxx, pyy, pxy, rho] = deal(ones(nx,ny));

w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
w_g = w;

cix = [0, 1, 0, -1, 0, 1, -1, -1, 1];
ciy = [0, 0, 1, 0, -1, 1, 1, -1, -1];

fneq = zeros(9,1);
ii = 2:nx-1; jj = 2:ny-1;


for l = 1:fpoints
    f(:,:,l) = w(l) * rho(:,:);
end
for l = 1:gpoints
    g(:,:,l) = w_g(l) * phi(:,:);
end

for t = 1:nsteps

    % phase field calc
    phi(ii,jj) = sum(g(ii,jj,:),3);

    % normal calculation and arrays
    for i = 2:nx-1
        for j = 2:ny-1
            grad_fix = 0;
            grad_fiy = 0;
            for k = 1:fpoints
                grad_fix = grad_fix + 3 * w(k) .* cix(k) .* phi(i+cix(k),j+ciy(k));
                grad_fiy = grad_fiy + 3 * w(k) .* ciy(k) .* phi(i+cix(k),j+ciy(k));
            end
            mod_grad(i,j) = sqrt(grad_fix.^2 + grad_fiy.^2);
            normx(i,j) = grad_fix ./ (mod_grad(i,j) + 1e-9);
            normy(i,j) = grad_fiy ./ (mod_grad(i,j) + 1e-9);                            
        end
    end

    % momenti
    for i = 2:nx-1
        for j = 2:ny-1
            ux(i,j) = ( ...
                f(i,j,2) - f(i,j,4) + f(i,j,6) - f(i,j,7) - f(i,j,8) + f(i,j,9) ...
            ) ./ rho(i,j);
            uy(i,j) = ( ...
                f(i,j,3) - f(i,j,5) + f(i,j,6) + f(i,j,7) - f(i,j,8) - f(i,j,9) ...
            ) ./ rho(i,j);   
            uu = 0.5* (ux(i,j).^2 + uy(i,j).^2) / cssq;
            rho(i,j) = sum(f(i,j,:),3); 
            for k = 1:fpoints
                udotc = (ux(i,j) * cix(k) + uy(i,j) * ciy(k)) / cssq;                    
                feq = w(k) * (rho(i,j) + rho(i,j) .* (udotc + 0.5 .* udotc.^2 - uu));
                fneq(k) = f(i,j,k) - feq;
            end
            pxx(i,j) = fneq(2) + fneq(4) + fneq(6) + fneq(7) + fneq(8) + fneq(9);
            pyy(i,j) = fneq(3) + fneq(5) + fneq(6) + fneq(7) + fneq(8) + fneq(9);
            pxy(i,j) = fneq(6) - fneq(7) + fneq(8) - fneq(9); 
        end 
    end

    % collision
    for i = 2:nx-1
        for j = 2:ny-1
             uu = 0.5 * (ux(i,j).^2 + uy(i,j).^2) / cssq;
             for k = 1:fpoints                                          
                 udotc = (ux(i,j) * cix(k) + uy(i,j) * ciy(k)) / cssq;
                 feq = w(k) * (rho(i,j) + rho(i,j) .* (udotc + 0.5 .* udotc.^2 - uu));
                 fneq = (w(k) / (2*cssq^2)) * ((cix(k) .* cix(k) - cssq) * pxx(i,j) + ...
                                               (ciy(k) .* ciy(k) - cssq) * pyy(i,j) + ...
                                                2 * cix(k) .* ciy(k) .* pxy(i,j));
                 f(i+cix(k),j+ciy(k),k) = feq + (1 - omega) * fneq;
             end
             for k = 1:gpoints
                 udotc = (ux(i,j) * cix(k) + uy(i,j) * ciy(k)) / cssq;
                 feq = w_g(k) .* phi(i,j) .* (1 + udotc);
                 Hi = sharp_c .* phi(i,j) .* (1 - phi(i,j)) .* (cix(k) .* normx(i,j) + ciy(k) .* normy(i,j));
                 g(i,j,k) = feq + w_g(k) .* Hi;
             end
        end
    end

    for k = 1:gpoints
        g(:,:,k) = circshift(g(:,:,k),[cix(k),ciy(k),0]);
    end

    % bcs
    for i = [1,nx]
        for j = [1,ny]
            for k=1:fpoints
                if(i+cix(k)>0 && j+ciy(k)>0)
                    f(i+cix(k),j+ciy(k),k) = rho(i,j) .* w(k);
                end
            end
            for k=1:gpoints
                if(i+cix(k)>0 && j+ciy(k)>0)
                    g(i+cix(k),j+ciy(k),k) = phi(i,j) .* w_g(k);
                end
            end
        end
    end
    phi(:,1) = phi(:,2);
    phi(:,ny) = phi(:,ny-1);
    phi(1,:) = phi(2,:);
    phi(nx,:) = phi(nx-1,:);

    if(mod(t,stamp) == 0)      
        imagesc(phi');
        axis xy;
        axis equal;
        colorbar;
        pause(0.01);
    end

    disp(t);
end