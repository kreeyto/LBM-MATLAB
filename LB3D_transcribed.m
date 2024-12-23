%{ 
TO DO:
1. remove x flattening;
2. record sim with full res;
3. zero external forces;
4. test diff phase field velocity sets;
%}

% D3Q19 
clc; clearvars; close all

%% parameters

% phase field velocity set
pfvs = "D3Q15";

% vis mode  
slicebool = 1;

vid = 1;
res = 1;
stamp = 100;

nonzero = 1e-9; % 10^-9

tau = 0.505;
cssq = 1/3;
omega = 1/tau;
sharp_c = 0.15*3;
sigma = 0.1;

[nx, ny, nz] = deal(150*res);
radius = 20*res;

nsteps = 10000; 

%=============
bool_ind = 1;
boolfun = bool_ind .* ones(nx,ny,nz);
%=============

fr = 15; 

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

%% arrays and variables

[ux, uy, uz, ...
 phi, normx, normy, normz, ...
 curvature, indicator, ...
 ffx, ffy, ffz, ...
 mod_grad, isfluid] = deal(zeros(nx,ny,nz));

[rho, pxx, pyy, pzz, ...
      pxy, pxz, pyz] = deal(ones(nx,ny,nz));

w = zeros(1,fpoints);
w_g = zeros(1,gpoints);

fneq = zeros(fpoints,1); 
isfluid(2:nx-1,2:ny-1,2:nz-1) = 1;

%% velocity set properties

w(1) = 1/3;
w(2:7) = 1/18;
w(8:19) = 1/36;

cix = [0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0];
ciy = [0, 0, 0, 1, -1, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 1, -1];
ciz = [0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1];

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

%% phase field init

for i = 2:nx-1
    for j = 2:ny-1
        for k = 2:nz-1
            Ri = sqrt((i-nx/2)^2/2.^2 + (j-ny/2)^2 + (k-nz/2)^2);
            phi(i,j,k) = 0.5 + 0.5 * tanh(2*(radius-Ri)/(3*res));
        end
    end
end

%% distribution function init

for i = 1:fpoints
    f(:,:,:,i) = w(i) * rho(:,:,:);
end

for i = 1:gpoints
    g(:,:,:,i) = w_g(i) * phi(:,:,:);
end

%% if volshow

if slicebool ~= 1
    hVol = volshow(phi, 'RenderingStyle', 'Isosurface');
    viewer = hVol.Parent;
    hFig = viewer.Parent;
end

%% video

if vid == 1
    currentFile = mfilename('fullpath');
    [currentDir, ~, ~] = fileparts(currentFile);
    videoPath = fullfile(currentDir, 'LBM3D.avi'); 
    obj = VideoWriter(videoPath);
    obj.FrameRate = fr; 
    open(obj);
end

%% simulation loop

try
for t = 1:nsteps

    % phase field calc
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz    
                if isfluid(i,j,k) == 1
                    phi(i,j,k) = sum(g(i,j,k,:),4);
                end
            end
        end
    end

    % normal calculation and arrays
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                if isfluid(i,j,k) == 1
                    grad_fix = 0; grad_fiy = 0; grad_fiz = 0;
                    for l = 1:fpoints
                        grad_fix = grad_fix + 3 * w(l) .* cix(l) .* phi(i+cix(l),j+ciy(l),k+ciz(l));
                        grad_fiy = grad_fiy + 3 * w(l) .* ciy(l) .* phi(i+cix(l),j+ciy(l),k+ciz(l));
                        grad_fiz = grad_fiz + 3 * w(l) .* ciz(l) .* phi(i+cix(l),j+ciy(l),k+ciz(l));
                    end
                    mod_grad(i,j,k) = sqrt(grad_fix.^2 + grad_fiy.^2 + grad_fiz.^2);
                    normx(i,j,k) = grad_fix ./ (mod_grad(i,j,k) + nonzero);
                    normy(i,j,k) = grad_fiy ./ (mod_grad(i,j,k) + nonzero);
                    normz(i,j,k) = grad_fiz ./ (mod_grad(i,j,k) + nonzero);
                    indicator(i,j,k) = sqrt(grad_fix.^2 + grad_fiy.^2 + grad_fiz.^2);
                end
            end
        end
    end

    % curvature
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                if isfluid(i,j,k) == 1
                    % surface tension force
                    curvature(i,j,k) = 0;
                    for l = 1:fpoints
                        curvature(i,j,k) = curvature(i,j,k) - 3 .* w(l) .* ...
                        (cix(l) .* normx(i+cix(l),j+ciy(l),k+ciz(l)) + ...
                         ciy(l) .* normy(i+cix(l),j+ciy(l),k+ciz(l)) + ...
                         ciz(l) .* normz(i+cix(l),j+ciy(l),k+ciz(l)) ...
                        );
                    end
                    ffx(i,j,k) = sigma .* curvature(i,j,k) .* normx(i,j,k) .* indicator(i,j,k) .* boolfun(i,j,k);
                    ffy(i,j,k) = sigma .* curvature(i,j,k) .* normy(i,j,k) .* indicator(i,j,k) .* boolfun(i,j,k);
                    ffz(i,j,k) = sigma .* curvature(i,j,k) .* normz(i,j,k) .* indicator(i,j,k) .* boolfun(i,j,k);
                end
            end
        end
    end

    % momenti
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                if isfluid(i,j,k) == 1
                    ux(i,j,k) = ( ...
                        f(i,j,k,2) - f(i,j,k,3) + f(i,j,k,8) - f(i,j,k,9) + f(i,j,k,10) - f(i,j,k,11) + f(i,j,k,14) - f(i,j,k,15) + f(i,j,k,16) - f(i,j,k,17) ...
                    ) ./ rho(i,j,k) + ffx(i,j,k) * 0.5 ./ rho(i,j,k);
                    uy(i,j,k) = ( ...
                        f(i,j,k,4) - f(i,j,k,5) + f(i,j,k,8) - f(i,j,k,9) + f(i,j,k,12) - f(i,j,k,13) - f(i,j,k,14) + f(i,j,k,15) + f(i,j,k,18) - f(i,j,k,19) ...
                    ) ./ rho(i,j,k) + ffy(i,j,k) * 0.5 ./ rho(i,j,k);
                    uz(i,j,k) = ( ...
                        f(i,j,k,6) - f(i,j,k,7) + f(i,j,k,10) - f(i,j,k,11) + f(i,j,k,12) - f(i,j,k,13) - f(i,j,k,16) + f(i,j,k,17) - f(i,j,k,18) + f(i,j,k,19) ...
                    ) ./ rho(i,j,k) + ffz(i,j,k) * 0.5 ./ rho(i,j,k);
                    uu = 0.5 * (ux(i,j,k).^2 + uy(i,j,k).^2 + uz(i,j,k).^2) / cssq;
                    rho(i,j,k) = sum(f(i,j,k,:),4);
                    for l = 1:fpoints
                        udotc = (ux(i,j,k) * cix(l) + uy(i,j,k) * ciy(l) + uz(i,j,k) * ciz(l)) / cssq;
                        HeF = (w(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5.*udotc.^2 - uu))) ...
                            .* ((cix(l) - ux(i,j,k)) .* ffx(i,j,k) + ...
                                (ciy(l) - uy(i,j,k)) .* ffy(i,j,k) + ...
                                (ciz(l) - uz(i,j,k)) .* ffz(i,j,k) ...
                               ) ./ (rho(i,j,k) .* cssq);
                        feq = w(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5.*udotc.^2 - uu)) - 0.5 .* HeF;
                        fneq(l) = f(i,j,k,l) - feq;
                    end
                    pxx(i,j,k) = fneq(2) + fneq(3) + fneq(8) + fneq(9) + fneq(10) + fneq(11) + fneq(14) + fneq(15) + fneq(16) + fneq(17);
                    pyy(i,j,k) = fneq(4) + fneq(5) + fneq(8) + fneq(9) + fneq(12) + fneq(13) + fneq(14) + fneq(15) + fneq(18) + fneq(19);
                    pzz(i,j,k) = fneq(6) + fneq(7) + fneq(10) + fneq(11) + fneq(12) + fneq(13) + fneq(16) + fneq(17) + fneq(18) + fneq(19);
                    pxy(i,j,k) = fneq(8) + fneq(9) - fneq(14) - fneq(15);
                    pxz(i,j,k) = fneq(10) + fneq(11) - fneq(16) - fneq(17);
                    pyz(i,j,k) = fneq(12) + fneq(13) - fneq(18) - fneq(19);
                end
            end
        end
    end

    % collision
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                if isfluid(i,j,k) == 1
                    uu = 0.5 * (ux(i,j,k).^2 + uy(i,j,k).^2 + uz(i,j,k).^2) / cssq;
                    for l = 1:fpoints
                        udotc = (ux(i,j,k) * cix(l) + uy(i,j,k) * ciy(l) + uz(i,j,k) * ciz(l)) / cssq;
                        feq = w(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5.*udotc.^2 - uu));
                        HeF = 0.5 .* (w(l) * (rho(i,j,k) + rho(i,j,k) .* (udotc + 0.5.*udotc.^2 - uu))) .* ...
                            ((cix(l) - ux(i,j,k)) .* ffx(i,j,k) + ...
                             (ciy(l) - uy(i,j,k)) .* ffy(i,j,k) + ...
                             (ciz(l) - uz(i,j,k)) .* ffz(i,j,k) ...
                            ) ./ (rho(i,j,k) .* cssq);
                        fneq = (cix(l) .* cix(l) - cssq) * pxx(i,j,k) + ...
                               (ciy(l) .* ciy(l) - cssq) * pyy(i,j,k) + ...
                               (ciz(l) .* ciz(l) - cssq) * pzz(i,j,k) + ...
                               2 * cix(l) .* ciy(l) .* pxy(i,j,k) + ...
                               2 * cix(l) .* ciz(l) .* pxz(i,j,k) + ...
                               2 * ciy(l) .* ciz(l) .* pyz(i,j,k);
                        f(i+cix(l),j+ciy(l),k+ciz(l),l) = feq + (1-omega) * (w(l) / (2*cssq^2)) * fneq + HeF;
                    end
                    for l = 1:gpoints
                        udotc = (ux(i,j,k) * cix(l) + uy(i,j,k) * ciy(l) + uz(i,j,k) * ciz(l)) / cssq;
                        feq = w_g(l) .* phi(i,j,k) .* (1 + udotc);
                        Hi = sharp_c .* phi(i,j,k) .* (1 - phi(i,j,k)) .* ...
                            (cix(l) .* normx(i,j,k) + ...
                             ciy(l) .* normy(i,j,k) + ...
                             ciz(l) .* normz(i,j,k)) .* boolfun(i,j,k); 
                        g(i,j,k,l) = feq + w_g(l) .* Hi;
                    end
                end
            end
        end
    end
    
    for l = 1:gpoints
        g(:,:,:,l) = circshift(g(:,:,:,l),[cix(l),ciy(l),ciz(l),0]);
    end

    % boundary conditions
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                if isfluid(i,j,k) == 0
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
    end
    phi(:,:,1) = phi(:,:,2);      
    phi(:,:,nz) = phi(:,:,nz-1); 
    phi(:,1,:) = phi(:,2,:); 
    phi(:,ny,:) = phi(:,ny-1,:);
    phi(1,:,:) = phi(2,:,:);   
    phi(nx,:,:) = phi(nx-1,:,:);
    % periodic x

    if mod(t,stamp) == 0      
        if slicebool == 1
            x = 1:nx; y = 1:ny; z = 1:nz;
            h = slice(x,y,z,phi,nx/2,[],[]); 
            shading interp; colorbar; axis tight; 
            title(['t = ', num2str(t)]);
        else
            hVol.Data = phi; 
        end
        drawnow;
        if vid == 1
            frame = getframe(gcf);
            writeVideo(obj, frame);
        end
    end

    disp(['tstep = ', num2str(t)]);

end
catch ME
    disp('Erro durante a simulação:');
    disp(ME.message);
end

if vid == 1
    close(obj);
end


