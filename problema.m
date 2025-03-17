clc; clearvars; close all

res = 1;

nx = 128*res;
ny = 64*res;

u_in = 0.2; 

tau = 0.55;
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

x_c = floor(nx/5);  
y_c = floor(ny/2);  
r   = 5; 

obstaculo = false(nx,ny);
for i = 1:nx
    for j = 1:ny
        dist = sqrt((i - x_c)^2 + (j - y_c)^2);
        if dist <= r 
            obstaculo(i,j) = true;
        end
    end
end

ux(1,:) = u_in;  
uy(1,:) = 0.0;

for j = 1:ny
    rho(1,j) = 1.0; 
end

for l = 1:fpoints
    f(:,:,l) = w(l) .* rho(:,:);
end

dx = 1; dy = 1;
vorticity_history = zeros(nx, ny, floor(nsteps/stamp));  
frame_index = 1;  

video_name = 'vorticity_simulation.mp4';
video_writer = VideoWriter(video_name, 'MPEG-4');
video_writer.FrameRate = 10; % Ajuste a taxa de quadros
open(video_writer);

for t = 1:nsteps

    for i = 2:nx-1
        for j = 2:ny-1
            if ~obstaculo(i,j)  
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
    end

    for i = 2:nx-1
        for j = 2:ny-1
            if obstaculo(i,j)
                f(i,j,2) = f(i,j,4);
                f(i,j,4) = f(i,j,2);
                f(i,j,3) = f(i,j,5);
                f(i,j,5) = f(i,j,3);
                f(i,j,6) = f(i,j,8);
                f(i,j,8) = f(i,j,6);
                f(i,j,7) = f(i,j,9);
                f(i,j,9) = f(i,j,7);
            end
        end
    end

    for j = 1:ny
        rho(1,j) = 1.0;     
        ux(1,j) = u_in;
        uy(1,j) = 0.0;
        
        uu = 0.5 * (u_in^2) / cssq;
        for k = 1:fpoints
            udotc = (ux(1,j) * cix(k) + uy(1,j) * ciy(k)) / cssq;
            feq = w(k) * (rho(1,j) + rho(1,j) * (udotc + 0.5*udotc^2 - uu));
            if (1+cix(k) > 1 && 1+cix(k) <= nx && j+ciy(k) > 0 && j+ciy(k) <= ny)
                f(1+cix(k), j+ciy(k), k) = feq;
            end
        end
    end   

    % Cálculo da vorticidade (curl)
    dudy = (ux(:,[2:end end]) - ux(:,[1 1:end-1])) / (2*dy);
    dvdx = (uy([2:end end],:) - uy([1 1:end-1],:)) / (2*dx);
    vorticity = dvdx - dudy;

    if mod(t, stamp) == 0  
        vorticity_history(:,:,frame_index) = vorticity;
        frame_index = frame_index + 1;
    end

    disp(['Passo: ', num2str(t)]);
end

for k = 1:frame_index-1
    imagesc(vorticity_history(:,:,k)'); 
    axis xy; axis equal; colorbar;
    colormap(jet);
    title(['Vorticidade no passo ', num2str(k * stamp)]);
    frame = getframe(gcf); 
    writeVideo(video_writer, frame); 
    pause(0.01);
end

close(video_writer);
disp(['Vídeo salvo como ', video_name]);
