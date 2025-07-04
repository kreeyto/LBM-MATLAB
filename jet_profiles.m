clc; clearvars; close all

base_path = fullfile(getenv('HOME'),'Desktop','MULTIC-TS-LBM','bin','D3Q19','32113');
files = dir(fullfile(base_path, '32113_uz*.bin'));
nfiles = numel(files);

[NX,NY,NZ] = deal(64,64,256);
xc = round((NX-1)/2);
yc = round((NY-1)/2);

uz_stack = nan(NZ,NY,NX,nfiles);
for f = 1:nfiles
    fname = fullfile(base_path, files(f).name);
    fid = fopen(fname,'rb');
    raw = fread(fid,'float32=>float32');
    fclose(fid);
    uz = reshape(raw,[NX,NY,NZ]);
    uz = permute(uz,[3,2,1]);  
    uz_stack(:,:,:,f) = uz;
end

uz_mean = mean(uz_stack, 4, 'omitnan');

%figure
%uz_plot = squeeze(uz_mean(:,:,xc+1)) / 0.05;
%imagesc(uz_plot);
%axis equal tight
%set(gca, 'YDir', 'normal');  

z_target = min(128,NZ-1); % twaek

uz_profile = uz_mean(z_target+1,:,xc+1);  
Uc = max(uz_profile);  
r = abs((0:NY-1)-yc);  
r_over_z = r / z_target;
uz_norm = uz_profile / Uc;

figure;
%plot(r_over_z,uz_norm,'-o','MarkerFaceColor','r');
half = yc+1:NY;  
r_over_z_half = r(half) / z_target;
uz_half = uz_profile(half) / Uc;
plot(r_over_z_half, uz_half, 'ko', 'MarkerFaceColor','r');
hold on;
%f_theory = @(r) exp(-75 * (r).^2);  
%r_fit = linspace(0,max(r_over_z),200);
%plot(r_fit,f_theory(r_fit),'b-','LineWidth',1.5,'DisplayName','$\exp(-7.5(r/z)^2)$');
grid on;
% export_centerline = [z_over_D(:), uz_norm(:)];
% writematrix(export_centerline, 'centerline_decay_32113.dat', 'Delimiter', 'tab');
