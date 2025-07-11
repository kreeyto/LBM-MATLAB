clc; clearvars; close all

% === diretórios de cada caso ===
base_paths = {
    fullfile(getenv('HOME'),'Desktop','MULTIC-TS-LBM','bin','D3Q19','100'),   % We = 100
    fullfile(getenv('HOME'),'Desktop','MULTIC-TS-LBM','bin','D3Q19','500'),   % We = 500
    fullfile(getenv('HOME'),'Desktop','MULTIC-TS-LBM','bin','D3Q19','2500')   % We = 2500
};

we_labels = {'$We=100$','$We=500$','$We=2500$'};
colors = {'r','g','b'};
z_target = 8 * 10;  % posição axial para todos os casos
U_inflow = 0.05;

% === dimensões do domínio ===
[NX,NY,NZ] = deal(64,64,256);
xc = round((NX-1)/2);
yc = round((NY-1)/2);

% === armazenamento de perfis e campos médios ===
r_profiles = cell(1,3);
uz_profiles = cell(1,3);
uz_planes = cell(1,3);  % para armazenar uz_mean(:,:,xc+1)

for k = 1:3
    files = dir(fullfile(base_paths{k}, '*_uz*.bin'));
    nfiles = numel(files);

    uz_stack = nan(NZ,NY,NX,nfiles);
    for f = 1:nfiles
        fname = fullfile(base_paths{k}, files(f).name);
        fid = fopen(fname,'rb');
        raw = fread(fid,'float32=>float32');
        fclose(fid);
        uz = reshape(raw,[NX,NY,NZ]);
        uz = permute(uz,[3,2,1]);  
        uz_stack(:,:,:,f) = uz;
    end

    uz_mean = mean(uz_stack,4,'omitnan');
    uz_mean = uz_mean / U_inflow;
    %uz_mean(uz_mean < 0) = 0;  
    uz_planes{k} = squeeze(uz_mean(:,:,xc+1));

    uz_profile = uz_mean(z_target+1,:,xc+1);  
    Uc = max(uz_profile);

    r = abs((0:NY-1)-yc);  
    r_over_z = r(yc+1:NY) / z_target;
    uz_half = uz_profile(yc+1:NY) / Uc;
    %valid = uz_half > 0;

    r_profiles{k} = r_over_z; % (valid)
    uz_profiles{k} = uz_half;
end

%% === três planos (z,y) lado a lado ===
figure('Name','time-averaged velocity planes')
for k = 1:3
    subplot(1,3,k)
    imagesc(uz_planes{k},[0 1])
    set(gca,'YDir','normal')
    colormap(hot)
    axis equal tight
    hold on
    yline(z_target+1,'w--','LineWidth',2.0)
    axis off
end

%% === plot sobreposto (self-similarity) ===
figure('Name','self-similarity','Position',[600 200 500 400]); hold on
for k = 1:3
    plot(r_profiles{k}, uz_profiles{k}, ...
        'LineWidth', 2, 'Color', colors{k}, 'DisplayName', we_labels{k});
end
xlabel('$r/z$','Interpreter','latex')
ylabel('$u_z/U_c$','Interpreter','latex')
legend('Interpreter','latex','Location','northeast')
grid on

output_dir = fullfile(pwd,'self_similarity_dat');
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

for k = 1:3
    data = [r_profiles{k}(:),uz_profiles{k}(:)];
    fname = fullfile(output_dir,sprintf('self_similarity_We%d.dat',str2double(regexp(we_labels{k},'\d+','match','once'))));
    writematrix(data, fname,'Delimiter','tab');
end