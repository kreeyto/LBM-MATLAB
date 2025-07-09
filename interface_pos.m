clc; clearvars; close all

base_path = fullfile(getenv('HOME'),'Desktop','MULTIC-TS-LBM','bin','D3Q19','577');

[NX,NY,NZ] = deal(64);

xc = round((NX-1)/2);
yc = round((NY-1)/2);
zc = round((NZ-1)/2);

files = dir(fullfile(base_path,'577_phi*.bin'));
nfiles = numel(files);
int_pos = nan(nfiles,1);
timesteps = nan(nfiles,1);
for f = 1:nfiles
    fname = fullfile(base_path, files(f).name);
    timesteps(f) = str2double(files(f).name(9:14));
    fid = fopen(fname,'rb');
    raw = fread(fid,'float32=>float32');
    fclose(fid);
    phi = reshape(raw,[NX,NY,NZ]);
    phi = permute(phi,[3,2,1]);  
    phi_line = squeeze(phi(xc+1,yc+1,:));  
    xvals = 0:NX-1;
    idx = find(phi_line >= 0.5,1,'first');
    if ~isempty(idx) && idx > 1
        x1 = xvals(idx-1);
        x2 = xvals(idx);
        p1 = phi_line(idx-1);
        p2 = phi_line(idx);
        alpha = (0.5-p1) / (p2-p1);
        int_pos(f) = x1 + alpha * (x2-x1);
    else
        int_pos(f) = NaN;
    end
end

[t_sorted, idx] = sort(timesteps);
z_interface_sorted = xc - int_pos(idx);

plot(t_sorted, z_interface_sorted, '-o', 'LineWidth', 1.5);
grid on;

%% export to tikz
export_data = [t_sorted(:), z_interface_sorted(:)];
writematrix(export_data, 'interface_position_vs_time.dat', 'Delimiter', 'tab');

%% findpeaks
y = z_interface_sorted;
t = t_sorted;
[~, locs] = findpeaks(y, t, 'MinPeakProminence', 0.05);
T_est = mean(diff(locs));
fprintf('periodo estimado findpeaks: %.3f tsteps\n', T_est);

