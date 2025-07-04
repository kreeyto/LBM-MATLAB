clc; clearvars; close all

base_path = fullfile(getenv('HOME'),'Desktop','MULTIC-TS-LBM','bin','D3Q19','577');

[NX,NY,NZ] = deal(64);

xc = round((NX-1)/2);
yc = round((NY-1)/2);
zc = round((NZ-1)/2);

files = dir(fullfile(base_path, '577_phi*.bin'));
nfiles = numel(files);
int_pos = nan(nfiles, 1);
timesteps = nan(nfiles, 1);
for f = 1:nfiles
    fname = fullfile(base_path, files(f).name);
    timesteps(f) = str2double(files(f).name(9:14));
    fid = fopen(fname,'rb');
    raw = fread(fid,'float32=>float32');
    fclose(fid);
    phi = reshape(raw,[NX,NY,NZ]);
    phi = permute(phi,[3,2,1]);  
    phi_line = squeeze(phi(:,yc+1,xc+1));  
    zvals = 0:NZ-1;
    idx = find(phi_line >= 0.5,1,'first');
    if ~isempty(idx) && idx > 1
        z1 = zvals(idx-1);
        z2 = zvals(idx);
        p1 = phi_line(idx-1);
        p2 = phi_line(idx);
        alpha = (0.5-p1) / (p2-p1);
        int_pos(f) = z1 + alpha * (z2-z1);
    else
        int_pos(f) = NaN;
    end
end

[t_sorted, idx] = sort(timesteps);
z_interface_sorted = int_pos(idx)-zc;

plot(t_sorted, z_interface_sorted, '-o', 'LineWidth', 1.5);
grid on;

%% export to tikz
%export_data = [t_sorted(:), z_interface_sorted(:)];
%writematrix(export_data, 'interface_position_vs_time.dat', 'Delimiter', 'tab');

%% findpeaks
y = z_interface_sorted;
t = t_sorted;
[~, locs] = findpeaks(y, t, 'MinPeakProminence', 0.05);
T_est = mean(diff(locs));
fprintf('periodo estimado findpeaks: %.3f tsteps\n', T_est);

%% fft
signal = z_interface_sorted(:);     
time = t_sorted(:);                 
signal = signal - mean(signal);
dt = mean(diff(time));            
N = length(signal);               
Y = abs(fft(signal));              
f = (0:N-1)/(N*dt);            
Y_half = Y(2:floor(N/2));
f_half = f(2:floor(N/2));
[~, idx_max] = max(Y_half);
f_dom = f_half(idx_max);          
T_fft = 1 / f_dom;                
fprintf('periodo estimado fft: %.4f tsteps\n', T_fft);

%% ajuste de oscilaçao
%{
valid = isfinite(z_interface_sorted);
t_fit = t_sorted(valid);
z_fit = z_interface_sorted(valid);

z0_guess = mean(z_fit);
A_guess  = (max(z_fit) - min(z_fit)) / 2;
T_guess  = T_est;  
omega_guess = 2 * pi / T_guess;
alpha_guess = 0.001;
phi_guess = 0;

params0 = [A_guess, omega_guess, alpha_guess, phi_guess, z0_guess];

osc_model = @(p, t) ...
    p(1) * cos(p(2)*t + p(4)) .* exp(-p(3)*t) + p(5); 

opts = optimoptions('lsqcurvefit','Display','off');
[param_fit, resnorm] = lsqcurvefit(osc_model, params0, t_fit, z_fit, [], [], opts);

A_fit     = param_fit(1);
omega_fit = param_fit(2);
alpha_fit = param_fit(3);
phi_fit   = param_fit(4);
z0_fit    = param_fit(5);

T_fit = 2*pi / omega_fit;
Q_fit = omega_fit / (2 * alpha_fit);

fprintf('\n--- Ajuste de oscilação amortecida ---\n');
fprintf('A      = %.4f\n', A_fit);
fprintf('omega  = %.4f\n', omega_fit);
fprintf('T      = %.4f timesteps\n', T_fit);
fprintf('alpha  = %.6f\n', alpha_fit);
fprintf('Q      = %.2f\n', Q_fit);

% Curva ajustada
z_fit_model = osc_model(param_fit, t_fit);

% Plot com ajuste
figure;
plot(t_fit, z_fit, 'k.-', 'DisplayName', 'Sinal original'); hold on;
plot(t_fit, z_fit_model, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Ajuste oscilatório');
xlabel('Timestep');
ylabel('z_{\text{interface}} - z_c');
title('Ajuste de oscilação amortecida');
legend;
grid on;
%}