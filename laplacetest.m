clc; clearvars; close all

% Caminhos dos arquivos
base_path = fullfile(getenv('HOME'),'Desktop','MULTIC-TS-LBM','bin','D3Q19','6969');
rho_file  = fullfile(base_path, '6969_rho005000.bin');  
phi_file  = fullfile(base_path, '6969_phi005000.bin');  

[NX,NY,NZ] = deal(64);
[xc,yc,zc] = deal(round((NX-1)/2), round((NY-1)/2), round((NZ-1)/2));

fid = fopen(rho_file,'rb');
rho = fread(fid, 'float32=>float32');
fclose(fid);
rho = reshape(rho, [NX,NY,NZ]);
rho = permute(rho, [3,2,1]);  

fid = fopen(phi_file,'rb');
phi = fread(fid, 'float32=>float32');
fclose(fid);
phi = reshape(phi, [NX,NY,NZ]);
phi = permute(phi, [3,2,1]); 

in_mask  = phi > 0.9;    
out_mask = phi < 0.1;   

cs2 = 1/3;
p_field = cs2 * rho;

p_in  = mean(p_field(in_mask));
p_out = mean(p_field(out_mask));
delta_p = p_in - p_out;

R = 9; 

sigma_est = (delta_p * R) / 2;

fprintf('--- Teste de Laplace com campo de fase (phi) ---\n');
fprintf('p_in     = %.6f\n', p_in);
fprintf('p_out    = %.6f\n', p_out);
fprintf('Δp       = %.6e\n', delta_p);
fprintf('Raio     = %.2f\n', R);
fprintf('Sigma ≈  %.6e (unidades LBM)\n', sigma_est);
