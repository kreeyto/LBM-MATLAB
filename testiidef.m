clc; clearvars; close all
 
% Configuração inicial
[nx, ny, nz] = deal(64);
fpoints = 19;
cssq = 1/3;
w = rand(1, fpoints);
cix = randi([-1, 1], 1, fpoints);
ciy = randi([-1, 1], 1, fpoints);
ciz = randi([-1, 1], 1, fpoints);
f = rand(nx, ny, nz, fpoints);
rho = rand(nx, ny, nz);
rho_loop = rho;
rho_vec = rho;
[ffx, ffy, ffz] = deal(rand(nx, ny, nz));
fneq_loop = zeros(fpoints,1,1); 
fneq_vec = fneq_loop;

% Inicialização de variáveis
[ux_loop, uy_loop, uz_loop] = deal(zeros(nx, ny, nz));
[ux_vec, uy_vec, uz_vec] = deal(zeros(nx, ny, nz));
[pxx_loop, pyy_loop, pzz_loop] = deal(zeros(nx, ny, nz));
[pxx_vec, pyy_vec, pzz_vec] = deal(zeros(nx, ny, nz));
[pxy_loop, pxz_loop, pyz_loop] = deal(zeros(nx, ny, nz));
[pxy_vec, pxz_vec, pyz_vec] = deal(zeros(nx, ny, nz));

% Abordagem com loops
for i = 2:nx-1
    for j = 2:ny-1
        for k = 2:nz-1
            ux_loop(i, j, k) = sum(f(i, j, k, [2, 16, 10, 8, 14])) - sum(f(i, j, k, [3, 11, 17, 15, 8]));
            uy_loop(i, j, k) = sum(f(i, j, k, [4, 8, 15, 18, 12])) - sum(f(i, j, k, [5, 14, 9, 13, 19]));
            uz_loop(i, j, k) = sum(f(i, j, k, [7, 16, 11, 18, 13])) - sum(f(i, j, k, [6, 10, 17, 12, 19]));
            ux_loop(i, j, k) = ux_loop(i, j, k) ./ rho_loop(i, j, k) + ffx(i, j, k) * 0.5 ./ rho_loop(i, j, k);
            uy_loop(i, j, k) = uy_loop(i, j, k) ./ rho_loop(i, j, k) + ffy(i, j, k) * 0.5 ./ rho_loop(i, j, k);
            uz_loop(i, j, k) = uz_loop(i, j, k) ./ rho_loop(i, j, k) + ffz(i, j, k) * 0.5 ./ rho_loop(i, j, k);
            uu_loop = 0.5 * (ux_loop(i, j, k).^2 + uy_loop(i, j, k).^2 + uz_loop(i, j, k).^2) / cssq;
            rho_loop(i, j, k) = sum(f(i, j, k, :), 4);
            for l = 1:fpoints
                udotc_loop = (ux_loop(i, j, k) * cix(l) + uy_loop(i, j, k) * ciy(l) + uz_loop(i, j, k) * ciz(l)) / cssq;
                HeF_loop = 0.5 .* (w(l) * (rho_loop(i, j, k) + rho_loop(i, j, k) .* (udotc_loop + 0.5 .* udotc_loop.^2 - uu_loop))) .* ...
                    ((cix(l) - ux_loop(i, j, k)) .* ffx(i, j, k) + ...
                     (ciy(l) - uy_loop(i, j, k)) .* ffy(i, j, k) + ...
                     (ciz(l) - uz_loop(i, j, k)) .* ffz(i, j, k)) ./ (rho_loop(i, j, k) .* cssq);
                feq_loop = w(l) * (rho(i, j, k) + rho_loop(i, j, k) .* (udotc_loop + 0.5 .* udotc_loop.^2 - uu_loop)) - HeF_loop;
                fneq_loop(l) = f(i, j, k, l) - feq_loop;
            end
            pxx_loop(i, j, k) = sum(fneq_loop([2, 3, 8, 9, 10, 11, 14, 15, 16, 17]));
            pyy_loop(i, j, k) = sum(fneq_loop([4, 5, 8, 9, 12, 13, 14, 15, 18, 19]));
            pzz_loop(i, j, k) = sum(fneq_loop([6, 7, 10, 11, 12, 13, 16, 17, 18, 19]));
            pxy_loop(i, j, k) = fneq_loop(8) + fneq_loop(9) - fneq_loop(14) - fneq_loop(15);
            pxz_loop(i, j, k) = fneq_loop(10) + fneq_loop(11) - fneq_loop(16) - fneq_loop(17);
            pyz_loop(i, j, k) = fneq_loop(12) + fneq_loop(13) - fneq_loop(18) - fneq_loop(19);
        end
    end
end

% Abordagem vetorizada
ux_vec(2:nx-1, 2:ny-1, 2:nz-1) = sum(f(2:nx-1, 2:ny-1, 2:nz-1, [2, 16, 10, 8, 14]), 4) - sum(f(2:nx-1, 2:ny-1, 2:nz-1, [3, 11, 17, 15, 8]), 4);
uy_vec(2:nx-1, 2:ny-1, 2:nz-1) = sum(f(2:nx-1, 2:ny-1, 2:nz-1, [4, 8, 15, 18, 12]), 4) - sum(f(2:nx-1, 2:ny-1, 2:nz-1, [5, 14, 9, 13, 19]), 4);
uz_vec(2:nx-1, 2:ny-1, 2:nz-1) = sum(f(2:nx-1, 2:ny-1, 2:nz-1, [7, 16, 11, 18, 13]), 4) - sum(f(2:nx-1, 2:ny-1, 2:nz-1, [6, 10, 17, 12, 19]), 4);
ux_vec(2:nx-1, 2:ny-1, 2:nz-1) = ux_vec(2:nx-1, 2:ny-1, 2:nz-1) ./ rho_vec(2:nx-1, 2:ny-1, 2:nz-1) + ffx(2:nx-1, 2:ny-1, 2:nz-1) * 0.5 ./ rho_vec(2:nx-1, 2:ny-1, 2:nz-1);
uy_vec(2:nx-1, 2:ny-1, 2:nz-1) = uy_vec(2:nx-1, 2:ny-1, 2:nz-1) ./ rho_vec(2:nx-1, 2:ny-1, 2:nz-1) + ffy(2:nx-1, 2:ny-1, 2:nz-1) * 0.5 ./ rho_vec(2:nx-1, 2:ny-1, 2:nz-1);
uz_vec(2:nx-1, 2:ny-1, 2:nz-1) = uz_vec(2:nx-1, 2:ny-1, 2:nz-1) ./ rho_vec(2:nx-1, 2:ny-1, 2:nz-1) + ffz(2:nx-1, 2:ny-1, 2:nz-1) * 0.5 ./ rho_vec(2:nx-1, 2:ny-1, 2:nz-1);

uu_vec = 0.5 * (ux_vec(2:nx-1,2:ny-1,2:nz-1).^2 + uy_vec(2:nx-1,2:ny-1,2:nz-1).^2 + uz_vec(2:nx-1,2:ny-1,2:nz-1).^2) / cssq;

rho_vec(2:nx-1,2:ny-1,2:nz-1) = sum(f(2:nx-1,2:ny-1,2:nz-1,:),4);
for l = 1:fpoints

    udotc_vec = (ux_vec(2:nx-1,2:ny-1,2:nz-1) * cix(l) + uy_vec(2:nx-1,2:ny-1,2:nz-1) * ciy(l) + uz_vec(2:nx-1,2:ny-1,2:nz-1) * ciz(l)) / cssq;

    HeF_vec = 0.5 .* (w(l) * (rho_vec(2:nx-1,2:ny-1,2:nz-1) + rho_vec(2:nx-1,2:ny-1,2:nz-1) .* (udotc_vec + 0.5.*udotc_vec.^2 - uu_vec))) ...
        .* ((cix(l) - ux_vec(2:nx-1,2:ny-1,2:nz-1)) .* ffx(2:nx-1,2:ny-1,2:nz-1) + ...
            (ciy(l) - uy_vec(2:nx-1,2:ny-1,2:nz-1)) .* ffy(2:nx-1,2:ny-1,2:nz-1) + ...
            (ciz(l) - uz_vec(2:nx-1,2:ny-1,2:nz-1)) .* ffz(2:nx-1,2:ny-1,2:nz-1) ...
           ) ./ (rho_vec(2:nx-1,2:ny-1,2:nz-1) .* cssq);

    feq_vec = w(l) * (rho_vec(2:nx-1,2:ny-1,2:nz-1) + rho_vec(2:nx-1,2:ny-1,2:nz-1) .* (udotc_vec + 0.5.*udotc_vec.^2 - uu_vec)) - HeF_vec;

    fneq_vec(l) = f(2:nx-1, 2:ny-1, 2:nz-1, l) - feq_vec;

end

pxx_vec(2:nx-1,2:ny-1,2:nz-1) = sum(fneq_vec(2:nx-1,2:ny-1,2:nz-1,[2,3,8,9,10,11,14,15,16,17]),4);
pyy_vec(2:nx-1,2:ny-1,2:nz-1) = sum(fneq_vec(2:nx-1,2:ny-1,2:nz-1,[4,5,8,9,12,13,14,15,18,19]),4);
pzz_vec(2:nx-1,2:ny-1,2:nz-1) = sum(fneq_vec(2:nx-1,2:ny-1,2:nz-1,[6,7,10,11,12,13,16,17,18,19]),4);
pxy_vec(2:nx-1,2:ny-1,2:nz-1) = fneq_vec(8) + fneq_vec(9) - fneq_vec(14) - fneq_vec(15);
pxz_vec(2:nx-1,2:ny-1,2:nz-1) = fneq_vec(10) + fneq_vec(11) - fneq_vec(16) - fneq_vec(17);
pyz_vec(2:nx-1,2:ny-1,2:nz-1) = fneq_vec(12) + fneq_vec(13) - fneq_vec(18) - fneq_vec(19);

% Comparação dos resultados
% IGUAL is_ux_identical = isequal(ux_loop, ux_vec);
% IGUAL is_uy_identical = isequal(uy_loop, uy_vec);
% IGUAL is_uz_identical = isequal(uz_loop, uz_vec);
% NEM FAZ SENTIDO is_uu_identical = isequal(uu_loop, uu_vec);
% IGUAL is_rho_identical = isequal(rho_loop, rho_vec);
% NEM FAZ SENTIDO is_udotc_identical = isequal(udotc_loop, udotc_vec);
% NEM FAZ SENTIDO is_HeF_identical = isequal(HeF_loop, HeF_vec);
%is_feq_identical = isequal(feq_loop, feq_vec);
%is_fneq_identical = isequal(fneq_loop, fneq_vec);
is_pxx_identical = isequal(pxx_loop, pxx_vec);
is_pyy_identical = isequal(pyy_loop, pyy_vec);
is_pzz_identical = isequal(pzz_loop, pzz_vec);
is_pxy_identical = isequal(pxy_loop, pxy_vec);
is_pxz_identical = isequal(pxz_loop, pxz_vec);
is_pyz_identical = isequal(pyz_loop, pyz_vec);

if is_pxx_identical && is_pyy_identical ...
    && is_pzz_identical && is_pxy_identical && is_pxz_identical && is_pyz_identical
    disp('Os resultados são idênticos!');
else
    disp('Os resultados são diferentes.');
    if ~is_pxx_identical, disp('Diferenças encontradas em "pxx".'); end
    if ~is_pyy_identical, disp('Diferenças encontradas em "pyy".'); end
    if ~is_pzz_identical, disp('Diferenças encontradas em "pzz".'); end
    if ~is_pxy_identical, disp('Diferenças encontradas em "pxy".'); end
    if ~is_pxz_identical, disp('Diferenças encontradas em "pxz".'); end
    if ~is_pyz_identical, disp('Diferenças encontradas em "pyz".'); end
end

