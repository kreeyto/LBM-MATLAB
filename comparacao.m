clc; clearvars; close all

nx = 10; 
ny = 10; 
nz = 10; 
fpoints = 19; 
isfluid = randi([0, 1], nx, ny, nz); 
f = rand(nx, ny, nz, fpoints); 

u_loop = zeros(nx, ny, nz);
v_loop = zeros(nx, ny, nz);
w_loop = zeros(nx, ny, nz);

for i = 1:nx
    for j = 1:ny
        for k = 1:nz
            if (isfluid(i, j, k) == 1)
                u_loop(i, j, k) = sum(f(i, j, k, [2, 16, 10, 8, 14])) - sum(f(i, j, k, [3, 11, 17, 15, 8]));
                v_loop(i, j, k) = sum(f(i, j, k, [4, 8, 15, 18, 12])) - sum(f(i, j, k, [5, 14, 9, 13, 19]));
                w_loop(i, j, k) = sum(f(i, j, k, [7, 16, 11, 18, 13])) - sum(f(i, j, k, [6, 10, 17, 12, 19]));
            end
        end
    end
end

u_vect = zeros(nx, ny, nz);
v_vect = zeros(nx, ny, nz);
w_vect = zeros(nx, ny, nz);

uindices1 = [2, 16, 10, 8, 14];
uindices2 = [3, 11, 17, 15, 8];
vindices1 = [4, 8, 15, 18, 12];
vindices2 = [5, 14, 9, 13, 19];
windices1 = [7, 16, 11, 18, 13];
windices2 = [6, 10, 17, 12, 19];

mask = (isfluid == 1);

usums = sum(f(:, :, :, uindices1), 4) - sum(f(:, :, :, uindices2), 4);
vsums = sum(f(:, :, :, vindices1), 4) - sum(f(:, :, :, vindices2), 4);
wsums = sum(f(:, :, :, windices1), 4) - sum(f(:, :, :, windices2), 4);
u_vect(mask) = usums(mask);
v_vect(mask) = vsums(mask);
w_vect(mask) = wsums(mask);

if isequal(u_loop, u_vect) && isequal(v_loop, v_vect) && isequal(w_loop, w_vect)
    disp('As duas versões produzem resultados idênticos.');
else
    disp('As duas versões produzem resultados diferentes.');
    disp('Diferenças em u:');
    disp(u_loop - u_vect);
    disp('Diferenças em v:');
    disp(v_loop - v_vect);
    disp('Diferenças em w:');
    disp(w_loop - w_vect);
end
