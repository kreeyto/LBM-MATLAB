clc; clearvars;
% ------------------------------
% Droplet Oscillation Period - Theoretical vs Measured (LBM)
% ------------------------------

% Parâmetros ajustáveis:
sigma = 0.1;       % Surface tension
rho = 1.0;          % Density (assumindo rho1 = rho2)
nu = 0.1667;          % Kinematic viscosity
R0 = 35;          % Droplet equilibrium radius (lu)

% Cálculo da frequência natural (Lamb, sem damping):
omega2_star_sq = (24 * sigma) / (5 * rho * R0^3);
omega2_star = sqrt(omega2_star_sq);

% Cálculo de alpha (Miller-Scriven):
alpha = (5 * sqrt(nu)) / (2 * rho * R0);

% Frequência corrigida (Miller-Scriven com damping):
omega2_damped = omega2_star - 0.5 * alpha * sqrt(omega2_star) + 0.25 * alpha^2;

% Período teórico (em timesteps LBM, assumindo dt = 1):
T_theory = 2 * pi / omega2_damped;

% Exemplo de valor "medido" vindo de um LBM:
T_LBM = 9500;  % <-- Substitua por seu valor real de simulação

% Erro percentual:
error_T = 100 * abs(T_LBM - T_theory) / T_theory;

% Exibir resultados:
fprintf('\nTable: Theoretical and Simulated Oscillation Period (n=2 mode)\n');
fprintf('---------------------------------------------------------------\n');
fprintf('%-10s %-10s %-10s %-10s %-10s\n', 'T (Theory)', 'T (LBM)', 'Error (%)', 'nu', 'rho');
fprintf('%-10.1f %-10.1f %-10.2f %-10.4f %-10.1f\n', T_theory, T_LBM, error_T, nu, rho);
