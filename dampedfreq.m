clc; clearvars;
    
tau = 0.55;
sigma = 0.02; % which approximates to 0.02     
radius = 11.33;        

rho_one = 1.0;
rho_two = rho_one;
nu_one = 1/3 * (tau - 0.5);      
nu_two = nu_one;

monte_bool = 1;

if monte_bool == true
    m = 2;
    n = 1;
    omega_undamped = sqrt( (m*(m+1)*(m-1)*(m+1))*sigma / (radius^3*(2*m+1)) );
    alpha = ((2*n+1)^2*(nu_one*nu_two*rho_one*rho_two)^(1/2)) / (2*radius*(n*rho_two+(n+1)*rho_one)*((nu_one*rho_one)^(1/2)+(nu_two*rho_two)^(1/2)));
else
    n = 2;
    omega_undamped = sqrt(( n*(n+1)*(n-1)*(n+2)*sigma ) / ( (n*rho_one + (n+1)*rho_two) * radius^3 ));
    alpha = ( (2*n+1)^2 * sqrt(nu_one*nu_two*rho_one*rho_two) ) / ...
            ( sqrt(2) * radius * (n*rho_one + (n+1)*rho_two) * (sqrt(nu_one*rho_one) + sqrt(nu_two*rho_two)) );
end
omega_damped = omega_undamped - 0.5 * alpha * sqrt(omega_undamped) + 1/4 * alpha^2;

T_theory = 2 * pi / omega_damped;
T_LBM = 937.5;  
error_T = 100 * abs(T_theory - T_LBM) / T_theory;

fprintf('\nTable: Theoretical and Simulated Oscillation Period (n=2 mode)\n');
fprintf('---------------------------------------------------------------\n');
fprintf('%-10s %-10s %-10s %-10s %-10s\n', 'T (Theory)', 'T (LBM)', 'Error (%)', 'nu', 'rho');
fprintf('%-10.1f %-10.1f %-10.2f %-10.4f %-10.1f\n', T_theory, T_LBM, error_T, nu_one, rho_one);
