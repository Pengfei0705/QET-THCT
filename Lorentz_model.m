clear
close all
clc
load Index_data;

f = v(30:150)';
omega = 2.*pi.*v(30:150)';
eps1_exp = epsilon_real(30:150)';
eps2_exp = epsilon_imag(30:150)';

% p = [eps_inf, wp, gamma, f1, w01, Gamma1]
init_guess = [2, 2e14, 1e14, ...
              0.4, 0.6e14, 2e14, ...
              0.2, 0.9e14, 2e14, ...
              0.3, 1.48e14, 2e14];
lb = zeros(size(init_guess));
ub = Inf(size(init_guess));

% ????????? + ???
residual_fun = @(p) ... 
    [real(...
        p(1) ...
        - p(2)^2 ./ (omega.^2 + 1i * p(3) * omega) ...
        + p(4) * p(5)^2 ./ (p(5)^2 - omega.^2 - 1i * p(6) * omega) ...
        + p(7) * p(8)^2 ./ (p(8)^2 - omega.^2 - 1i * p(9) * omega) ...
        + p(10) * p(11)^2 ./ (p(11)^2 - omega.^2 - 1i * p(12) * omega) ...
    ) - eps1_exp;
    
    imag(...
        p(1) ...
        - p(2)^2 ./ (omega.^2 + 1i * p(3) * omega) ...
        + p(4) * p(5)^2 ./ (p(5)^2 - omega.^2 - 1i * p(6) * omega) ...
        + p(7) * p(8)^2 ./ (p(8)^2 - omega.^2 - 1i * p(9) * omega) ...
        + p(10) * p(11)^2 ./ (p(11)^2 - omega.^2 - 1i * p(12) * omega) ...
    ) - eps2_exp
];

% ??
options = optimoptions('lsqnonlin', 'Display', 'off', ...
    'MaxFunctionEvaluations', 1e5, 'FunctionTolerance', 1e-12);
problem = createOptimProblem('lsqnonlin', ...
    'x0', init_guess, ...
    'objective', residual_fun, ...
    'lb', lb, 'ub', ub, ...
    'options', options);
ms = MultiStart('UseParallel', false, 'Display', 'iter');
[p_fit, resnorm] = run(ms, problem, 30);

% ==============================
% ??????
% ==============================
eps_fit = p_fit(1) - p_fit(2)^2 ./ (omega.^2 + 1i * p_fit(3) * omega);
for j = 0:2
    fj = p_fit(4 + j*3);
    w0j = p_fit(5 + j*3);
    Gammaj = p_fit(6 + j*3);
    eps_fit = eps_fit + fj * w0j^2 ./ (w0j^2 - omega.^2 - 1i * Gammaj * omega);
end

% ========================================
% ???
% ========================================
figure;
subplot(2,1,1)
plot(f*1e-12, eps1_exp, 'b.', f*1e-12, real(eps_fit), 'r-');
ylabel('Re[\epsilon]');
legend('Measured', 'Fitted');
title('Drude-Lorentz (5 poles) - Real Part');

subplot(2,1,2)
plot(f*1e-12, eps2_exp, 'b.', f*1e-12, imag(eps_fit), 'r-');
ylabel('Im[\epsilon]');
legend('Measured', 'Fitted');
xlabel('Frequency (THz)');
title('Drude-Lorentz (5 poles) - Imag Part');


%% Evaluation
residual_real = real(eps_fit) - eps1_exp;
residual_imag = imag(eps_fit) - eps2_exp;

% ???? MSE
mse_real = mean(residual_real.^2);
mse_imag = mean(residual_imag.^2);

% ????? RMSE
rmse_real = sqrt(mse_real);
rmse_imag = sqrt(mse_imag);

% ???? R² (Real)
SS_res_real = sum(residual_real.^2);
SS_tot_real = sum( (eps1_exp - mean(eps1_exp)).^2 );
R2_real = 1 - SS_res_real / SS_tot_real;

% ???? R² (Imag)
SS_res_imag = sum(residual_imag.^2);
SS_tot_imag = sum( (eps2_exp - mean(eps2_exp)).^2 );
R2_imag = 1 - SS_res_imag / SS_tot_imag;

% ????
fprintf('???????\n');
fprintf('Real part:    MSE = %.4e, RMSE = %.4e, R² = %.4f\n', mse_real, rmse_real, R2_real);
fprintf('Imaginary part: MSE = %.4e, RMSE = %.4e, R² = %.4f\n', mse_imag, rmse_imag, R2_imag);



