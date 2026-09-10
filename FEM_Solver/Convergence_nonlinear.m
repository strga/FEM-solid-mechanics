clear all;
close all;
clc;

% Load data
residual_data = load('/Users/strihavka/Documents/FEM_Diplom/FEM/Static_beam_diplom/residual_R_log_case5.txt');   % [iteration, ||R||]
du_data       = load('/Users/strihavka/Documents/FEM_Diplom/FEM/Static_beam_diplom/residual_du_log_case5.txt');  % [iteration, ||Δu||]

% Extract columns
iters_res = residual_data(:,1);
res_norm = residual_data(:,2);
iters_du = du_data(:,1);
du_norm = du_data(:,2);

% Plot on semilog scale
figure;
semilogy(iters_res, res_norm, 'b-o', 'DisplayName', '||R||', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
semilogy(iters_du, du_norm, 'r--s', 'DisplayName', '||\Deltau||', 'MarkerSize', 10, 'LineWidth', 2);
hold off;

xlabel('Iterace');
ylabel('conv');
%title('Newton–Raphson Convergence');
legend('Location', 'northeast', 'FontSize', 14);
grid on;
