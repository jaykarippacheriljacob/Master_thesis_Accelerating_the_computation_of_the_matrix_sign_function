clear;
clc;
close all;

%% Adding paths for accessing the functions
addpath(fullfile(pwd, 'Combo_LR_def_precond'));
addpath(fullfile(pwd, 'Combo_LR_def_quad_sketched_FOM'));
addpath(fullfile(pwd, 'Combo_LR_def_quad_trunc_sketched_FOM'));
addpath(fullfile(pwd, 'Combo_LR_def_quad_restarted_arnoldi'));

%% Define test parameters
rng(2130); % setting random seed generator for reproducibility

A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
N = size(A, 2); % Size of the matrix
gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12),gamma5hat);
A = Gamma5*A;

b = randn(N, 1); % Generate a random N x 1 vector

m = 5; % Define the number of critical eigenvalues

s = 500; % Sketch matrix row dimension

k1 = 3; % No. of iterations for the Krylov's subspace to be used in pre-conditioning polynomial

trunc = 10; % Truncate orthogonalization to the last 'trunc' vector

max_iter = 50; % Maximum no. of iterations for the restart of the Arnoldi decomposition
% Set tolerance level
tol = 1e-10;
% Set Error minimum decay rate for convergence
min_decay = 0.95;

% Load the exact result
loadedData = load('exact_result.mat', 'exact_result');
exact_result = loadedData.exact_result;  % Extract the value from the structure

% Define the range of k values
k_values = 10:10:150;

% Initialize arrays to store relative errors
rel_err_quad_sketchFOM = zeros(length(k_values), 1);
rel_err_quad_truncsketchFOM = zeros(length(k_values), 1);
rel_err_precond_polyleft = zeros(length(k_values), 1);
rel_err_precond_polyright = zeros(length(k_values), 1);
rel_err_quad_restarted_arnoldi = zeros(length(k_values), 1);

% Initializing array to store the cost
cost_quad_restarted_arnoldi = zeros(length(k_values), 1);
cost_quad_sketchFOM = zeros(length(k_values), 1);
cost_precond_polyleft = zeros(length(k_values), 1);
cost_precond_polyright = zeros(length(k_values), 1);
cost_quad_truncsketchFOM = zeros(length(k_values), 1);

% Loop over the range of k values
for i = 1:length(k_values)
    k = k_values(i);

    %% Calculation of f(A)b using combination of LR deflation and Quadrature based restarted Arnoldi process
    [iter, fA_b] = combo_lr_quad_restarted_arnoldi(A, b, m, k, max_iter, tol, min_decay);
    rel_err_quad_restarted_arnoldi(i) = norm(exact_result - fA_b) / norm(exact_result);
    cost_quad_restarted_arnoldi(i) = 1 + iter*2*k;

    %% Calculation of f(A)b using combination of LR deflation and Quadrature based sketched FOM scheme process
    fA_b = combo_lr_quad_sketchFOM(A, b, m, k, s);
    rel_err_quad_sketchFOM(i) = norm(exact_result - fA_b) / norm(exact_result);
    cost_quad_sketchFOM(i) = 1 + 2*k;

    %% Calculation of f(A)b using combination of LR deflation and Quadrature based truncated sketched FOM scheme process
    fA_b = combo_lr_quad_truncsketchFOM(A, b, m, k, s, trunc);
    rel_err_quad_truncsketchFOM(i) = norm(exact_result - fA_b) / norm(exact_result);
    cost_quad_truncsketchFOM(i) = 1 + 2*k;

    %% Calculation of f(A)b using combination of LR deflation and left polynomially preconditioned Arnoldi process
    fA_b = combo_lr_def_precond_polyleft(A, b, m, k, k1);
    rel_err_precond_polyleft(i) = norm(exact_result - fA_b) / norm(exact_result);
    cost_precond_polyleft(i) = 2*k1 + 1 + 2*k1 + k * (2 + 2*k1 + 2*k1);

    %% Calculation of f(A)b using combination of LR deflation and right polynomially preconditioned Arnoldi process
    fA_b = combo_lr_def_precond_polyright(A, b, m, k, k1);
    rel_err_precond_polyright(i) = norm(exact_result - fA_b) / norm(exact_result);
    cost_precond_polyright(i) = 2*k1 + 1 + 2*k1 + k * (2 + 2*k1 + 2*k1);
end

%% Plotting the relative errors
figure;
semilogy(cost_quad_sketchFOM, rel_err_quad_sketchFOM, 'r o', 'DisplayName', 'Quad. based sketched FOM');
% for i = 1:length(k_values)
%     text(cost_quad_sketchFOM(i), rel_err_quad_sketchFOM(i), sprintf('%.1f', k_values(i)), 'Color', 'red', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
% end

hold on;
semilogy(cost_quad_truncsketchFOM, rel_err_quad_truncsketchFOM, 'y ^', 'DisplayName', 'Quad. based trun.sketched FOM');

semilogy(cost_precond_polyleft, rel_err_precond_polyleft, 'g *', 'DisplayName', 'Left precond. Arnoldi');
% for i = 1:length(k_values)
%     text(cost_precond_polyleft(i), rel_err_precond_polyleft(i), sprintf('%.1f', k_values(i)), 'Color', 'red', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
% end

semilogy(cost_precond_polyright, rel_err_precond_polyright, 'b ^', 'DisplayName', 'Right precond. Arnoldi');
% for i = 1:length(k_values)
%     text(cost_precond_polyright(i), rel_err_precond_polyright(i), sprintf('%.1f', k_values(i)), 'Color', 'red', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
% end

semilogy(cost_quad_restarted_arnoldi, rel_err_quad_restarted_arnoldi, 'k s', 'DisplayName', 'Quad. based restarted Arnoldi');
% for i = 1:length(k_values)
%     text(cost_quad_restarted_arnoldi(i), rel_err_quad_restarted_arnoldi(i), sprintf('%.1f', k_values(i)), 'Color', 'red', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
% end
hold off;

xlabel('Cost of computation {A*b}');
ylabel('Relative Error');
title('Relative Error vs Cost of computation');
legend('show');
grid on;
