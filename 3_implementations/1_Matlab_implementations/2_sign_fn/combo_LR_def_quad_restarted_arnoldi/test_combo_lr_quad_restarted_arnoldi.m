% clear;
% clc;
% close all;
%% Define test parameters
rng(2130); % setting random seed generator for reproductibility

A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
N = size(A, 2); % Size of the matrix
gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12),gamma5hat);
A = Gamma5*A;

b = randn(N, 1); % Generate a random N x 1 vector

m = 4; % Define the number of critical eigenvalues

k = 15; % No. of iterations for the krylov's subspace
max_iter = 50; % Maximum no.of iterations for the restart of the Arnoldi decomposition

% Set tolerance level
tol = 1e-10;
% Set Error minimum decay rate for convergence
min_decay = 0.95;

%% Verification of the result with exact calculation of inverse square root

start = cputime;

% exact_result = (A*(inv(sqrtm(full(A * A)))))*b;
% Save the value to exact_result.mat file
% save('exact_result.mat', 'exact_result');
% Load the value from the file
loadedData = load('exact_result.mat', 'exact_result');
exact_result = loadedData.exact_result;  % Extract the value from the structure

finish = cputime;
disp(['Time taken for direct calculation of f(A)b = ', num2str(finish - start), ' s']);

%% Calculation of f(A)b using combination of LR deflation and Quadrature based restarted arnoldi process

start = cputime;

[iter, fA_b] = combo_lr_quad_restarted_arnoldi(A, b, m, k, max_iter, tol, min_decay);
cost_quad_restarted_arnoldi = 18 + 4 * k + m^2 + 3 * iter + 4 * m * iter + m^2 * iter;
display(cost_quad_restarted_arnoldi);

finish = cputime;
disp(['Time taken by combination of LR deflation and Quadrature based restarted arnoldi = ', num2str(finish - start), ' s']);

% Display the relative error
rel_err = norm(exact_result - fA_b) / norm(exact_result);
disp(['Relative Error between exact and combination of LR deflation, Quad. based restarted arnoldi f(A)b: ' num2str(rel_err)]);
