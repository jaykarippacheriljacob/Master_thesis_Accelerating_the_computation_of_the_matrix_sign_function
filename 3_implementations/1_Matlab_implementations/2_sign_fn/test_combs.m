clear;
clc;
close all;
%% Adding paths for accessing the functions
addpath("C:\Users\jkjbt\Documents\GitHub\Master_thesis_Accelerating_the_computation_of_the_matrix_sign_function\3_implementations\1_Matlab_implementations\2_sign_fn\Combo_LR_def_precond");
addpath("C:\Users\jkjbt\Documents\GitHub\Master_thesis_Accelerating_the_computation_of_the_matrix_sign_function\3_implementations\1_Matlab_implementations\2_sign_fn\Combo_LR_def_quad_sketched_FOM");
addpath("C:\Users\jkjbt\Documents\GitHub\Master_thesis_Accelerating_the_computation_of_the_matrix_sign_function\3_implementations\1_Matlab_implementations\2_sign_fn\combo_LR_def_quad_restarted_arnoldi");

%% Define test parameters
rng(2130); % setting random seed generator for reproductibility

A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
N = size(A, 2); % Size of the matrix
gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12),gamma5hat);
A = Gamma5*A;

b = randn(N, 1); % Generate a random N x 1 vector

m = 5; % Define the number of critical eigenvalues

k = 130; % No. of iterations for the krylov's subspace
s = 500; % Sketch matrix row dimension

k1 = 3; % No. of iterations for the krylov's subspace to be used in pre-conditioning polynomial

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

%% Calculation of f(A)b using combination of LR deflation and Quadrature based sketched FOM scheme process

start = cputime;

fA_b = combo_lr_quad_sketchFOM(A, b, m, k, s);

finish = cputime;
disp(['Time taken by combination of LR deflation and Quadrature based sketched FOM scheme = ', num2str(finish - start), ' s']);

% Display the relative error
rel_err = norm(exact_result - fA_b) / norm(exact_result);
disp(['Relative Error between exact and combination of LR deflation, Quad. based sketched FOM scheme f(A)b: ' num2str(rel_err)]);

%% Calculation of f(A)b using combination of LR deflation and left polynomially preconditioned arnoldi process

start = cputime;

fA_b = combo_lr_def_precond_polyleft(A, b, m, k, k1);

finish = cputime;
disp(['Time taken by combination of LR deflation and left polynomially preconditioned arnoldi process = ', num2str(finish - start), ' s']);

% Display the relative error
rel_err = norm(exact_result - fA_b) / norm(exact_result);
disp(['Relative Error between exact and combination of LR deflation and left polynomially preconditioned arnoldi process f(A)b: ' num2str(rel_err)]);

%% Calculation of f(A)b using combination of LR deflation and right polynomially preconditioned arnoldi process

start = cputime;

fA_b = combo_lr_def_precond_polyright(A, b, m, k, k1);

finish = cputime;
disp(['Time taken by combination of LR deflation and right polynomially preconditioned arnoldi process = ', num2str(finish - start), ' s']);

% Display the relative error
rel_err = norm(exact_result - fA_b) / norm(exact_result);
disp(['Relative Error between exact and combination of LR deflation and right polynomially preconditioned arnoldi process f(A)b: ' num2str(rel_err)]);

%% Calculation of f(A)b using combination of LR deflation and Quadrature based restarted arnoldi process

start = cputime;

fA_b = combo_lr_quad_restarted_arnoldi(A, b, m, k, max_iter, tol, min_decay);

finish = cputime;
disp(['Time taken by combination of LR deflation and Quadrature based restarted arnoldi = ', num2str(finish - start), ' s']);

% Display the relative error
rel_err = norm(exact_result - fA_b) / norm(exact_result);
disp(['Relative Error between exact and combination of LR deflation, Quad. based restarted arnoldi f(A)b: ' num2str(rel_err)]);

