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

m = 5; % Define the number of critical eigenvalues

start = cputime;

k = 50; % No. of iterations for the krylov's subspace
k1 = 3; % No. of iterations for the krylov's subspace to be used in pre-conditioning polynomial

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

%% Calculation of f(A)b using combination of LR deflation and left polynomially preconditioned arnoldi process
fA_b = combo_lr_def_precond_polyleft(A, b, m, k, k1);

finish = cputime;
disp(['Time taken by combination of LR deflation and left polynomially preconditioned arnoldi process = ', num2str(finish - start), ' s']);

% Display the relative error
rel_err = norm(exact_result - fA_b) / norm(exact_result);
disp(['Relative Error between exact and combination of LR deflation and left polynomially preconditioned arnoldi process f(A)b: ' num2str(rel_err)]);


%% Calculation of f(A)b using combination of LR deflation and right polynomially preconditioned arnoldi process
fA_b = combo_lr_def_precond_polyright(A, b, m, k, k1);

finish = cputime;
disp(['Time taken by combination of LR deflation and right polynomially preconditioned arnoldi process = ', num2str(finish - start), ' s']);

% Display the relative error
rel_err = norm(exact_result - fA_b) / norm(exact_result);
disp(['Relative Error between exact and combination of LR deflation and right polynomially preconditioned arnoldi process f(A)b: ' num2str(rel_err)]);
