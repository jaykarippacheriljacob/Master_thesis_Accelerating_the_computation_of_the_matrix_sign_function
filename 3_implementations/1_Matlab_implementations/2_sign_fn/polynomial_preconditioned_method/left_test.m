clear;
clc;
close all;
%% Define test parameters
rng(2130); % setting random seed generator for reproductibility

A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
N = size(A, 2); % Size of the matrix
gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12),gamma5hat);
A = Gamma5*A;

b = randn(N, 1); % Generate a random N x 1 vector

m = 23; % No. of iterations for the krylov's subspace

% A_sqr = A * A;
% Ab = A * b;

start = cputime;
%% Finding the Preconditioned polynomial using Newton's interpolation using divided differences
m1 = 4;
% p_A = pre_condi_poly(A, b, m1);
% p_A = pre_condi_poly(A_sqr, m1);
% theta = ritz_value(A_sqr, m1);
% p_A = eval_pre_condi_poly(A_sqr, theta, m1);

%% Calculation of left polynomially preconditioned arnoldi process for A^{-1/2}
% fA_b = left_precondi_poly_fom(p_A, A, b, m);
% fA_b = left_precondi_poly_fom(p_A, A_sqr, Ab, m);
fA_b = left_precondi_poly_fom(A, b, m, m1);

finish = cputime;
disp(['Time taken by left polynomially preconditioned arnoldi process = ', num2str(finish - start), ' s']);
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

%% Display the relative error
rel_err = norm(exact_result - fA_b) / norm(exact_result);
disp(['Relative Error between exact and left polynomially preconditioned arnoldi f(A)b: ' num2str(rel_err)]);
