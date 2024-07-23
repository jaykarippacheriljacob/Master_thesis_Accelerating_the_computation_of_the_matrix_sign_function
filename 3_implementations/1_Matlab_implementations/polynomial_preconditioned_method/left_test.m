clear;
clc;
close all;
%% Define test parameters
rng(2130); % setting random seed generator for reproductibility
% N = 50; % Size of the matrix
% A = rand(N) + 10.0 * eye(N); % Generate a random N x N matrix
% A = gallery('poisson', 50);
% d = eig(A); % Compute the eigen values of the generated matrix
A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
N = size(A, 2); % Size of the matrix
gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12),gamma5hat);
A = Gamma5*A;
% A = A - 0.8 *speye(N);
% d = eig(full(A)); % Compute the eigen values of the generated matrix
b = randn(N, 1); % Generate a random N x 1 vector
% b = ones(N, 1);
m = 25; % No. of iterations for the krylov's subspace

% plot(real(d), imag(d), '*'); % plot the real vs imaginary part of the eigen values

A_sqr = A * A;
Ab = A * b;

start = cputime;
%% Finding the Preconditioned polynomial using Newton's interpolation using divided differences
m1 = 5;
% p_A = pre_condi_poly(A, b, m1);
p_A = pre_condi_poly(A_sqr, Ab, m1);

%% Calculation of left polynomially preconditioned arnoldi process for A^{-1/2}
% fA_b = left_precondi_poly_fom(p_A, A, b, m);
fA_b = left_precondi_poly_fom(p_A, A_sqr, Ab, m);

finish = cputime;
disp(['Time taken by left polynomially preconditioned arnoldi process = ', num2str(finish - start), ' s']);
%% Verification of the result with exact calculation of inverse square root
start = cputime;

% exact_result = sqrtm(inv(full(A)))*b;
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
