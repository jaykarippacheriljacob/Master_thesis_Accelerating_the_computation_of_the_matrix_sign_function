clear;
clc;
close all;

% Define test parameters
% N = 50; % Size of the matrix
% A = rand(N) + 10.0 * eye(N); % Generate a random N x N matrix

% A = gallery('poisson', 50);
A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
N = size(A, 2); % Size of the matrix
% A = A - 0.9 *speye(N);

% d = eigs(A);
% d = eig(full(A)); % Compute the eigen values of the generated matrix
b = randn(N, 1); % Generate a random N x 1 vector
% b = ones(N, 1);
m = 200; % No. of iterations for the krylov's subspace
max_iter = 1000; % Maximum no.of iterations for the restart of the Arnoldi decomposition

% plot(real(d), imag(d), '*'); % plot the real vs imaginary part of the eigen values

% A_sqr = A * A;
% Ab = A * b;

start = cputime;

% Call the Quadrature based restarted arnoldi function
quadrature_approximation = Quadrature_based_restarted_arnoldi(A, b, m, max_iter);

finish = cputime;
disp(['Time taken by Quadrature based restarted arnoldi = ', num2str(finish - start), ' s']);

start = cputime;

% Compute f(A)x directly using the sign function
%exact_result = (A*(inv(sqrtm(full(A * A)))))*b;
exact_result = inv(sqrtm(full(A))) * b;

finish = cputime;
disp(['Time taken without Quadrature based restarted arnoldi = ', num2str(finish - start), ' s']);

% Display the relative error
rel_err = norm(exact_result - quadrature_approximation) / norm(exact_result);
disp(['Relative Error: ' num2str(rel_err)]);

% Set tolerance level
tol = 1e-10;
    
% Check if the residual is sufficiently small
if rel_err < tol
    disp('Quadrature based restarted arnoldi verified.');
else
    disp('Quadrature based restarted arnoldi not verified.');
end
