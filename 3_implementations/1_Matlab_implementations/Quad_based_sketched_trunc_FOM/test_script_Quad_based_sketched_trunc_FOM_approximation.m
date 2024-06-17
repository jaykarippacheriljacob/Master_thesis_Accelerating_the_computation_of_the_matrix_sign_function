clear;
clc;
close all;

% Define test parameters
% N = 50; % Size of the matrix
% A = rand(N) + 10.0 * eye(N); % Generate a random N x N matrix
% A = gallery('poisson', 50);
% d = eig(A); % Compute the eigen values of the generated matrix
A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
N = size(A, 2); % Size of the matrix
% A = A - 0.8 *speye(N);
% d = eig(full(A)); % Compute the eigen values of the generated matrix
b = randn(N, 1); % Generate a random N x 1 vector
% b = ones(N, 1);
m = 100; % No. of iterations for the krylov's subspace
s = 101; % Sketch matrix row dimension

% plot(real(d), imag(d), '*'); % plot the real vs imaginary part of the eigen values

A_sqr = A * A;
Ab = A * b;

start = cputime;

% Call the Sketched GMRES approximation function
Quadrature_based_sketched_FOM_approximation = Quadrature_based_sketched_FOM(A_sqr, Ab, m, s);
% gmres_approximation = sketched_GMRES(A, b, m, s);

finish = cputime;
disp(['Time taken by Quadrature based sketched FOM scheme = ', num2str(finish - start), ' s']);

start = cputime;

% Compute f(A)x directly using the sign function
exact_result = (A*(inv(sqrtm(full(A * A)))))*b;
% exact_result = sqrtm(inv(full(A)))*b;

finish = cputime;
disp(['Time taken without Quadrature based sketched FOM scheme = ', num2str(finish - start), ' s']);

% Display the relative error
rel_err = norm(exact_result - Quadrature_based_sketched_FOM_approximation) / norm(exact_result);
disp(['Relative Error: ' num2str(rel_err)]);

% Set tolerance level
tol = 1e-10;
    
% Check if the residual is sufficiently small
if rel_err < tol
    disp('Quadrature based sketched FOM method verified.');
else
    disp('Quadrature based sketched FOM not verified.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start = cputime;

[x, ~, ~, ~, ~] = sketched_gmres_invsqrtm(A_sqr, Ab, m, s);
% [x, ~, ~, ~, ~] = sketched_gmres_invsqrtm(A, b, m, s);

finish = cputime;
disp(['Time taken by Quadrature based sketched FOM scheme = ', num2str(finish - start), ' s']);

% Display the relative error
rel_err1 = norm(exact_result - x) / norm(exact_result);
disp(['Relative Error: ' num2str(rel_err1)]);