clear;
clc;
close all;

% Define test parameters
%N = 50; % Size of the matrix
% A = rand(N) + 10.0 * eye(N); % Generate a random N x N matrix
% A = gallery('poisson', 50);
%d = eig(A); % Compute the eigen values of the generated matrix
A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
N = size(A, 2); % Size of the matrix
% A = A - 0.8 *speye(N);
% d = eig(full(A)); % Compute the eigen values of the generated matrix
b = randn(N, 1); % Generate a random N x 1 vector
% b = ones(N, 1);
m = 19; % No. of iterations for the krylov's subspace
s = 20; % Sketch matrix row dimension
trunc = 3; % Truncate orthogonalization to the last '3' vector

% plot(real(d), imag(d), '*'); % plot the real vs imaginary part of the eigen values

start = cputime;

% Call the Sketched FOM approximation function
fom_approximation = sketched_truncated_fom(A, b, m, s, trunc);

finish = cputime;
disp(['Time taken by Sketched FOM scheme = ', num2str(finish - start), ' s']);

start = cputime;

% Compute f(A)x directly using the sign function
exact_result = (A*(inv(sqrtm(full(A * A)))))*b;

finish = cputime;
disp(['Time taken without Sketched FOM scheme = ', num2str(finish - start), ' s']);

% Display the relative error
rel_err = norm(exact_result - fom_approximation) / norm(exact_result);
disp(['Relative Error: ' num2str(rel_err)]);

% Set tolerance level
tol = 1e-10;
    
% Check if the residual is sufficiently small
if rel_err < tol
    disp('Sketched FOM method verified.');
else
    disp('Sketched FOM not verified.');
end
