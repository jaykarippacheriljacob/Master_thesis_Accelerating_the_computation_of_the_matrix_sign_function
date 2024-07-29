clear;
clc;
close all;

% Define test parameters
rng(2130);
%N = 50; % Size of the matrix
% A = rand(N) + 10.0 * eye(N); % Generate a random N x N matrix
% A = gallery('poisson', 50);
%d = eig(A); % Compute the eigen values of the generated matrix

A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
N = size(A, 2); % Size of the matrix
gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12),gamma5hat);
A = Gamma5*A;

% A = A - 0.8 *speye(N);
% d = eig(full(A)); % Compute the eigen values of the generated matrix
b = randn(N, 1); % Generate a random N x 1 vector
% b = ones(N, 1);
m = 170; % No. of iterations for the krylov's subspace
trunc = 3; % Truncate orthogonalization to the last '3' vector

% plot(real(d), imag(d), '*'); % plot the real vs imaginary part of the eigen values

A_sqr = A * A;
Ab = A * b;

start = cputime;

% Call the FOM approximation function
% fom_approximation = truncated_fom(A, b, m, trunc);
fom_approximation = truncated_fom(A_sqr, Ab, m, trunc);

finish = cputime;
disp(['Time taken by truncated FOM scheme = ', num2str(finish - start), ' s']);

start = cputime;

% Compute f(A)x directly using the sign function
% exact_result = (A*(inv(sqrtm(full(A * A)))))*b;
% Save the value to exact_result.mat file
% save('exact_result.mat', 'exact_result');
% Load the value from the file
loadedData = load('exact_result.mat', 'exact_result');
exact_result = loadedData.exact_result;  % Extract the value from the structure

finish = cputime;
disp(['Time taken without FOM scheme = ', num2str(finish - start), ' s']);

% Display the relative error
rel_err = norm(exact_result - fom_approximation) / norm(exact_result);
disp(['Relative Error: ' num2str(rel_err)]);

% Set tolerance level
tol = 1e-10;
    
% Check if the residual is sufficiently small
if rel_err < tol
    disp('FOM method verified.');
else
    disp('FOM not verified.');
end
