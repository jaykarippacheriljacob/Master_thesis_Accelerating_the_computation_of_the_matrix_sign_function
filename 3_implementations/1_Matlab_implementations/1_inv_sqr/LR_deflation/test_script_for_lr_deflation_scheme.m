clear;
clc;
close all;

% check what are the possibilities of having zeros in the sign function
% calucation.

% Define test parameters

rng(2130); % setting random seed generator for reproductibility
% n = 20; % Size of the matrix
% A = rand(n) + 10.0 * eye(n); % Generate a random N x N matrix
% d = eig(A); % Compute the eigen values of the generated matrix
A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
n = size(A, 2); % Size of the matrix
gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(n/12),gamma5hat);
A = Gamma5*A;
% d = eigs(A); % Compute the eigen values of the generated matrix
% plot(real(d), imag(d), '*'); % plot the real vs imaginary part of the eigen values

% Define vector x
%x = [1; 2; 3];
x = rand(n,1);

% Define the number of critical eigenvalues
m = 10;

% Define the Krylov subspace dimension
k = 100;

% Verify compute_eigenvectors function of a matrix.
verify_compute_eigenvectors(A, m, 1);

start = cputime;

% Compute f(A)x using LR_deflation_scheme
result_lr = lr_deflation_scheme(A, x, m, k);

finish = cputime;
disp(['Time taken by lr deflation scheme = ', num2str(finish - start), ' s']);

start = cputime;

% Compute f(A)x directly using the sign function
% exact_result = (A*(inv(sqrtm(full(A * A)))))*x;
% exact_result = sqrtm(inv(full(A)))*x;
% Save the value to exact_result.mat file
% save('exact_result.mat', 'exact_result');
% Load the value from the file
loadedData = load('exact_result.mat', 'exact_result');
exact_result = loadedData.exact_result;  % Extract the value from the structure

finish = cputime;
disp(['Time taken without deflation scheme = ', num2str(finish - start), ' s']);

rel_err = norm(exact_result - result_lr)/norm(result_lr);
disp(['Relative error = ', num2str(rel_err)])

% Set tolerance level
tol = 1e-10;
    
% Check if the residual is sufficiently small
if rel_err < tol
    disp('LR deflation method verified.');
else
    disp('LR deflation not verified.');
end