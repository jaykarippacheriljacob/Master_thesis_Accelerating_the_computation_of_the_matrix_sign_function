clear;
clc;
close all;

% check what are the possibilities of having zeros in the sign function
% calucation.

% Set the seed to a specific value, for example, 42
%rng(49);
% n = 50;
% A = rand(n) + 3.0 * eye(n);
% d = eig(A);
% plot(real(d), imag(d), '*');

A = read_matrix('4x4x4x4b6.0000id3n1.mat');
n = size(A, 2);
d = eigs(A);
plot(real(d), imag(d), '*');

% Define vector x
%x = [1; 2; 3];
x = rand(n,1);

% Define the number of critical eigenvalues
m = 5;

% Define the Krylov subspace dimension
k = 10;

start = cputime;

% Compute f(A)x using LR_deflation_scheme
result_lr = lr_deflation_scheme(A, x, m, k);

finish = cputime;
disp(['Time taken by lr deflation scheme = ', num2str(finish - start), ' s']);

start = cputime;

% Compute f(A)x directly using the sign function
result_direct_1 = (A*(inv(sqrtm(full(A * A)))))*x;

finish = cputime;
disp(['Time taken without deflation scheme = ', num2str(finish - start), ' s']);

rel_err = norm(result_direct_1 - result_lr)/norm(result_lr);
disp(['Relative error = ', num2str(rel_err)])

% Set tolerance level
tol = 1e-10;
    
% Check if the residual is sufficiently small
if rel_err < tol
    disp('LR deflation method verified.');
else
    disp('LR deflation not verified.');
end