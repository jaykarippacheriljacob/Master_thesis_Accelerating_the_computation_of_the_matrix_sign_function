clear;
clc;

% check what are the possibilities of having zeros in the sign function
% calucation.

% Set the seed to a specific value, for example, 42
%rng(49);
n = 50;
A = rand(n) + 3.0 * eye(n);
d = eig(A);
plot(real(d), imag(d), '*');

% Define vector x
%x = [1; 2; 3];
x = rand(n,1);

% Define the number of critical eigenvalues
m = 5;

% Define the Krylov subspace dimension
k = 10;

% Compute f(A)x using LR_deflation_scheme
result_lr = lr_deflation_scheme(A, x, m, k);

% Compute f(A)x directly using the sign function
%result_direct_1 = (A*((A^2)^(-1/2)))*x;
result_direct_1 = (A*(inv(sqrtm(A * A))))*x;

%S = compute_sign_function(A);
% disp(S);
%result_direct_2 = S * x;

%trial = sign(A);
% disp(trial);

rel_err = norm(result_direct_1 - result_lr)/norm(result_lr);
disp(rel_err)

% Set tolerance level
tol = 1e-10;
    
% Check if the residual is sufficiently small
if rel_err < tol
    disp('LR deflation method verified.');
else
    disp('LR deflation not verified.');
end