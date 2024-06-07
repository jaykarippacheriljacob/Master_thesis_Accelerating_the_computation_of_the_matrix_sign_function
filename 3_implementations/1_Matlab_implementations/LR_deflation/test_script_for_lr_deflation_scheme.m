clear
clc

% check what are the possibilities of having zeros in the sign function
% calucation.

% Set the seed to a specific value, for example, 42
rng(49);
n = 50;
A = rand(n) + 10.0 * eye(n);
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

% % disp(A*((A^2)^(-1/2)))
% % Display results
% disp("Approximation using LR_deflation_scheme:");
% disp(result_lr);
% 
% %disp("Direct computation of f(A)x (1):");
% %disp(result_direct_1);
% 
% disp("Direct computation of f(A)x (2):");
% disp(result_direct_2);
% 
% % Compare results
% if isequal(round(result_lr, 6), round(result_direct_2, 6))
%     disp("Results match.");
% else
%     disp("Results do not match.");
% end
