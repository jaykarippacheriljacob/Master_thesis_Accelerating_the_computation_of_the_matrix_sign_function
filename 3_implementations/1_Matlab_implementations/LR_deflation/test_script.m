% Define matrix A
A = [4 1 2; 2 5 3; 1 1 3];

% Define vector x
x = [1; 2; 3];

% Define the number of critical eigenvalues
m = 3;

% Define the Krylov subspace dimension
k = 5;

% Compute f(A)x using LR_deflation_scheme
result_lr = lr_deflation_scheme(A, x, m, k);

% Compute f(A)x directly using the sign function
result_direct_1 = (A*((A^2)^(-1/2)))*x;

S = compute_sign_function(A);
result_direct_2 = S * x;

% Display results
disp("Approximation using LR_deflation_scheme:");
disp(result_lr);

disp("Direct computation of f(A)x (1):");
disp(result_direct_1);

disp("Direct computation of f(A)x (2):");
disp(result_direct_2);

% Compare results
if isequal(round(result_lr, 6), round(result_direct_1, 6))
    disp("Results match.");
else
    disp("Results do not match.");
end
