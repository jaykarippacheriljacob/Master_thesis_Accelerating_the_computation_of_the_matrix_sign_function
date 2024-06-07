% Test script to compare Schur matrices obtained from the custom function and built-in function

% Generate a random matrix A
A = rand(5);

% Number of critical eigenvalues
m = 3;

% Compute Schur matrix using the custom function
[Sm_custom, Tm_custom, lambda_custom] = computeSchur(A, m);

% Compute Schur matrix using built-in function
[Sm_builtin, Tm_builtin] = schur(A,"complex");

% Compare Schur matrices
disp('Comparison of Schur matrices:');
disp('Custom function Schur matrix:');
disp(Sm_custom);
disp('Built-in function Schur matrix:');
disp(Sm_builtin);

% Compare upper triangular matrices
disp('Comparison of upper triangular matrices:');
disp('Custom function upper triangular matrix:');
disp(Tm_custom);
disp('Built-in function upper triangular matrix:');
disp(Tm_builtin);
