% Test script
% Generate a test matrix
Dm = rand(1000);

% Profiling memory usage
profile -memory on
f_Dm = compute_sign_function_diag(Dm);
profile off
profview

% Profiling computation time
profile on
f_Dm = compute_sign_function_diag(Dm);
profile off
profview
