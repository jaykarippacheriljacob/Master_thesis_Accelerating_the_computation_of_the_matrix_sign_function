clear;
clc;
close all;

%% Define test parameters
rng(2130);

A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
N = size(A, 2); % Size of the matrix
gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12),gamma5hat);
A = Gamma5*A;

b = randn(N, 1); % Generate a random N x 1 vector

m = 19; % No. of iterations for the krylov's subspace
s = 20; % Sketch matrix row dimension
trunc = 3; % Truncate orthogonalization to the last '3' vector

A_sqr = A * A;
Ab = A * b;

start = cputime;

% Call the Sketched FOM approximation function
fom_approximation = sketched_truncated_fom(A_sqr, Ab, m, s, trunc);

finish = cputime;
disp(['Time taken by Sketched FOM scheme = ', num2str(finish - start), ' s']);

start = cputime;

% Compute f(A)x directly using the sign function
% exact_result = (A*(inv(sqrtm(full(A * A)))))*b;
% Save the value to exact_result.mat file
% save('exact_result.mat', 'exact_result');
% Load the value from the file
loadedData = load('exact_result.mat', 'exact_result');
exact_result = loadedData.exact_result;  % Extract the value from the structure

finish = cputime;
disp(['Time taken by Sketched FOM scheme with s "', num2str(s), '" & m "', num2str(m), '" and trunc "', num2str(trunc),'" = ', num2str(finish - start), ' s']);

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
