clear;
clc;
close all;

%% Define test parameters
rng(2130); % setting random seed generator for reproductibility

A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
N = size(A, 2); % Size of the matrix
gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12),gamma5hat);
A = Gamma5*A;

b = randn(N, 1); % Generate a random N x 1 vector

m = 80; % No. of iterations for the krylov's subspace
s = 100; % Sketch matrix row dimension

A_sqr = A * A;
Ab = A * b;

start = cputime;

% Call the Sketched GMRES approximation function
Quadrature_based_sketched_FOM_approximation = Quadrature_based_sketched_FOM(A_sqr, Ab, m, s);

finish = cputime;
disp(['Time taken by Quadrature based sketched FOM scheme = ', num2str(finish - start), ' s']);

start = cputime;

% Compute f(A)x directly using the sign function
% exact_result = (A*(inv(sqrtm(full(A * A)))))*b;
% exact_result = sqrtm(inv(full(A)))*b;
% Save the value to exact_result.mat file
% save('exact_result.mat', 'exact_result');
% Load the value from the file
loadedData = load('exact_result.mat', 'exact_result');
exact_result = loadedData.exact_result;  % Extract the value from the structure

finish = cputime;
disp(['Time taken without Quadrature based sketched FOM scheme = ', num2str(finish - start), ' s']);

% Display the relative error
rel_err = norm(exact_result - Quadrature_based_sketched_FOM_approximation) / norm(exact_result);
disp(['Relative Error: ' num2str(rel_err)]);

% Set tolerance level
tol = 1e-10;
    
% Check if the residual is sufficiently small
if rel_err < tol
    disp('Quadrature based sketched FOM method verified.');
else
    disp('Quadrature based sketched FOM not verified.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start = cputime;

[x, ~, ~, ~, ~] = sketched_gmres_invsqrtm(A_sqr, Ab, m, s);
% [x, ~, ~, ~, ~] = sketched_gmres_invsqrtm(A, b, m, s);

finish = cputime;
disp(['Time taken by Quadrature based sketched FOM scheme = ', num2str(finish - start), ' s']);

% Display the relative error
rel_err1 = norm(exact_result - x) / norm(exact_result);
disp(['Relative Error: ' num2str(rel_err1)]);