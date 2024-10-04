clear;
clc;
close all;

%% Define test parameters
% rng(2130);
% 
% A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
% N = size(A, 2); % Size of the matrix
% gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
% Gamma5 = kron(speye(N/12),gamma5hat);
% A = Gamma5*A;
% 
% b = randn(N, 1); % Generate a random N x 1 vector
% 
% m = 30; % No. of iterations for the krylov's subspace
% max_iter = 50; % Maximum no.of iterations for the restart of the Arnoldi decomposition
% 
% % Set tolerance level
% tol = 1e-10;
% % Set Error minimum decay rate for convergence
% min_decay = 0.95;
% % Number of target eigenvalues for implicit deflation
% thick_number = 5;
% 
% A_sqr = A * A;
% Ab = A * b;

%%
% Initialize 2D Laplacian
N = 100;
e = ones(N,1);
A = (N+1)^2*gallery('poisson',N);
s = eigs(A,1,'SM');
A = A/s;
b = kron(e,e);
b = b/norm(b);

m = 50; % No. of iterations for the krylov's subspace
max_iter = 20; % Maximum no.of iterations for the restart of the Arnoldi decomposition

% Set tolerance level
tol = 1e-14;
% Set Error minimum decay rate for convergence
min_decay = 0.95;
% Number of target eigenvalues for implicit deflation
thick_number = 5;

% A_sqr = A * A;
% Ab = A * b;

start = cputime;

% Call the Quadrature based restarted arnoldi function
[quadrature_approximation, ~, ~] = Quad_based_imp_rest_arnoldi(A, b, m, max_iter, thick_number, tol, min_decay);

finish = cputime;
disp(['Time taken by Quadrature based restarted arnoldi = ', num2str(finish - start), ' s']);

start = cputime;

% Compute f(A)x directly using the sign function
% exact_result = (A*(inv(sqrtm(full(A * A)))))*b;
% Save the value to exact_result.mat file
% save('exact_result.mat', 'exact_result');
% Load the value from the file
% loadedData = load('exact_result.mat', 'exact_result');
% exact_result = loadedData.exact_result;  % Extract the value from the structure

% Load exact solution
load exact_solutions
exact_result = exact_poisson_invsqrt;
finish = cputime;
disp(['Time taken without Quadrature based restarted arnoldi = ', num2str(finish - start), ' s']);

% Display the relative error
rel_err = norm(exact_result - quadrature_approximation) / norm(exact_result);
disp(['Relative Error: ' num2str(rel_err)]);
    
% Check if the residual is sufficiently small
if rel_err < tol
    disp('Quadrature based restarted arnoldi verified.');
else
    disp('Quadrature based restarted arnoldi not verified.');
end
