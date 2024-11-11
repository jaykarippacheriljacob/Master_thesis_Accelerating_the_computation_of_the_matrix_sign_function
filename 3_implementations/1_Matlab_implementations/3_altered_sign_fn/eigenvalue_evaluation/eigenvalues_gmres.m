% Set matrix size and load data
A = read_matrix('../Matrix_A/periodic_L4_b3.55_k0.137n0_1.mat');
% A = read_matrix('../Matrix_A/periodic_L8_b3.55_k0.137n0_1.mat');
N = size(A, 2);

gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12),gamma5hat);
A = Gamma5*A;

m = 64; %no.of critical eigenvalues

% Start the timing
start = cputime;

% Setting up of the gmres solver
tol = 1e-6;
maxit = 200;

% Generate the Ritz values (example for m=10, adjust as needed)
m1 = 60; % Krylov subspace dimension
theta = ritz_value(A, m1); % Replace this with your actual function to compute Ritz values
precond = @(x) eval_precon_poly_inv(A, x, theta, m1);

Afun = @(x) gmres(A, x, [], tol, maxit, precond);
% Afun = @(x) gmres(A, x, [], tol, maxit);

[V, D] = eigs(Afun, N, m, 'smallestabs');

% Output results
finish = cputime;
disp('Eigenvalues:');
disp(diag(D));
% disp(diag(D));
disp(['Computation time: ', num2str(finish - start), ' seconds']);

start = cputime;
[Vm, D1] = eigs(A, m, 'smallestabs');
finish = cputime;
disp('Eigenvalues:');
disp(diag(D1));
disp(['Computation time: ', num2str(finish - start), ' seconds']);