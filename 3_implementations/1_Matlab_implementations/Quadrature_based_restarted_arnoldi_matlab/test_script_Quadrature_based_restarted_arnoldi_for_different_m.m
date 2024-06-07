% clear;
% clc;
% close all;

% Define test parameters
% N = 50; % Size of the matrix
% A = rand(N) + 10.0 * eye(N); % Generate a random N x N matrix

% A = gallery('poisson', 50);
A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
N = size(A, 2); % Size of the matrix
% A = A - 0.9 *speye(N);

% d = eigs(A);
% d = eig(full(A)); % Compute the eigen values of the generated matrix
b = randn(N, 1); % Generate a random N x 1 vector
% b = ones(N, 1);
% m = 70; % No. of iterations for the krylov's subspace
max_iter = 1000; % Maximum no.of iterations for the restart of the Arnoldi decomposition

% plot(real(d), imag(d), '*'); % plot the real vs imaginary part of the eigen values

% A_sqr = A * A;
% Ab = A * b;

start = cputime;

% Compute f(A)x directly using the sign function
% exact_result = (A*(inv(sqrtm(full(A * A)))))*b;
exact_result = inv(sqrtm(full(A))) * b;

finish = cputime;
disp(['Time taken without Quadrature based restarted arnoldi = ', num2str(finish - start), ' s']);

% Initialize the figure and set up the plot
figure;
hold on;
ylabel('log(Relative Error)');
xlabel('m');
title('log(Relative Error) vs. m');
grid on;

for m = 450:25:900
    start = cputime;
    
    % Call the Quadrature based restarted arnoldi function
    % quadrature_approximation = Quadrature_based_restarted_arnoldi(A_sqr, Ab, m, max_iter);
    quadrature_approximation = Quadrature_based_restarted_arnoldi(A, b, m, max_iter);
    
    finish = cputime;
    disp(['Time taken by Quadrature based restarted arnoldi with " m " - ', num2str(m), ' = ', num2str(finish - start), ' s']);
    
    % Compute relative error
    rel_err = norm(exact_result - quadrature_approximation) / norm(exact_result);

    % Plot the point
    plot(m, log10(rel_err), 'bo');
    % plot(m, rel_err, 'bo');
end

hold off;