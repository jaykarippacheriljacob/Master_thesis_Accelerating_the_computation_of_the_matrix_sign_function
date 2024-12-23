clear;
clc;
close all;

% Define test parameters
rng(2130); % setting random seed generator for reproductibility
% N = 500; % Size of the matrix
% A = rand(N) + 10.0 * eye(N); % Generate a random N x N matrix

% A = gallery('poisson', 50);

A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
N = size(A, 2); % Size of the matrix
gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12),gamma5hat);
A = Gamma5*A; 

% A = A - 0.8 *speye(N);

% d = eigs(A);
% d = eig(full(A)); % Compute the eigen values of the generated matrix
% b = randn(N, 1); % Generate a random N x 1 vector

% b = ones(N, 1);

% Initialize 2D Laplacian
% N = 100;
% e = ones(N,1);
% A = (N+1)^2*gallery('poisson',N);
% s = eigs(A,1,'SM');
% A = A/s;
% b = kron(e,e);
% b = b/norm(b);

% A = generateSparseMatrix(N);
b = randn(N, 1); % Generate a random N x 1 vector


% m = 70; % No. of iterations for the krylov's subspace
max_iter = 50; % Maximum no.of iterations for the restart of the Arnoldi decomposition

% Set tolerance level
tol = 1e-10;
% Set Error minimum decay rate for convergence
min_decay = 0.95;
% Truncate orthogonalization to the last 'trunc' vector
trunc = 5;

% plot(real(d), imag(d), '*'); % plot the real vs imaginary part of the eigen values

A_sqr = A * A;
Ab = A * b;

% start = cputime;

% Compute f(A)x directly using the sign function
% exact_result = (A*(inv(sqrtm(full(A * A)))))*b;
% exact_result = inv(sqrtm(full(A))) * b;
% Save the value to exact_result.mat file
% save('exact_result.mat', 'exact_result');
% Load the value from the file
loadedData = load('exact_result.mat', 'exact_result');
exact_result = loadedData.exact_result;  % Extract the value from the structure

% finish = cputime;
% disp(['Time taken without Quadrature based restarted arnoldi = ', num2str(finish - start), ' s']);

% Initialize the figure and set up the plot
figure;
hold on;
ylabel('log(Relative Error)');
xlabel('Iteration');
title('log(Relative Error) vs. Iteration for various m');
grid on;

% Store the colors for different curves
colors = lines(length(10:10:100)); % Generate distinct colors for each curve
legendEntries = cell(length(10:10:100), 1);
colorIndex = 1;

for m = 10:10:100
    start = cputime;
    
    % Call the Quadrature based restarted arnoldi function
    % [quadrature_approximation, iter, f] = Quad_based_restarted_arnoldi_trunc(A, b, m, max_iter, tol, min_decay, trunc);
    [quadrature_approximation, iter, f] = Quad_based_restarted_arnoldi_trunc(A_sqr, Ab, m, max_iter, tol, min_decay, trunc);
    
    finish = cputime;
    disp(['Time taken by Quadrature based restarted arnoldi with " m " - ', num2str(m), ' = ', num2str(finish - start), ' s']);
    
    % Compute relative error for each iteration
    rel_err = arrayfun(@(k) norm(exact_result - f(:, k)) / norm(exact_result), 1:iter);
    % rel_err = zeros(1, iter-1);
    % for k = 2:iter
    %     rel_err(k) = norm(f(:, k-1) - f(:, k))/norm(f(:, k-1));
    % end

    % Plot the relative error for this m
    plot(1:iter, log10(rel_err), '-o', 'Color', colors(colorIndex, :), 'MarkerIndices', 1:iter);
    
    % Update legend entries
    legendEntries{colorIndex} = ['m = ', num2str(m), ', Iterations = ', num2str(iter)];
    colorIndex = colorIndex + 1;
    
end

legend(legendEntries, 'Location', 'best');
hold off;