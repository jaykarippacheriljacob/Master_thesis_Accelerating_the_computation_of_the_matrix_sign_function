clear;
clc;
close all;

% Define test parameters
rng(2130800); % setting random seed generator for reproductibility

% Read the input matrix from a file
A = read_matrix('4x4x4x4b6.0000id3n1.mat');
N = size(A, 2); % Size of the matrix
b = randn(N, 1); % Generate a random N x 1 vector

max_iter = 50; % Maximum number of iterations for the restart of the Arnoldi decomposition

% Load the exact result from the file
loadedData = load('exact_result.mat', 'exact_result');
exact_result = loadedData.exact_result;  % Extract the value from the structure

% Initialize the figure and set up the plot
figure;
hold on;
ylabel('log(Relative Error)');
xlabel('Iteration');
title('log(Relative Error) vs. Iteration for various m');
grid on;

% Store the colors for different curves
colors = lines(length(60:2:100)); % Generate distinct colors for each curve
legendEntries = cell(length(60:5:100), 1);
colorIndex = 1;

for m = 60:2:100
    start = cputime;
    
    % Call the Quadrature based restarted arnoldi function
    [quadrature_approximation, iter, f] = Quadrature_based_restarted_arnoldi(A, b, m, max_iter);
    
    finish = cputime;
    disp(['Time taken by Quadrature based restarted arnoldi with " m " - ', num2str(m), ' = ', num2str(finish - start), ' s']);
    
    % Compute relative error for each iteration
    rel_err = arrayfun(@(k) norm(exact_result - f(:, k)) / norm(exact_result), 1:iter);

    % Plot the relative error for this m
    plot(1:iter, log10(rel_err), 'Color', colors(colorIndex, :));
    
    % Update legend entries
    legendEntries{colorIndex} = ['m = ', num2str(m)];
    colorIndex = colorIndex + 1;
end

legend(legendEntries, 'Location', 'best');
hold off;
