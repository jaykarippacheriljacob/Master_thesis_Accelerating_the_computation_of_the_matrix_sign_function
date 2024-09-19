clear;
clc;
close all;

%% Select which methods to test
do_lr_deflation = true;

%% Adding paths for accessing the functions
addpath(fullfile(pwd, 'LR_deflation'));

%% Define test parameters
rng(2130); % setting random seed generator for reproducibility

A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
N = size(A, 2); % Size of the matrix
gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12),gamma5hat);
A = Gamma5*A;

b = randn(N, 1); % Generate a random N x 1 vector

m = 5; % Define the number of critical eigenvalues

% Compute f(A)x directly using the sign function
% exact_result = (A*(inv(sqrtm(full(A * A)))))*b;
% Save the value to exact_result.mat file
% save('exact_result.mat', 'exact_result');
% Load the exact result
loadedData = load('exact_result.mat', 'exact_result');
exact_result = loadedData.exact_result;  % Extract the value from the structure

% Define the range of k values
k_values = 10:10:150;

% Initialize arrays to store relative errors
rel_err_lr_deflation = zeros(length(k_values), 1);

% New way to avoid to loop over k_values
if do_lr_deflation
    start = cputime;

    % Compute f(A)x using LR_deflation_scheme
    fA_b = lr_deflation_scheme(A, b, m, k_values);

    finish = cputime;
    disp(['Time taken by lr deflation scheme = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values)
        rel_err_lr_deflation(i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end 

%% Plotting the relative errors
figure;
if do_lr_deflation
    semilogy(k_values, rel_err_lr_deflation, 'r-o', 'DisplayName', 'LR Deflation');
    hold on;
end
hold off;

xlabel('# mvms ');
ylabel('Relative Error');
title('Relative Error vs #mvms');
legend('show');
grid on;
