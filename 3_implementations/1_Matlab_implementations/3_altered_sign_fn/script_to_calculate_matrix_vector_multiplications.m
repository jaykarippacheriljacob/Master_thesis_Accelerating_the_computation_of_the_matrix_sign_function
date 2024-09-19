clear;
clc;
close all;

%% Select which methods to test
do_lr_deflation = true;
do_left_precondi_poly_fom = true;
do_right_precondi_poly_fom = true;
do_quad_based_sketched_fom = true;

%% Adding paths for accessing the functions
addpath(fullfile(pwd, 'LR_deflation'));
addpath(fullfile(pwd, 'left_polynomial_preconditioned_FOM'));
addpath(fullfile(pwd, 'right_polynomial_preconditioned_FOM'));
addpath(fullfile(pwd, 'Quad_based_sketched_FOM'));

%% Define test parameters
rng(2130); % setting random seed generator for reproducibility

A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
N = size(A, 2); % Size of the matrix
gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12),gamma5hat);
A = Gamma5*A;

b = randn(N, 1); % Generate a random N x 1 vector

m = 5; % Define the number of critical eigenvalues

k1 = 3; % No. of iterations for the Krylov's subspace to be used in pre-conditioning polynomial

s = 500; % Sketch matrix row dimension

%% Compute f(A)x directly using the sign function
% exact_result = (A*(inv(sqrtm(full(A * A)))))*b;
% Save the value to exact_result.mat file
% save('exact_result.mat', 'exact_result');
% Load the exact result
loadedData = load('exact_result.mat', 'exact_result');
exact_result = loadedData.exact_result;  % Extract the value from the structure

%% Define the range of k values
k_values = 10:10:150;

%% Initialize arrays to store relative errors
rel_err_lr_deflation = zeros(length(k_values), 1);
rel_err_left_precondi_poly_fom = zeros(length(k_values), 1);
rel_err_right_precondi_poly_fom = zeros(length(k_values), 1);
rel_err_quad_based_sketched_fom = zeros(length(k_values), 1);

%% Invoking various functions to compute the product of the sign matrix of A and  b.
if do_lr_deflation
    start = cputime;

    % Compute f(A)x using LR_deflation
    fA_b = lr_deflation(A, b, m, k_values);

    finish = cputime;
    disp(['Time taken by lr deflation = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values)
        rel_err_lr_deflation(i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end 

if do_left_precondi_poly_fom
    start = cputime;

    % Compute f(A)x using left_precondi_poly_FOM
    fA_b = left_precondi_poly_FOM(A, b, k_values, k1);

    finish = cputime;
    disp(['Time taken by Left preconditioned FOM = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values)
        rel_err_left_precondi_poly_fom(i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end

if do_left_precondi_poly_fom
    start = cputime;

    % Compute f(A)x using right_precondi_poly_FOM
    fA_b = right_precondi_poly_FOM(A, b, k_values, k1);

    finish = cputime;
    disp(['Time taken by Right preconditioned FOM = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values)
        rel_err_right_precondi_poly_fom(i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end

if do_quad_based_sketched_fom
    start = cputime;

    % Compute f(A)x using Quadrature_based_sketched_FOM
    fA_b = Quadrature_based_sketched_FOM(A, b, k_values, s);

    finish = cputime;
    disp(['Time taken by Quadrature based sketched FOM = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values)
        rel_err_quad_based_sketched_fom (i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end

%% Plotting the relative errors wrt the no.of matrix mvms
figure;
if do_lr_deflation
    semilogy(k_values, rel_err_lr_deflation, 'r-o', 'DisplayName', 'LR Deflation');
    hold on;
end

if do_left_precondi_poly_fom
    semilogy(k_values, rel_err_left_precondi_poly_fom, 'b-o', 'DisplayName', 'Left preconditioned FOM');
    hold on;
end

if do_right_precondi_poly_fom
    semilogy(k_values, rel_err_right_precondi_poly_fom, 'g-^', 'DisplayName', 'Right preconditioned FOM');
    hold on;
end

if do_quad_based_sketched_fom
    semilogy(k_values, rel_err_quad_based_sketched_fom, 'k-^', 'DisplayName', 'Quadrature based sketched FOM');
    hold on;
end
hold off;

xlabel('# mvms ');
ylabel('Relative Error');
title('Relative Error vs #mvms');
legend('show');
grid on;
