clear;
clc;
close all;

%% Adding paths for accessing the functions

addpath(fullfile(pwd, 'Matrix_A'));
addpath(fullfile(pwd, 'LR_deflation'));

addpath(fullfile(pwd, 'Left_polynomial_preconditioned_FOM'));
addpath(fullfile(pwd, 'Right_polynomial_preconditioned_FOM'));

addpath(fullfile(pwd, 'Quad_based_sketched_trun_FOM'));

addpath(fullfile(pwd, 'Quad_based_Expl_restarted_arnoldi'));
addpath(fullfile(pwd, 'Quad_based_Impl_restarted_arnoldi'));

addpath(fullfile(pwd, 'Combo_LR_def_LPoly_precond'));
% addpath(fullfile(pwd, 'Combo_LR_def_LPoly_precond_Copy'));
addpath(fullfile(pwd, 'Combo_LR_def_RPoly_precond'));

addpath(fullfile(pwd, 'Combo_LR_def_quad_sketched_trun_FOM'));

addpath(fullfile(pwd, 'Combo_LR_def_quad_expl_rest_arnoldi'));

addpath(fullfile(pwd, 'Combo_Lp_precond_quad_Impl_rest_arnoldi'));
addpath(fullfile(pwd, 'Combo_Rp_precond_quad_Impl_rest_arnoldi'));

%% Select which methods to test
do_lr_deflation = false;

do_left_precondi_poly_fom = false;
do_right_precondi_poly_fom = false;

do_quad_based_sketched_fom = false;
do_quad_based_sketched_trun_fom = false;

do_quad_based_Expl_restarted_arnoldi = true;
do_quad_based_Impl_restarted_arnoldi = true;

do_combo_LR_def_LPoly_precond = false;
do_combo_LR_def_RPoly_precond = false;

do_combo_LR_def_quad_sketched_trun_FOM = false;

do_combo_LR_def_quad_expl_rest_arnoldi = true;

do_combo_Lp_precond_quad_Impl_rest_arnoldi = true;
do_combo_Rp_precond_quad_Impl_rest_arnoldi = false;

%% Select the matrix to be tested
do_4x4_Herm = true;
do_4x4_Non_Herm = false;
do_8x4_Non_Herm = false;
do_16x4_Non_Herm = false;

value_present = true; % Read the .mat file containing the exact solution

%% Loading the Matrix A and defining the vector b for the Lattice size 4^4, Hermitian matrix
rng(2130); % setting random seed generator for reproducibility

if do_4x4_Herm
    A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
    N = size(A, 2); % Size of the matrix

    % Use eigs to find the smallest eigenvalue
    smallest_eigenvalue = eigs(A, 1, 'smallestreal');
    shift = 0.99 * real(smallest_eigenvalue);
    A = A - shift*speye(N);

    gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
    Gamma5 = kron(speye(N/12),gamma5hat);
    A = Gamma5*A;

    b = randn(N, 1); % Generate a random N x 1 vector
    
elseif do_4x4_Non_Herm
    A = read_matrix('periodic_L4_b3.55_k0.137n0_1.mat'); % Read the input matrix from a file.
    N = size(A, 2); % Size of the matrix

    % Use eigs to find the smallest eigenvalue
    smallest_eigenvalue = eigs(A, 1, 'smallestreal');
    shift = 0.99 * real(smallest_eigenvalue);
    A = A - shift*speye(N);

    gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
    Gamma5 = kron(speye(N/12),gamma5hat);
    A = Gamma5*A;

    b = randn(N, 1); % Generate a random N x 1 vector

elseif do_8x4_Non_Herm
    A = read_matrix('periodic_L8_b3.55_k0.137n0_1.mat'); % Read the input matrix from a file.
    N = size(A, 2); % Size of the matrix
    gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
    Gamma5 = kron(speye(N/12),gamma5hat);
    A = Gamma5*A;

    b = randn(N, 1); % Generate a random N x 1 vector

elseif do_16x4_Non_Herm
    A = read_matrix('periodic_L16_b3.55_k0.137n0_1.mat'); % Read the input matrix from a file.
    N = size(A, 2); % Size of the matrix
    gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
    Gamma5 = kron(speye(N/12),gamma5hat);
    A = Gamma5*A;

    b = randn(N, 1); % Generate a random N x 1 vector
end


%% Define the range of k values
k_values = 20:10:150;
% k_values = 50;

%% Define test parameters

m = [0, 2, 4, 8, 16, 32, 64]; % Define the number of critical eigenvalues
% m = 8;

k1 = 8; % No. of iterations for the Krylov's subspace to be used in pre-conditioning polynomial

s = 500; % Sketch matrix row dimension

trunc = 19; % Truncate orthogonalization to the last 'trunc' vector

min_decay = 0.95; % Set Error minimum decay rate for convergence

tol = 1e-10; % Set tolerance level

max_iter = 50; % Maximum no.of restarts for the Arnoldi process

thick_number = 8; % Number of target eigenvalues for implicit deflation

k2 = 1; % No. of times the preconditioned Arnoldi process has to be exceuted.

%% Initialize arrays to store relative errors
rel_err_lr_deflation = zeros(length(k_values), 1);

rel_err_left_precondi_poly_fom = zeros(length(k_values), 1);
rel_err_right_precondi_poly_fom = zeros(length(k_values), 1);

rel_err_quad_based_sketched_trun_fom = zeros(length(k_values), 1);

rel_err_quad_based_Expl_restarted_arnoldi = zeros(length(k_values), 1);
rel_err_quad_based_Impl_restarted_arnoldi = zeros(length(k_values), 1);

rel_err_combo_LR_def_LPoly_precond = zeros(length(k_values) * length(m), 1);
rel_err_combo_LR_def_RPoly_precond = zeros(length(k_values) * length(m), 1);

rel_err_combo_LR_def_quad_sketched_trun_FOM = zeros(length(k_values) * length(m), 1);

rel_err_combo_LR_def_quad_expl_rest_arnoldi = zeros(length(k_values) * length(m), 1);

rel_err_combo_Lp_precond_quad_Impl_rest_arnoldi = zeros(length(k_values), 1);
rel_err_combo_Rp_precond_quad_Impl_rest_arnoldi = zeros(length(k_values), 1);


%% Compute f(A)x directly using the sign function
% Load the exact result
%% For the Lattice size 4^4, Hermitian matrix
if value_present
    if do_4x4_Herm
        loadedData = load('4x4_Herm.mat', 'exact_result');
        exact_result = loadedData.exact_result;  % Extract the value from the structure

   elseif do_4x4_Non_Herm
        loadedData = load('4x4_Non_Herm.mat', 'exact_result');
        exact_result = loadedData.exact_result;  % Extract the value from the structure
    end

else
    exact_result = (A*(inv(sqrtm(full(A * A)))))*b;
    % Save the value to exact_result.mat file
    disp('Saving...');
    if do_4x4_Herm
        save('4x4_Herm.mat', 'exact_result');
        disp('Saved the exact result to 4x4_Herm.mat');

    elseif do_4x4_Non_Herm
        save('4x4_Non_Herm.mat', 'exact_result');
        disp('Saved the exact result to 4x4_Non_Herm.mat');

    elseif do_8x4_Non_Herm
        save('8x4_Non_Herm.mat', 'exact_result');
        disp('Saved the exact result to 8x4_Non_Herm.mat');
        
    elseif do_16x4_Non_Herm
        save('16x4_Non_Herm.mat', 'exact_result');
        disp('Saved the exact result to 16x4_Non_Herm.mat');
    end
end

%% Invoking various functions to compute the product of the sign matrix of A and  b.
%% LR_deflation
if do_lr_deflation
    start = cputime;

    % Compute f(A)x using LR_deflation
    [fA_b, mvms_lr_deflation] = lr_deflation(A, b, m, k_values);

    finish = cputime;
    disp(['Time taken by lr deflation = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values) * length(m)
        rel_err_lr_deflation(i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end 

%% left_precondi_poly_FOM
if do_left_precondi_poly_fom
    start = cputime;

    % Compute f(A)x using left_precondi_poly_FOM
    [fA_b, mvms_left_precondi_poly_fom] = left_precondi_poly_FOM(A, b, k_values, k1);

    finish = cputime;
    disp(['Time taken by Left preconditioned FOM = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values)
        rel_err_left_precondi_poly_fom(i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end

%% right_precondi_poly_FOM
if do_right_precondi_poly_fom
    start = cputime;

    % Compute f(A)x using right_precondi_poly_FOM
    [fA_b, mvms_right_precondi_poly_fom] = right_precondi_poly_FOM(A, b, k_values, k1);

    finish = cputime;
    disp(['Time taken by Right preconditioned FOM = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values)
        rel_err_right_precondi_poly_fom(i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end

%% Quadrature_based_sketched_trun_FOM
if do_quad_based_sketched_trun_fom
    start = cputime;

    % Compute f(A)x using Quadrature_based_sketched_trun_FOM
    [fA_b, mvms_quad_based_sketched_trun_fom] = Quadrature_based_sketched_trun_FOM(A, b, k_values, s, trunc);

    finish = cputime;
    disp(['Time taken by Quadrature based sketched truncated FOM = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values)
        rel_err_quad_based_sketched_trun_fom (i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end

%% Quad_based_Expl_restarted_arnoldi
if do_quad_based_Expl_restarted_arnoldi
    start = cputime;
    mvms_quad_based_Expl_restarted_arnoldi = zeros(length(k_values), 1);
    for i = 1:length(k_values)
        [quadrature_approximation, ~, ~, mvms_quad_based_Expl_restarted_arnoldi(i)] = Quad_based_Exp_restarted_arnoldi(A, b, k_values(i), max_iter, tol, min_decay);
        rel_err_quad_based_Expl_restarted_arnoldi(i) = norm(exact_result - quadrature_approximation) / norm(exact_result);
    end
    finish = cputime;
    disp(['Time taken by Quadrature based Explicit restarted Arnoldi = ', num2str(finish - start), ' s']);
end

%% Quad_based_Impl_restarted_arnoldi
if do_quad_based_Impl_restarted_arnoldi
    start = cputime;
    mvms_quad_based_Impl_restarted_arnoldi = zeros(length(k_values), 1);
    for i = 1:length(k_values)
        [quadrature_approximation, ~, ~, mvms_quad_based_Impl_restarted_arnoldi(i)] = Quad_based_imp_rest_arnoldi(A, b,k_values(i), max_iter, thick_number, tol, min_decay);
        rel_err_quad_based_Impl_restarted_arnoldi(i) = norm(exact_result - quadrature_approximation) / norm(exact_result);
    end
    finish = cputime;
    disp(['Time taken by Quadrature based Implicit restarted Arnoldi = ', num2str(finish - start), ' s']);
end

%% combo_LR_def_LPoly_precond
if do_combo_LR_def_LPoly_precond
    start = cputime;

    % Compute f(A)x using combo_LR_def_LPoly_precond
    [fA_b, mvms_combo_LR_def_LPoly_precond] = combo_LR_def_LPoly_precond(A, b, m, k_values, k1);

    finish = cputime;
    disp(['Time taken by Combination of LR deflation and Left preconditioned FOM = ', num2str(finish - start), ' s']);
    exact_result = fA_b(:, end);
    % Loop over the range of k values
    for i = 1:length(k_values) * length(m)
        rel_err_combo_LR_def_LPoly_precond(i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end

%% combo_LR_def_RPoly_precond
if do_combo_LR_def_RPoly_precond
    start = cputime;

    % Compute f(A)x using combo_LR_def_RPoly_precond
    [fA_b, mvms_combo_LR_def_RPoly_precond] = combo_LR_def_RPoly_precond(A, b, m, k_values, k1);

    finish = cputime;
    disp(['Time taken by Combination of LR deflation and Right preconditioned FOM = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values) * length(m)
        rel_err_combo_LR_def_RPoly_precond(i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end

%% combo_LR_def_quad_sketched_trun_FOM
if do_combo_LR_def_quad_sketched_trun_FOM
    start = cputime;

    % Compute f(A)x using combo_LR_def_quad_sketched_trun_FOM
    [fA_b, mvms_combo_LR_def_quad_sketched_trun_FOM] = combo_LR_def_quad_sketched_trun_FOM(A, b, m, k_values, s, trunc);

    finish = cputime;
    disp(['Time taken by Combination of LR deflation and Quadrature based sketched truncated FOM = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values) * length(m)
        rel_err_combo_LR_def_quad_sketched_trun_FOM (i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end

%% combo_LR_def_quad_expl_rest_arnoldi
if do_combo_LR_def_quad_expl_rest_arnoldi
    start = cputime;

    % Compute f(A)x using combo_LR_def_quad_expl_rest_arnoldi
    [fA_b, mvms_combo_LR_def_quad_expl_rest_arnoldi] = combo_LR_def_quad_expl_rest_arnoldi(A, b, m, k_values, max_iter, tol, min_decay);

    finish = cputime;
    disp(['Time taken by Combination of LR deflation and Quadrature based Explicit Restarted Arnoldi = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values) * length(m)
        rel_err_combo_LR_def_quad_expl_rest_arnoldi (i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end

%% combo_Lp_precond_quad_Impl_rest_arnoldi
if do_combo_Lp_precond_quad_Impl_rest_arnoldi
    start = cputime;

    % Compute f(A)x using combo_Lp_precond_quad_Impl_rest_arnoldi
    [fA_b, mvms_combo_Lp_precond_quad_Impl_rest_arnoldi] = Combo_Lp_precond_quad_Impl_rest_arnoldi(A, b, k_values, k1, k2, max_iter, thick_number, tol, min_decay);
    % [fA_b] = Combo_Lp_precond_quad_Impl_rest_arnoldi(A, b, k_values, k1, k2, max_iter, thick_number, tol, min_decay);
    finish = cputime;
    disp(['Time taken by Combination of Left polynomial preconditioning and Quadrature based Implicit restarted Arnoldi = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values)
        rel_err_combo_Lp_precond_quad_Impl_rest_arnoldi (i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end

%% combo_Rp_precond_quad_Impl_rest_arnoldi
if do_combo_Rp_precond_quad_Impl_rest_arnoldi
    start = cputime;

    % Compute f(A)x using combo_Rp_precond_quad_Impl_rest_arnoldi
    [fA_b, mvms_combo_Rp_precond_quad_Impl_rest_arnoldi] = Combo_Rp_precond_quad_Impl_rest_arnoldi(A, b, k_values, k1, k2, max_iter, thick_number, tol, min_decay);

    finish = cputime;
    disp(['Time taken by Combination of Right polynomial preconditioning and Quadrature based Implicit restarted Arnoldi = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values)
        rel_err_combo_Rp_precond_quad_Impl_rest_arnoldi (i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end

%% Plotting the relative errors wrt the no.of matrix mvms
figure;
if do_lr_deflation
    j = 1;
    for i = 1:length(m)
        % Create a dynamic display name that includes the value of m(i)
        display_name = sprintf('LR Deflation (m = %d)', m(i));

        semilogy(mvms_lr_deflation(j:j+length(k_values)-1), rel_err_lr_deflation(j:j+length(k_values)-1), 'c-o', 'DisplayName', display_name);
        hold on;
        j = j+length(k_values);
    end
end

if do_left_precondi_poly_fom
    semilogy(mvms_left_precondi_poly_fom, rel_err_left_precondi_poly_fom, 'r-o', 'DisplayName', 'Left preconditioned FOM');
    hold on;
end

if do_right_precondi_poly_fom
    semilogy(mvms_right_precondi_poly_fom, rel_err_right_precondi_poly_fom, 'b-o', 'DisplayName', 'Right preconditioned FOM');
    hold on;
end

if do_quad_based_sketched_trun_fom
    semilogy(mvms_quad_based_sketched_trun_fom, rel_err_quad_based_sketched_trun_fom, '-o', 'DisplayName', 'Quadrature based sketched truncated FOM');
    hold on;
end

if do_quad_based_Expl_restarted_arnoldi
    semilogy(mvms_quad_based_Expl_restarted_arnoldi, rel_err_quad_based_Expl_restarted_arnoldi, 'k-o', 'DisplayName', 'Quadrature based Explicit restarted Arnoldi');
    hold on;
end

if do_quad_based_Impl_restarted_arnoldi
    semilogy(mvms_quad_based_Impl_restarted_arnoldi, rel_err_quad_based_Impl_restarted_arnoldi, '-o', 'DisplayName', 'Quadrature based Implicit restarted Arnoldi', 'Color', [1, 0.5, 0]);
    hold on;
end

if do_combo_LR_def_LPoly_precond
    j = 1;
    for i = 1:length(m)
        % Create a dynamic display name that includes the value of m(i)
        display_name = sprintf('Combination of LR deflation and Left preconditioned FOM (m = %d)', m(i));
        
        % Plot using the dynamic display name
        semilogy(mvms_combo_LR_def_LPoly_precond(j:j+length(k_values)-1), ...
                 rel_err_combo_LR_def_LPoly_precond(j:j+length(k_values)-1), ...
                 'r-*', 'DisplayName', display_name);
        hold on;
        
        % Update index j for the next segment of data
        j = j + length(k_values);
    end
end

if do_combo_LR_def_RPoly_precond
    j = 1;
    for i = 1:length(m)
        % Create a dynamic display name that includes the value of m(i)
        display_name = sprintf('Combination of LR deflation and Right preconditioned FOM (m = %d)', m(i));
        
        % Plot using the dynamic display name
        semilogy(mvms_combo_LR_def_RPoly_precond(j:j+length(k_values)-1), ...
                 rel_err_combo_LR_def_RPoly_precond(j:j+length(k_values)-1), ...
                 'b-*', 'DisplayName', display_name);
        hold on;
        
        % Update index j for the next segment of data
        j = j + length(k_values);
    end
end

if do_combo_LR_def_quad_sketched_trun_FOM
    j = 1;
    for i = 1:length(m)
        % Create a dynamic display name that includes the value of m(i)
        display_name = sprintf('Combination of LR deflation and Quadrature based sketched truncated FOM (m = %d)', m(i));
        
        % Plot using the dynamic display name and magenta color
        semilogy(mvms_combo_LR_def_quad_sketched_trun_FOM(j:j+length(k_values)-1), ...
                 rel_err_combo_LR_def_quad_sketched_trun_FOM(j:j+length(k_values)-1), ...
                 'm-*', 'DisplayName', display_name);
        hold on;
        
        % Update index j for the next segment of data
        j = j + length(k_values);
    end
end

if do_combo_LR_def_quad_expl_rest_arnoldi
    j = 1;
    for i = 1:length(m)
        % Create a dynamic display name that includes the value of m(i)
        display_name = sprintf('Combination of LR deflation and Quadrature based Explicit restarted Arnoldi (m = %d)', m(i));
        
        % Plot using the dynamic display name and black color
        semilogy(mvms_combo_LR_def_quad_expl_rest_arnoldi(j:j+length(k_values)-1), ...
                 rel_err_combo_LR_def_quad_expl_rest_arnoldi(j:j+length(k_values)-1), ...
                 'k-*', 'DisplayName', display_name);
        hold on;
        
        % Update index j for the next segment of data
        j = j + length(k_values);
    end
end


if do_combo_Lp_precond_quad_Impl_rest_arnoldi
    semilogy(mvms_combo_Lp_precond_quad_Impl_rest_arnoldi, rel_err_combo_Lp_precond_quad_Impl_rest_arnoldi, '-*', 'DisplayName', 'Combination of Left polynomial preconditioning and Quadrature based Implicit restarted Arnoldi', 'Color', [1, 0.5, 0]);
    % semilogy(k_values, rel_err_combo_Lp_precond_quad_Impl_rest_arnoldi, 'g-^', 'DisplayName', 'Combination of Left polynomial preconditioning and Quadrature based Implicit restarted Arnoldi');
    hold on;
end

if do_combo_Rp_precond_quad_Impl_rest_arnoldi
    semilogy(mvms_combo_Rp_precond_quad_Impl_rest_arnoldi, rel_err_combo_Rp_precond_quad_Impl_rest_arnoldi, '-*', 'DisplayName', 'Combination of Right polynomial preconditioning and Quadrature based Implicit restarted Arnoldi', 'Color', [0.6, 0.4, 0.2]);
    hold on;
end
hold off;

xlabel('# mvms ');
ylabel('Relative Error');
title('Relative Error vs #mvms');
legend('show');
grid on;
