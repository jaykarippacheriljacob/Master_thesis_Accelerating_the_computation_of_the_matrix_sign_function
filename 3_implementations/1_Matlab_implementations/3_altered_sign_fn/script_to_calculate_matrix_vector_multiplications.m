clear;
clc;
close all;

%% Adding paths for accessing the functions

addpath(fullfile(pwd, 'Matrix_A'));
addpath(fullfile(pwd, 'LR_deflation'));

addpath(fullfile(pwd, 'Left_polynomial_preconditioned_FOM'));
addpath(fullfile(pwd, 'Right_polynomial_preconditioned_FOM'));

addpath(fullfile(pwd, 'Quad_based_sketched_FOM'));
addpath(fullfile(pwd, 'Quad_based_sketched_trun_FOM'));

addpath(fullfile(pwd, 'Quad_based_Expl_restarted_arnoldi'));
addpath(fullfile(pwd, 'Quad_based_Impl_restarted_arnoldi'));

addpath(fullfile(pwd, 'Combo_LR_def_LPoly_precond'));
addpath(fullfile(pwd, 'Combo_LR_def_RPoly_precond'));

addpath(fullfile(pwd, 'Combo_LR_def_quad_sketched_FOM'));
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

do_quad_based_Expl_restarted_arnoldi = false;
do_quad_based_Impl_restarted_arnoldi = false;

do_combo_LR_def_LPoly_precond = false;
do_combo_LR_def_RPoly_precond = false;

do_combo_LR_def_quad_sketched_FOM = false;
do_combo_LR_def_quad_sketched_trun_FOM = false;

do_combo_LR_def_quad_expl_rest_arnoldi = false;

do_combo_Lp_precond_quad_Impl_rest_arnoldi = true;
do_combo_Rp_precond_quad_Impl_rest_arnoldi = false;

%% Select the matrix to be tested
do_4x4_Herm = false;
do_4x4_Non_Herm = true;
do_8x4_Non_Herm = false;
do_16x4_Non_Herm = false;

value_present = true; % Read the .mat file containing the exact solution

%% Loading the Matrix A and defining the vector b for the Lattice size 4^4, Hermitian matrix
rng(2130); % setting random seed generator for reproducibility

if do_4x4_Herm
    A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
    N = size(A, 2); % Size of the matrix
    gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
    Gamma5 = kron(speye(N/12),gamma5hat);
    A = Gamma5*A;

    b = randn(N, 1); % Generate a random N x 1 vector
    
elseif do_4x4_Non_Herm
    A = read_matrix('periodic_L4_b3.55_k0.137n0_1.mat'); % Read the input matrix from a file.
    N = size(A, 2); % Size of the matrix
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

%% Define test parameters

m = 50; % Define the number of critical eigenvalues

k1 = 3; % No. of iterations for the Krylov's subspace to be used in pre-conditioning polynomial

s = 500; % Sketch matrix row dimension

trunc = 19; % Truncate orthogonalization to the last 'trunc' vector

min_decay = 0.95; % Set Error minimum decay rate for convergence

tol = 1e-8; % Set tolerance level

max_iter = 50; % Maximum no.of restarts for the Arnoldi process

thick_number = 5; % Number of target eigenvalues for implicit deflation

k2 = 1; % No. of times the preconditioned Arnoldi process has to be exceuted.

%% Initialize arrays to store relative errors
rel_err_lr_deflation = zeros(length(k_values), 1);

rel_err_left_precondi_poly_fom = zeros(length(k_values), 1);
rel_err_right_precondi_poly_fom = zeros(length(k_values), 1);

rel_err_quad_based_sketched_fom = zeros(length(k_values), 1);
rel_err_quad_based_sketched_trun_fom = zeros(length(k_values), 1);

rel_err_quad_based_Expl_restarted_arnoldi = zeros(length(k_values), 1);
rel_err_quad_based_Impl_restarted_arnoldi = zeros(length(k_values), 1);

rel_err_combo_LR_def_LPoly_precond = zeros(length(k_values), 1);
rel_err_combo_LR_def_RPoly_precond = zeros(length(k_values), 1);

rel_err_combo_LR_def_quad_sketched_FOM = zeros(length(k_values), 1);
rel_err_combo_LR_def_quad_sketched_trun_FOM = zeros(length(k_values), 1);

rel_err_combo_LR_def_quad_expl_rest_arnoldi = zeros(length(k_values), 1);

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
    for i = 1:length(k_values)
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

%% Quadrature_based_sketched_FOM
if do_quad_based_sketched_fom
    start = cputime;

    % Compute f(A)x using Quadrature_based_sketched_FOM
    [fA_b, mvms_quad_based_sketched_fom] = Quadrature_based_sketched_FOM(A, b, k_values, s);

    finish = cputime;
    disp(['Time taken by Quadrature based sketched FOM = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values)
        rel_err_quad_based_sketched_fom (i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
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
    
    % Loop over the range of k values
    for i = 1:length(k_values)
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
    for i = 1:length(k_values)
        rel_err_combo_LR_def_RPoly_precond(i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end

%% combo_LR_def_quad_sketched_FOM
if do_combo_LR_def_quad_sketched_FOM
    start = cputime;

    % Compute f(A)x using combo_LR_def_quad_sketched_FOM
    [fA_b, mvms_combo_LR_def_quad_sketched_FOM] = combo_LR_def_quad_sketched_FOM(A, b, m, k_values, s);

    finish = cputime;
    disp(['Time taken by Combination of LR deflation and Quadrature based sketched FOM = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values)
        rel_err_combo_LR_def_quad_sketched_FOM(i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
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
    for i = 1:length(k_values)
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
    for i = 1:length(k_values)
        rel_err_combo_LR_def_quad_expl_rest_arnoldi (i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
end

%% combo_Lp_precond_quad_Impl_rest_arnoldi
if do_combo_Lp_precond_quad_Impl_rest_arnoldi
    start = cputime;

    % Compute f(A)x using combo_Lp_precond_quad_Impl_rest_arnoldi
    % [fA_b, mvms_combo_Lp_precond_quad_Impl_rest_arnoldi] = Combo_Lp_precond_quad_Impl_rest_arnoldi(A, b, k_values, k1, k2, max_iter, thick_number, tol, min_decay);
    [fA_b] = Combo_Lp_precond_quad_Impl_rest_arnoldi(A, b, k_values, k1, k2, max_iter, thick_number, tol, min_decay);
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
    semilogy(mvms_lr_deflation, rel_err_lr_deflation, 'r-o', 'DisplayName', 'LR Deflation');
    hold on;
end

if do_left_precondi_poly_fom
    semilogy(mvms_left_precondi_poly_fom, rel_err_left_precondi_poly_fom, 'b-o', 'DisplayName', 'Left preconditioned FOM');
    hold on;
end

if do_right_precondi_poly_fom
    semilogy(mvms_right_precondi_poly_fom, rel_err_right_precondi_poly_fom, 'g-^', 'DisplayName', 'Right preconditioned FOM');
    hold on;
end

if do_quad_based_sketched_fom
    semilogy(mvms_quad_based_sketched_fom, rel_err_quad_based_sketched_fom, 'k-^', 'DisplayName', 'Quadrature based sketched FOM');
    hold on;
end

if do_quad_based_sketched_trun_fom
    semilogy(mvms_quad_based_sketched_trun_fom, rel_err_quad_based_sketched_trun_fom, 'k-o', 'DisplayName', 'Quadrature based sketched truncated FOM');
    hold on;
end

if do_quad_based_Expl_restarted_arnoldi
    semilogy(mvms_quad_based_Expl_restarted_arnoldi, rel_err_quad_based_Expl_restarted_arnoldi, 'b-o', 'DisplayName', 'Quadrature based Explicit restarted Arnoldi');
    hold on;
end

if do_quad_based_Impl_restarted_arnoldi
    semilogy(mvms_quad_based_Impl_restarted_arnoldi, rel_err_quad_based_Impl_restarted_arnoldi, 'r-^', 'DisplayName', 'Quadrature based Implicit restarted Arnoldi');
    hold on;
end

if do_combo_LR_def_LPoly_precond
    semilogy(mvms_combo_LR_def_LPoly_precond, rel_err_combo_LR_def_LPoly_precond, 'r-o', 'DisplayName', 'Combination of LR deflation and Left preconditioned FOM');
    hold on;
end

if do_combo_LR_def_RPoly_precond
    semilogy(mvms_combo_LR_def_RPoly_precond, rel_err_combo_LR_def_RPoly_precond, 'r-^', 'DisplayName', 'Combination of LR deflation and Right preconditioned FOM');
    hold on;
end

if do_combo_LR_def_quad_sketched_FOM
    semilogy(mvms_combo_LR_def_quad_sketched_FOM, rel_err_combo_LR_def_quad_sketched_FOM, 'y-o', 'DisplayName', 'Combination of LR deflation and Quadrature based sketched FOM');
    hold on;
end

if do_combo_LR_def_quad_sketched_trun_FOM
    semilogy(mvms_combo_LR_def_quad_sketched_trun_FOM, rel_err_combo_LR_def_quad_sketched_trun_FOM, 'y-^', 'DisplayName', 'Combination of LR deflation and Quadrature based sketched truncated FOM');
    hold on;
end

if do_combo_LR_def_quad_expl_rest_arnoldi
    semilogy(mvms_combo_LR_def_quad_expl_rest_arnoldi, rel_err_combo_LR_def_quad_expl_rest_arnoldi, 'g-^', 'DisplayName', 'Combination of LR deflation and Quadrature based Explicit restarted Arnoldi');
    hold on;
end

if do_combo_Lp_precond_quad_Impl_rest_arnoldi
    semilogy(mvms_combo_Lp_precond_quad_Impl_rest_arnoldi, rel_err_combo_Lp_precond_quad_Impl_rest_arnoldi, 'g-^', 'DisplayName', 'Combination of Left polynomial preconditioning and Quadrature based Implicit restarted Arnoldi');
    hold on;
end

if do_combo_Rp_precond_quad_Impl_rest_arnoldi
    semilogy(mvms_combo_Rp_precond_quad_Impl_rest_arnoldi, rel_err_combo_Rp_precond_quad_Impl_rest_arnoldi, 'g-^', 'DisplayName', 'Combination of Right polynomial preconditioning and Quadrature based Implicit restarted Arnoldi');
    hold on;
end
hold off;

xlabel('# mvms ');
ylabel('Relative Error');
title('Relative Error vs #mvms');
legend('show');
grid on;
