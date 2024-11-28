 clear;
clc;
close all;

%% Adding paths for accessing the functions and Matrices

addpath(fullfile(pwd, 'Matrix_A'));

%% Select which methods to test
do_lr_deflation = false;

do_left_precondi_poly_arnoldi = false;
do_right_precondi_poly_arnoldi = false;

do_quad_based_sketched_trun_arnoldi = false;

do_quad_based_restarted_arnoldi = false;
do_quad_based_Impl_restarted_arnoldi = false;

do_combo_LR_def_LPoly_precond = true;
do_combo_LR_def_RPoly_precond = true;

do_combo_LR_def_quad_sketched_trun_arnoldi = true;

do_combo_LR_def_quad_rest_arnoldi = true;

do_combo_Lp_precond_quad_Impl_rest_arnoldi = true;
do_combo_Rp_precond_quad_Impl_rest_arnoldi = true;

%% Select the matrix to be tested
do_4x4_Herm = false;
do_4x4_Non_Herm = false;
do_8x4_Non_Herm = true;
do_16x4_Non_Herm = false;

value_present = true; % Read the .mat file containing the exact solution
do_plot = false; % To plot the graphs for the mvms
do_save = true; % To save the R.E, mvms, k_values and m

%% Loading the Matrix A and defining the vector b for the Lattice size 4^4, Hermitian matrix
rng(2130); % setting random seed generator for reproducibility

if do_4x4_Herm
    A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
    N = size(A, 2); % Size of the matrix
    A(end, 1) = A(end, 1)+ 1e-7;

    % Use eigs to find the smallest eigenvalue
    smallest_eigenvalue = eigs(A, 1, 'smallestreal');
    shift = 0.99 * real(smallest_eigenvalue);
    A = A - shift*speye(N);

    gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
    Gamma5 = kron(speye(N/12),gamma5hat);
    A = Gamma5*A;

    % Calculate the conjugate transpose of A
    A_conjugate_transpose = A';
    
    % Check if A is equal to its conjugate transpose
    if isequal(A, A_conjugate_transpose)
        disp('The matrix is Hermitian.');
    else
        disp('The matrix is not Hermitian.');
    end

    b = randn(N, 1); % Generate a random N x 1 vector
    
elseif do_4x4_Non_Herm
    A = read_matrix('periodic_L4_b3.55_k0.137n0_1.mat'); % Read the input matrix from a file.
    N = size(A, 2); % Size of the matrix

    % Use eigs to find the smallest eigenvalue
    % smallest_eigenvalue = eigs(A, 1, 'smallestreal');
    % shift = 0.99 * real(smallest_eigenvalue);
    % A = A - shift*speye(N);

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
k_values = 10:10:150;
% k_values = 50;

%% Define test parameters

m = [0, 2, 4, 8, 16, 64, 128]; % Define the number of critical eigenvalues
% m = [0, 64];

k1 = 4; % No. of iterations for the Krylov's subspace to be used in pre-conditioning polynomial

s = 300; % Sketch matrix row dimension

trunc = 2; % Truncate orthogonalization to the last 'trunc' vector

min_decay = 0.95; % Set Error minimum decay rate for convergence

tol = 1e-12; % Set tolerance level

max_iter = 50; % Maximum no.of restarts for the Arnoldi process

thick_number = [2, 4, 8, 10]; % Number of target eigenvalues for implicit deflation

k2 = 2; % No. of times the preconditioned Arnoldi process has to be exceuted.

%% Initialize arrays to store relative errors
rel_err_lr_deflation = zeros(length(k_values)* length(m), 1);

rel_err_left_precondi_poly_arnoldi = zeros(length(k_values), 1);
rel_err_right_precondi_poly_arnoldi = zeros(length(k_values), 1);

rel_err_quad_based_sketched_trun_arnoldi = zeros(length(k_values), 1);

rel_err_quad_based_restarted_arnoldi = zeros(length(k_values), 1);
rel_err_quad_based_Impl_restarted_arnoldi = zeros(length(k_values), 1);

rel_err_combo_LR_def_LPoly_precond = zeros(length(k_values) * length(m), 1);
rel_err_combo_LR_def_RPoly_precond = zeros(length(k_values) * length(m), 1);

rel_err_combo_LR_def_quad_sketched_trun_arnoldi = zeros(length(k_values) * length(m), 1);

rel_err_combo_LR_def_quad_rest_arnoldi = zeros(length(k_values) * length(m), 1);

rel_err_combo_Lp_precond_quad_Impl_rest_arnoldi = zeros(length(k_values)*length(thick_number), 1);
rel_err_combo_Rp_precond_quad_Impl_rest_arnoldi = zeros(length(k_values)*length(thick_number), 1);


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
    elseif do_8x4_Non_Herm
        loadedData = load('8x4_Non_Herm.mat', 'exact_result');
        exact_result = loadedData.exact_result;  % Extract the value from the structure
    end

else
    % Save the value to exact_result.mat file
    disp('Saving...');
    if do_4x4_Herm
        exact_result = (A*(inv(sqrtm(full(A * A)))))*b;
        save('4x4_Herm.mat', 'exact_result');
        disp('Saved the exact result to 4x4_Herm.mat');

    elseif do_4x4_Non_Herm
        exact_result = (A*(inv(sqrtm(full(A * A)))))*b;
        save('4x4_Non_Herm.mat', 'exact_result');
        disp('Saved the exact result to 4x4_Non_Herm.mat');

    elseif do_8x4_Non_Herm
        addpath(fullfile(pwd, 'Combo_LR_def_LPoly_precond'));
        [exact_result, ~] = combo_LR_def_LPoly_precond(A, b, 2000, 200, 8);
        save('8x4_Non_Herm.mat', 'exact_result');
        disp('Saved the exact result to 8x4_Non_Herm.mat');
        rmpath(fullfile(pwd, 'Combo_LR_def_LPoly_precond'));
        
    elseif do_16x4_Non_Herm
        save('16x4_Non_Herm.mat', 'exact_result');
        disp('Saved the exact result to 16x4_Non_Herm.mat');
    end
end

%% Invoking various functions to compute the product of the sign matrix of A and  b.
%% LR_deflation
if do_lr_deflation
    addpath(fullfile(pwd, 'LR_deflation'));

    start = cputime;

    % Compute f(A)x using LR_deflation
    [fA_b, mvms_lr_deflation] = lr_deflation(A, b, m, k_values);

    finish = cputime;
    disp(['Time taken by lr deflation = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values) * length(m)
        rel_err_lr_deflation(i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
    rmpath(fullfile(pwd, 'LR_deflation'));
end 

%% left_precondi_poly_arnoldi
if do_left_precondi_poly_arnoldi
    addpath(fullfile(pwd, 'Left_polynomial_preconditioned_arnoldi'));

    start = cputime;

    % Compute f(A)x using left_precondi_poly_arnoldi
    [fA_b, mvms_left_precondi_poly_arnoldi] = left_precondi_poly_arnoldi(A, b, k_values, k1);

    finish = cputime;
    disp(['Time taken by Left preconditioned Arnoldi = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values)
        rel_err_left_precondi_poly_arnoldi(i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
    rmpath(fullfile(pwd, 'Left_polynomial_preconditioned_arnoldi'));
end

%% right_precondi_poly_arnoldi
if do_right_precondi_poly_arnoldi
    addpath(fullfile(pwd, 'Right_polynomial_preconditioned_arnoldi'));

    start = cputime;

    % Compute f(A)x using right_precondi_poly_arnoldi
    [fA_b, mvms_right_precondi_poly_arnoldi] = right_precondi_poly_arnoldi(A, b, k_values, k1);

    finish = cputime;
    disp(['Time taken by Right preconditioned Arnoldi = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values)
        rel_err_right_precondi_poly_arnoldi(i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
    rmpath(fullfile(pwd, 'Right_polynomial_preconditioned_arnoldi'));
end

%% Quadrature_based_sketched_trun_arnoldi
if do_quad_based_sketched_trun_arnoldi
    addpath(fullfile(pwd, 'Quad_based_sketched_trun_arnoldi'));
    start = cputime;

    % Compute f(A)x using Quadrature_based_sketched_trun_arnoldi
    [fA_b, mvms_quad_based_sketched_trun_arnoldi] = Quadrature_based_sketched_trun_arnoldi(A, b, k_values, s, trunc);

    finish = cputime;
    disp(['Time taken by Quadrature based sketched truncated Arnoldi = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values)
        rel_err_quad_based_sketched_trun_arnoldi (i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
    rmpath(fullfile(pwd, 'Quad_based_sketched_trun_arnoldi'));
end

%% Quad_based_restarted_arnoldi
if do_quad_based_restarted_arnoldi
    addpath(fullfile(pwd, 'Quad_based_restarted_arnoldi'));
    start = cputime;
    mvms_quad_based_restarted_arnoldi = zeros(length(k_values), 1);
    for i = 1:length(k_values)
        [quadrature_approximation, ~, ~, mvms_quad_based_restarted_arnoldi(i)] = Quad_based_restarted_arnoldi(A, b, k_values(i), max_iter, tol, min_decay);
        rel_err_quad_based_restarted_arnoldi(i) = norm(exact_result - quadrature_approximation) / norm(exact_result);
    end
    finish = cputime;
    disp(['Time taken by Quadrature based restarted Arnoldi = ', num2str(finish - start), ' s']);
    rmpath(fullfile(pwd, 'Quad_based_restarted_arnoldi'));
end

%% Quad_based_Impl_restarted_arnoldi
if do_quad_based_Impl_restarted_arnoldi
    addpath(fullfile(pwd, 'Quad_based_Impl_restarted_arnoldi'));

    start = cputime;
    mvms_quad_based_Impl_restarted_arnoldi = zeros(length(k_values), 1);
    for i = 1:length(k_values)
        [quadrature_approximation, ~, ~, mvms_quad_based_Impl_restarted_arnoldi(i)] = Quad_based_imp_rest_arnoldi(A, b,k_values(i), max_iter, thick_number, tol, min_decay);
        rel_err_quad_based_Impl_restarted_arnoldi(i) = norm(exact_result - quadrature_approximation) / norm(exact_result);
    end
    finish = cputime;
    disp(['Time taken by Quadrature based Implicit restarted Arnoldi = ', num2str(finish - start), ' s']);
    rmpath(fullfile(pwd, 'Quad_based_Impl_restarted_arnoldi'));
end

%% combo_LR_def_LPoly_precond
if do_combo_LR_def_LPoly_precond
    addpath(fullfile(pwd, 'Combo_LR_def_LPoly_precond'));

    start = cputime;

    % Compute f(A)x using combo_LR_def_LPoly_precond
    [fA_b, mvms_combo_LR_def_LPoly_precond] = combo_LR_def_LPoly_precond(A, b, m, k_values, k1);

    finish = cputime;
    disp(['Time taken by Combination of LR deflation and Left preconditioned Arnoldi = ', num2str(finish - start), ' s']);
    % Loop over the range of k values
    for i = 1:length(k_values) * length(m)
        rel_err_combo_LR_def_LPoly_precond(i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
    rmpath(fullfile(pwd, 'Combo_LR_def_LPoly_precond'));
end

%% combo_LR_def_RPoly_precond
if do_combo_LR_def_RPoly_precond
    addpath(fullfile(pwd, 'Combo_LR_def_RPoly_precond'));

    start = cputime;

    % Compute f(A)x using combo_LR_def_RPoly_precond
    [fA_b, mvms_combo_LR_def_RPoly_precond] = combo_LR_def_RPoly_precond(A, b, m, k_values, k1);

    finish = cputime;
    disp(['Time taken by Combination of LR deflation and Right preconditioned Arnoldi = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values) * length(m)
        rel_err_combo_LR_def_RPoly_precond(i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
    rmpath(fullfile(pwd, 'Combo_LR_def_RPoly_precond'));
end

%% combo_LR_def_quad_sketched_trun_arnoldi
if do_combo_LR_def_quad_sketched_trun_arnoldi
    addpath(fullfile(pwd, 'Combo_LR_def_quad_sketched_trun_arnoldi'));

    start = cputime;

    % Compute f(A)x using combo_LR_def_quad_sketched_trun_arnoldi
    [fA_b, mvms_combo_LR_def_quad_sketched_trun_arnoldi] = combo_LR_def_quad_sketched_trun_arnoldi(A, b, m, k_values, s, trunc);

    finish = cputime;
    disp(['Time taken by Combination of LR deflation and Quadrature based sketched truncated Arnoldi = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values) * length(m)
        rel_err_combo_LR_def_quad_sketched_trun_arnoldi (i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
    rmpath(fullfile(pwd, 'Combo_LR_def_quad_sketched_trun_arnoldi'));
end

%% combo_LR_def_quad_rest_arnoldi
if do_combo_LR_def_quad_rest_arnoldi
    addpath(fullfile(pwd, 'Combo_LR_def_quad_rest_arnoldi'));

    start = cputime;

    % Compute f(A)x using combo_LR_def_quad_rest_arnoldi
    [fA_b, mvms_combo_LR_def_quad_rest_arnoldi, restarts_combo_LR_def_quad_rest_arnoldi] = combo_LR_def_quad_rest_arnoldi(A, b, m, k_values, max_iter, tol, min_decay);

    finish = cputime;
    disp(['Time taken by Combination of LR deflation and Quadrature based Restarted Arnoldi = ', num2str(finish - start), ' s']);
    
    % Loop over the range of k values
    for i = 1:length(k_values) * length(m)
        rel_err_combo_LR_def_quad_rest_arnoldi (i) = norm(exact_result - fA_b(:,i)) / norm(exact_result);
    end
    rmpath(fullfile(pwd, 'Combo_LR_def_quad_rest_arnoldi'));
end

%% combo_Lp_precond_quad_Impl_rest_arnoldi
if do_combo_Lp_precond_quad_Impl_rest_arnoldi
    addpath(fullfile(pwd, 'Combo_Lp_precond_quad_Impl_rest_arnoldi'));
    
    fA_b = [];
    mvms_combo_Lp_precond_quad_Impl_rest_arnoldi = [];
    restarts_combo_Lp_precond_quad_Impl_rest_arnoldi = [];

    start = cputime;

    for j = 1:length(thick_number)
        no = thick_number(j);
        
        % Call the function and get results
        [fAb, mvms, r] = Combo_Lp_precond_quad_Impl_rest_arnoldi(A, b, k_values, k1, k2, max_iter, no, tol, min_decay);
        
        % Concatenate results
        fA_b = [fA_b, fAb]; % Concatenate fA_b
        mvms_combo_Lp_precond_quad_Impl_rest_arnoldi = [mvms_combo_Lp_precond_quad_Impl_rest_arnoldi; mvms]; % Append mvms
        restarts_combo_Lp_precond_quad_Impl_rest_arnoldi = [restarts_combo_Lp_precond_quad_Impl_rest_arnoldi; r];
    end

    % Timing
    finish = cputime;
    disp(['Time taken by Combination of Left polynomial preconditioning and Quadrature based Implicit restarted Arnoldi = ', num2str(finish - start), ' s']);

    % Compute relative error for each value of fA_b
    for j = 1:length(thick_number)
        for i = 1:length(k_values)
            rel_err_combo_Lp_precond_quad_Impl_rest_arnoldi((j-1)*length(k_values)+i) = norm(exact_result - fA_b(:,(j-1)*length(k_values)+i)) / norm(exact_result);
        end
    end
    rmpath(fullfile(pwd, 'Combo_Lp_precond_quad_Impl_rest_arnoldi'));
end


%% combo_Rp_precond_quad_Impl_rest_arnoldi
if do_combo_Rp_precond_quad_Impl_rest_arnoldi
    addpath(fullfile(pwd, 'Combo_Rp_precond_quad_Impl_rest_arnoldi'));

    fA_b = [];
    mvms_combo_Rp_precond_quad_Impl_rest_arnoldi = [];
    restarts_combo_Rp_precond_quad_Impl_rest_arnoldi = [];
    
    start = cputime;

    for j = 1:length(thick_number)
        no = thick_number(j);
        
        % Call the function and get results
        [fAb, mvms] = Combo_Rp_precond_quad_Impl_rest_arnoldi(A, b, k_values, k1, k2, max_iter, no, tol, min_decay);
        
        % Concatenate results
        fA_b = [fA_b, fAb]; % Concatenate fA_b
        mvms_combo_Rp_precond_quad_Impl_rest_arnoldi = [mvms_combo_Rp_precond_quad_Impl_rest_arnoldi; mvms]; % Append mvms
        restarts_combo_Rp_precond_quad_Impl_rest_arnoldi = [restarts_combo_Rp_precond_quad_Impl_rest_arnoldi; r];
    end

    % Timing
    finish = cputime;
    disp(['Time taken by Combination of Right polynomial preconditioning and Quadrature based Implicit restarted Arnoldi = ', num2str(finish - start), ' s']);

    % Compute relative error for each value of fA_b
    for j = 1:length(thick_number)
        for i = 1:length(k_values)
            rel_err_combo_Rp_precond_quad_Impl_rest_arnoldi((j-1)*length(k_values)+i) = norm(exact_result - fA_b(:,(j-1)*length(k_values)+i)) / norm(exact_result);
        end
    end
    rmpath(fullfile(pwd, 'Combo_Rp_precond_quad_Impl_rest_arnoldi'));
end


%% Plotting the relative errors wrt the no.of matrix mvms
if do_plot

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
    
    if do_left_precondi_poly_arnoldi
        semilogy(mvms_left_precondi_poly_arnoldi, rel_err_left_precondi_poly_arnoldi, 'r-o', 'DisplayName', 'Left preconditioned Arnoldi');
        hold on;
    end
    
    if do_right_precondi_poly_arnoldi
        semilogy(mvms_right_precondi_poly_arnoldi, rel_err_right_precondi_poly_arnoldi, 'b-o', 'DisplayName', 'Right preconditioned Arnoldi');
        hold on;
    end
    
    if do_quad_based_sketched_trun_arnoldi
        semilogy(mvms_quad_based_sketched_trun_arnoldi, rel_err_quad_based_sketched_trun_arnoldi, '-o', 'DisplayName', 'Quadrature based sketched truncated Arnoldi');
        hold on;
    end
    
    if do_quad_based_restarted_arnoldi
        semilogy(mvms_quad_based_restarted_arnoldi, rel_err_quad_based_restarted_arnoldi, 'k-o', 'DisplayName', 'Quadrature based restarted Arnoldi');
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
            display_name = sprintf('Combination of LR deflation and Left preconditioned Arnoldi (m = %d)', m(i));
            
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
            display_name = sprintf('Combination of LR deflation and Right preconditioned Arnoldi (m = %d)', m(i));
            
            % Plot using the dynamic display name
            semilogy(mvms_combo_LR_def_RPoly_precond(j:j+length(k_values)-1), ...
                     rel_err_combo_LR_def_RPoly_precond(j:j+length(k_values)-1), ...
                     'b-*', 'DisplayName', display_name);
            hold on;
            
            % Update index j for the next segment of data
            j = j + length(k_values);
        end
    end
    
    if do_combo_LR_def_quad_sketched_trun_arnoldi
        j = 1;
        for i = 1:length(m)
            % Create a dynamic display name that includes the value of m(i)
            display_name = sprintf('Combination of LR deflation and Quadrature based sketched truncated Arnoldi (m = %d)', m(i));
            
            % Plot using the dynamic display name and magenta color
            semilogy(mvms_combo_LR_def_quad_sketched_trun_arnodli(j:j+length(k_values)-1), ...
                     rel_err_combo_LR_def_quad_sketched_trun_arnoldi(j:j+length(k_values)-1), ...
                     'm-*', 'DisplayName', display_name);
            hold on;
            
            % Update index j for the next segment of data
            j = j + length(k_values);
        end
    end
    
    if do_combo_LR_def_quad_rest_arnoldi
        j = 1;
        for i = 1:length(m)
            % Create a dynamic display name that includes the value of m(i)
            display_name = sprintf('Combination of LR deflation and Quadrature based restarted Arnoldi (m = %d)', m(i));
            
            % Plot using the dynamic display name and black color
            semilogy(mvms_combo_LR_def_quad_rest_arnoldi(j:j+length(k_values)-1), ...
                     rel_err_combo_LR_def_quad_rest_arnoldi(j:j+length(k_values)-1), ...
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
end

%% Save the results {R.E, mvms, k_values and m}
if do_save

    if do_lr_deflation
        save(fullfile('results', 'lr_deflation_results.mat'), 'k_values', 'rel_err_lr_deflation', 'mvms_lr_deflation', 'm');
    end
    
    if do_left_precondi_poly_arnoldi
        save(fullfile('results', 'left_precondi_poly_arnoldi_results.mat'), 'k_values', 'rel_err_left_precondi_poly_arnoldi', 'mvms_left_precondi_poly_arnoldi', 'm');
    end
    
    if do_right_precondi_poly_arnoldi
        save(fullfile('results', 'right_precondi_poly_arnoldi_results.mat'), 'k_values', 'rel_err_right_precondi_poly_arnoldi', 'mvms_right_precondi_poly_arnoldi', 'm');
    end
    
    if do_quad_based_sketched_trun_arnoldi
        save(fullfile('results', 'quad_based_sketched_trun_arnoldi_results.mat'), 'k_values', 'rel_err_quad_based_sketched_trun_arnoldi', 'mvms_quad_based_sketched_trun_arnoldi', 'm');
    end
    
    if do_quad_based_restarted_arnoldi
        save(fullfile('results', 'quad_based_restarted_arnoldi_results.mat'), 'k_values', 'rel_err_quad_based_restarted_arnoldi', 'mvms_quad_based_restarted_arnoldi', 'm');
    end
    
    if do_combo_LR_def_LPoly_precond
        save(fullfile('results', 'combo_LR_def_LPoly_precond_results.mat'), 'k_values', 'rel_err_combo_LR_def_LPoly_precond', 'mvms_combo_LR_def_LPoly_precond', 'm');
    end
    
    if do_combo_LR_def_RPoly_precond
        save(fullfile('results', 'combo_LR_def_RPoly_precond_results.mat'), 'k_values', 'rel_err_combo_LR_def_RPoly_precond', 'mvms_combo_LR_def_RPoly_precond', 'm');
    end
    
    if do_combo_LR_def_quad_sketched_trun_arnoldi
        save(fullfile('results', 'combo_LR_def_quad_sketched_trun_arnoldi_results.mat'), 'k_values', 'rel_err_combo_LR_def_quad_sketched_trun_arnoldi', 'mvms_combo_LR_def_quad_sketched_trun_arnoldi', 'm');
    end
    
    if do_combo_LR_def_quad_rest_arnoldi
        save(fullfile('results', 'combo_LR_def_quad_rest_arnoldi_results.mat'), 'k_values', 'rel_err_combo_LR_def_quad_rest_arnoldi', 'mvms_combo_LR_def_quad_rest_arnoldi', 'm', 'restarts_combo_LR_def_quad_rest_arnoldi');
    end
    
    if do_combo_Lp_precond_quad_Impl_rest_arnoldi
        m = thick_number;
        save(fullfile('results', 'combo_Lp_precond_quad_Impl_rest_arnoldi_results.mat'), 'k_values', 'rel_err_combo_Lp_precond_quad_Impl_rest_arnoldi', 'mvms_combo_Lp_precond_quad_Impl_rest_arnoldi', 'm', 'restarts_combo_Lp_precond_quad_Impl_rest_arnoldi');
    end
    
    if do_combo_Rp_precond_quad_Impl_rest_arnoldi
        m = thick_number;
        save(fullfile('results', 'combo_Rp_precond_quad_Impl_rest_arnoldi_results.mat'), 'k_values', 'rel_err_combo_Rp_precond_quad_Impl_rest_arnoldi', 'mvms_combo_Rp_precond_quad_Impl_rest_arnoldi', 'm', 'restarts_combo_Rp_precond_quad_Impl_rest_arnoldi');
    end
end