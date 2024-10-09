function [Y, cost] = combo_LR_def_quad_expl_rest_arnoldi(A, x, m, k_values, max_iter, tol, min_decay)
    %% Combination of LR-deflation and Quadrture based Explicit restarted arnoldi approximation for f(A)b.
    % Input:
    %      A         - n x n matrix
    %      x         - n x 1 vector
    %      m         - No. of critical values for which defation has to be undergone.
    %      k_values  - vector of increasing Krylov subspace dimensions.
    %      max_iter  - Maximum no.of iterations for the restart of the Arnoldi decomposition
    %      tol       - Set tolerance for stopping criteria
    %      min_decay - the decay rate of error after each iteration.
    % Output:
    %      Y         - Approximation of f(A)b for k_values. We return the deflated 
    %                  left prec. Arnoldi matrix approximation for all
    %                  values in k_values dimensions.
    %      cost      - No.of matrix{A} vector multiplications
    
    addpath(fullfile(pwd, 'Combo_LR_def_quad_expl_rest_arnoldi'));

    %% Step 1: Compute left and right eigenvectors
    [Rm, Lm, Dm] = compute_eigenvectors(A, m);

    %% Step 2: Compute f(A) for critical eigenvalues
    f_Tm = compute_sign_function_diag(Dm);

    %% Step 3: Compute x_ominus = (1 − Rm * Lm' ) * x 
    % ' -> transpose conjugate
    % .' -> transpose
    x_ominus = compute_x_ominus(Rm, Lm, x);

    %% Step 4: Compute f(A)x_ominus using the restarted arnoldi process for k_values krylov subspace dimensions.
    no_k = length(k_values);
    n = size(A,2);
    fA_x_ominus = zeros(n,no_k);
    cost = zeros(no_k, 1);

    for l=1:no_k  %l-th column holds approx. for subspace dimension k(l)
        [fA_x_ominus(:,l), ~, ~,cost(l)] = Quad_based_Expl_restarted_arnoldi(A, x_ominus, k_values(l), max_iter, tol, min_decay);
    end

    %% Step 6: Compute the approximation to f(A)x

    Y = Rm * (f_Tm * (Lm' * x)) + fA_x_ominus;
end
