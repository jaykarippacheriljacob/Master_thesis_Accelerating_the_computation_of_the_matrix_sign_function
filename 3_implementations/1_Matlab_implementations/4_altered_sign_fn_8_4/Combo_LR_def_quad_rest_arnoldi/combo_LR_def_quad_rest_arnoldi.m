function [fA_b, cost, restarts] = combo_LR_def_quad_rest_arnoldi(A, x, m_values, k_values, max_iter, tol, min_decay)
    %% Combination of LR-deflation and Quadrture based Explicit restarted arnoldi approximation for f(A)b.
    % Input:
    %      A         - n x n matrix.
    %      x         - n x 1 vector.
    %      m_values - vector of increasing no. of critical values for which defation has to be undergone.
    %      k_values  - vector of increasing Krylov subspace dimensions.
    %      max_iter  - Maximum no.of iterations for the restart of the Arnoldi decomposition.
    %      tol       - Set tolerance for stopping criteria.
    %      min_decay - the decay rate of error after each iteration.
    % Output:
    %      Y         - Approximation of f(A)b for k_values. We return the deflated 
    %                  left prec. Arnoldi matrix approximation for all
    %                  values in k_values dimensions.
    %      cost      - No.of matrix{A} vector multiplications.
    %      restarts  - No.of restarts completed undergone to converge
    
    addpath(fullfile(pwd, 'Combo_LR_def_quad_rest_arnoldi'));

    %% Initializing variables
    no_k = length(k_values);
    no_m = length(m_values);
    n = size(A,2);
    fA_x_ominus = zeros(n,no_k);
    cost = zeros(no_k*no_m, 1);
    restarts = zeros(no_k*no_m, 1);
    fA_b = zeros(n, no_k*no_m);
    j = 1;
    
    %% Step 1: Compute left and right eigenvectors
    mmax = max(m_values);
    [Rmmax, Lmmax, Dmmax] = compute_eigenvectors(A, mmax);

    for i = 1:length(m_values)
        m = m_values(i);
        Rm = Rmmax(:, 1:m);
        Lm = Lmmax(:, 1:m);
        Dm = Dmmax(1:m, 1:m);

        %% Step 2: Compute f(A) for critical eigenvalues
        f_Tm = compute_sign_function_diag(Dm);
    
        %% Step 3: Compute x_ominus = (1 − Rm * Lm' ) * x 
        % ' -> transpose conjugate
        % .' -> transpose
        x_ominus = compute_x_ominus(Rm, Lm, x);

        %% Step 4: Compute f(A)x_ominus using the restarted arnoldi process for k_values krylov subspace dimensions.
        for l=1:no_k  %l-th column holds approx. for subspace dimension k(l)
            [fA_x_ominus(:,l), restarts(j+l-1), ~,cost(j+l-1)] = Quad_based_restarted_arnoldi(A, x_ominus, k_values(l), max_iter, tol, min_decay);
        end

        %% Step 6: Compute the approximation to f(A)x
    
        fA_b(:,j:(j+no_k-1)) = Rm * (f_Tm * (Lm' * x)) + fA_x_ominus;
        j = j+no_k;
    end
    restarts = restarts - 1;
end
