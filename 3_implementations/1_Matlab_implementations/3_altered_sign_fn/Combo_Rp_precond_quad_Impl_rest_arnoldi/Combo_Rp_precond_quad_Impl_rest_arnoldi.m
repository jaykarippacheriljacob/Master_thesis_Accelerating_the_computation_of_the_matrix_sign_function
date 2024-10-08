function [fA_b] = Combo_Rp_precond_quad_Impl_rest_arnoldi(A, b, k_values, m1, m2, max_iter, thick_num, tol, min_decay)
    % Quadrature-based restarted Arnoldi approximation for f(A)b.
    % Input: 
    %      A         - N x N matrix
    %      b         - N x 1 vector
    %      k_values  - each restart cycle consists of m Arnoldi iterations
    %      m1        - no of interpolaiton points for the prec. polynomial 
    %                  which has degree m1-1
    %      m2        - No. of times the preconditioned Arnoldi process has to be exceuted.
    %      max_iter  - Maximum no.of restart cycles
    %      thick_num - Number of target eigenvalues for implicit deflation
    %      tol       - Set tolerance for stopping criteria
    %      min_decay - the decay rate of error after each iteration.
    % Output: 
    %      fA_b      - Approximation of f(A)b for k_values. We return the 
    %                  lright prec. Restarted Implicit Arnoldi matrix approximation for all
    %                  values in k_values dimensions.
    
    addpath(fullfile(pwd, 'Combo_Rp_precond_quad_Impl_rest_arnoldi'));

    %% Defining initial parameters
    b = A * b;

    %% Step 1: determine interpolation nodes for the polynomials
    theta = ritz_value(A, m1);
    
    %% Step 2: Generate basis Vm of Km(A, b)
    c = eval_pre_condi_poly(A, b, theta, m1);
    c_norm = norm(c);
    v = c / c_norm;

    %% Step 3: Compute f(A)b using the restarted arnoldi process for k_values krylov subspace dimensions.
    no_k = length(k_values);
    n = size(A,2);
    fA_b = zeros(n,no_k);
    for l=1:no_k  %l-th column holds approx. for subspace dimension k(l)
        [fA_b(:,l), ~, ~] = Quad_based_impr_rest_arnoldi(A, v, c_norm, theta, k_values(l), m1, m2, max_iter, thick_num, tol, min_decay);
    end
end