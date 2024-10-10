function [fA_b, cost] = Combo_Lp_precond_quad_Impl_rest_arnoldi(A, b, k_values, m1, m2, max_iter, thick_num, tol, min_decay)
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
    %                  left prec. Restarted Implicit Arnoldi matrix approximation for all
    %                  values in k_values dimensions.
    %      cost      - No.of matrix{A} vector multiplications
    
    addpath(fullfile(pwd, 'Combo_Lp_precond_quad_Impl_rest_arnoldi'));

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
    kmax = max(k_values);
    cost = zeros(length(k_values), 1);

    H = [];
    V = zeros(length(v), kmax);
    V(:, 1) = v;

    [v, H, V,eta] = left_precondi_Arnoldi_process(A, m1, kmax, 1, V, H, theta);

    for l=1:no_k  %l-th column holds approx. for subspace dimension k(l)
        [fA_b(:,l), iter, ~] = Quad_based_impl_rest_arnoldi(A, v, V(:, 1:k_values(l)), H(1:k_values(l), 1:k_values(l)), eta, c_norm, theta, k_values(l), m1, m2, max_iter, thick_num, tol, min_decay);
        temp_m2 = m2 - 1;
        cost(l) = 1 + 2*(m1-1) + 2*(m1-1) + k_values(l)*(2 + 2*(m1-1) + 2*(m1-1));
        % disp([num2str(cost(l)), ', ', num2str(1)]);
        for i = 2:iter
            if temp_m2 ~= 0
                cost(l) = cost(l) + k_values(l)*(2 + 2*(m1-1) + 2*(m1-1));
                temp_m2 = temp_m2 - 1;
                % disp([num2str(cost(l)), ', ', num2str(i)]);
            else
                cost(l) = cost(l) + 2*k_values(l);
                % disp([num2str(cost(l)), ', ', num2str(i)]);
            end
        end
    end
end