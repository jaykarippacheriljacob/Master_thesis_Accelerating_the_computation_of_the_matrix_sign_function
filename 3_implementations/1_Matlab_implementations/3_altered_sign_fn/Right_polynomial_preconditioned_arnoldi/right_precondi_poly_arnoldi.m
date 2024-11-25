function [fA_b, cost] = right_precondi_poly_arnoldi(A, b, k_values, k1)
    %% Right polynomial preconditioning approximation for f(A)b.
    % Input:
    %      A        - n x n matrix
    %      b        - n x 1 vector
    %      k_values - vector of increasing Krylov subspace dimensions. 
    %      k1       - No. of iterations for the krylov's subspace to be used in
    %                 pre-conditioning polynomial Arnoldi
    % Output: 
    %      fA_b     - Approximation of f(A)b for k_values. We return the 
    %                 right prec. Arnoldi matrix approximation for all
    %                 values in k_values dimensions.
    %      cost     - No.of matrix{A} vector multiplications

    %% Step 1: Approximate the matrix function on x_ominus using the preconditioned arnoldi process
    %          for the various values in k
    %          We do the Arnoldi process only once and then extract the Arnoldi approximations 
    %          for the different dimensions in k_values

    addpath(fullfile(pwd, 'Right_polynomial_preconditioned_FOM'));

    n = size(A,2);
    no_k = length(k_values);
    fA_b = zeros(n,no_k);
    kmax = max(k_values);
    cost = zeros(no_k, 1);

    [V,H,beta] = right_precondi_Arnoldi_process(A, b, kmax, k1);
    
    for l=1:no_k  %l-th column holds approx. for subspace dimension k(l)
        fA_b(:,l) = right_precondi_Arnoldi_approx(V,H,beta,k_values(l));
        cost(l) = 2*k1 + 1 + k_values(l)*(2 + 2*(k1-2) + 2*(k1-1));
    end
end