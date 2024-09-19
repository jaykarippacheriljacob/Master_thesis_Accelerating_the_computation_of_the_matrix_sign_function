function fA_b = right_precondi_poly_fom(A, b, k_values, k1)
    % Input: 
    %      p_A      - preconditioned polynomial based on Newton interpolation of
    %                 divided differences.
    %      A        - N x N matrix
    %      b        - N x 1 vector
    %      k_values - vector of increasing subspace dimensions. We return the deflated 
    %                 left prec. Arnoldi matrix approximation for all values in m 
    %      k1       - No. of iterations for the krylov's subspace to be used in
    %                 pre-conditioning polynomial Arnoldi
    % Output: 
    %      fA_b     - f(A)b

    %% Step 1: Approximate the matrix function on x_ominus using the preconditioned arnoldi process
    %          for the various values in k
    %          We do the Arnoldi process only once and then extract the Arnoldi approximations 
    %          for the different dimensions in k_values

    addpath(fullfile(pwd, 'right_polynomial_preconditioned_FOM'));

    n = size(A,2);
    no_k = length(k_values);
    fA_b = zeros(n,no_k);
    kmax = max(k_values);

    [V,H,beta] = right_precondi_Arnoldi_process(A, b, kmax, k1);
    
    for l=1:no_k  %l-th column holds approx. for subspace dimension k(l)
        fA_b(:,l) = right_precondi_Arnoldi_approx(V,H,beta,k_values(l));
    end
end