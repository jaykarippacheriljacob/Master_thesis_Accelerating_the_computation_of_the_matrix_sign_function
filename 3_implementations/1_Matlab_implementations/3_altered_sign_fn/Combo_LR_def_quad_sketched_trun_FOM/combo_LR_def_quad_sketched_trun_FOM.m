function [fA_b] = combo_LR_def_quad_sketched_trun_FOM(A, x, m, k_values, s, trunc, tol)
    % Input:
    %      A        - n x n matrix
    %      x        - n x 1 vector
    %      m        - No. of critical values for which defation has to be undergone.
    %      k_values - vector of increasing subspace dimensions. We return the deflated 
    %                 left prec. Arnoldi matrix approximation for all values in m 
    %      s        - sketch matrix row dimension
    %      trunc    - Truncate orthogonalization to the last 'trunc' vector
    %      tol      - parameter for setting the tolerance of the quadrature
    %                 rule
    % Output: 
    %      fA_b     - Approximation of f(A)b for k_values

    addpath(fullfile(pwd, 'Combo_LR_def_quad_sketched_trun_FOM'));
    
    %% Step 1: Compute left and right eigenvectors
    [Rm, Lm, Dm] = compute_eigenvectors(A, m);

    %% Step 2: Compute f(A) for critical eigenvalues
    f_Tm = compute_sign_function_diag(Dm);

    %% Step 3: Compute x_ominus = (1 − Rm * Lm' ) * x 
    % ' -> transpose conjugate
    % .' -> transpose
    x_ominus = compute_x_ominus(Rm, Lm, x);

    %% Step 4: Construct a sketched basis for the Krylov subspace using the arnoldi process
    %          We do the sketeched Arnoldi process only once and then extract the sketched Arnoldi approximations 
    %          for the different dimensions in k_values 
    
    n = size(A,2);
    no_k = length(k_values);
    fA_x_ominus = zeros(n,no_k);
    kmax = max(k_values);

    % Setting the tolerance of the quadrature rule if not provided.
    if nargin < 7
        tol = 1e-10;
    end

    [H, V, ~, SV, SAV, Sb] = QS_trun_FOM_process(A, x_ominus, kmax, s, trunc);
    % Checking whether the arnoldi relation is fulfilled.
    % A*Vm - Vm*Hm ~= h(m+1,m)*q(m+1)*e(m).' *******
    verify_Arnoldi_4_sqr(A, V, H);
    
    %% Step 5: Compute f(A)x_ominus
    for l=1:no_k  %l-th column holds approx. for subspace dimension k(l)
        fA_x_ominus(:,l) = QS_trun_FOM_approx(V, H, SV, SAV, Sb, k_values(l), tol);
    end
    
    %% Step 6: Compute the approximation to f(A)x
    fA_b = Rm * (f_Tm * (Lm' * x))*ones(1,no_k) + fA_x_ominus;
end
