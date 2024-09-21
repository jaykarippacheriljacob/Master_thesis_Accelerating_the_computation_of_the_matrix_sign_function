function [fA_b] = combo_LR_def_quad_sketched_FOM(A, x, m, k_values, s, tol)
    % Input:
    %      A        - n x n matrix
    %      x        - n x 1 vector
    %      m        - No. of critical values for which defation has to be undergone.
    %      k_values - vector of increasing subspace dimensions. We return the deflated 
    %                 left prec. Arnoldi matrix approximation for all values in m 
    %      s        - sketch matrix row dimension
    %      tol      - parameter for setting the tolerance of the quadrature
    %                 rule
    % Output: 
    %      fA_b     - Approximation of f(A)b for k_values

    addpath(fullfile(pwd, 'Combo_LR_def_quad_sketched_FOM'));
    
    %% Step 1: Compute left and right eigenvectors
    [Rm, Lm, Dm] = compute_eigenvectors(A, m);

    %% Step 2: Compute f(A) for critical eigenvalues
    f_Tm = compute_sign_function_diag(Dm);

    %% Step 3: Compute x_ominus = (1 − Rm * Lm' ) * x 
    % ' -> transpose conjugate
    % .' -> transpose
    x_ominus = compute_x_ominus(Rm, Lm, x);

    %% Step 4: Construct an orthonormal basis for the Krylov subspace using the arnoldi process
    %          We do the Arnoldi process only once and then extract the Arnoldi approximations 
    %          for the different dimensions in k_values 
    
    n = size(A,2);
    no_k = length(k_values);
    fA_x_ominus = zeros(n,no_k);
    kmax = max(k_values);

    % Setting the tolerance of the quadrature rule if not provided.
    if nargin < 6
        tol = 1e-10;
    end

    [H, V, ~, SV, SAV, Sb] = QS_FOM_process(A, x_ominus, kmax, s);
    
    %% Step 5: Compute f(A)x_ominus
    for l=1:no_k  %l-th column holds approx. for subspace dimension k(l)
        fA_x_ominus(:,l) = QS_FOM_approx(V, H, SV, SAV, Sb, k_values(l), tol);
    end
    
    %% Step 6: Compute the approximation to f(A)x
    fA_b = Rm * (f_Tm * (Lm' * x))*ones(1,no_k) + fA_x_ominus;
end
