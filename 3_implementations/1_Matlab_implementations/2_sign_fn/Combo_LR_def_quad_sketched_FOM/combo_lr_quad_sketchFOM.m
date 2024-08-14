function Y = combo_lr_quad_sketchFOM(A, x, m, k, s)
    % Input:
    %      A - n x n matrix
    %      x - n x 1 vector
    %      m - No. of critical values for which defation has to be undergone.
    %      k - no. of iterations for the krylov's subspace, m < min(N, s)
    %      s - sketch matrix row dimension
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

    fA_x_ominus = Quadrature_based_sketched_FOM(A, x_ominus, k, s);

    %% Step 6: Compute the approximation to f(A)x

    Y = Rm * (f_Tm * (Lm' * x)) + fA_x_ominus;
end
