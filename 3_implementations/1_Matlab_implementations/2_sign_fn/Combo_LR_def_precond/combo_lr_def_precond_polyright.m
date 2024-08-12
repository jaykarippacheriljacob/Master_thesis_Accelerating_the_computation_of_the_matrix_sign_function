function Y = combo_lr_def_precond_polyright(A, x, m, k, k1)
    % Input:
    %      A - n x n matrix
    %      x - n x 1 vector
    %      m - No. of critical values for which defation has to be undergone.
    %      k - no. of iterations for the krylov's subspace, m < n
    %      k1 - No. of iterations for the krylov's subspace to be used in
    %           pre-conditioning polynomial Arnoldi
    %% Add the path that contains the function for better operation of the function
    addpath("C:\Users\jkjbt\Documents\GitHub\Master_thesis_Accelerating_the_computation_of_the_matrix_sign_function\3_implementations\1_Matlab_implementations\2_sign_fn\Combo_LR_def_precond");

    %% Step 1: Compute left and right eigenvectors

    [Rm, Lm, Dm] = compute_eigenvectors(A, m);

    %% Step 2: Compute f(A) for critical eigenvalues

    f_Tm = compute_sign_function_diag(Dm);

    %% Step 3: Compute x_ominus = (1 − Rm * Lm' ) * x 

    % ' -> transpose conjugate
    % .' -> transpose
    x_ominus = compute_x_ominus(Rm, Lm, x);

    %% Step 4: Construct an orthonormal basis for the Krylov subspace using the arnoldi process

    fA_x_ominus = right_precondi_poly_fom(A, x_ominus, k, k1);

    %% Step 6: Compute the approximation to f(A)x

    Y = Rm * (f_Tm * (Lm' * x)) + fA_x_ominus;
end
