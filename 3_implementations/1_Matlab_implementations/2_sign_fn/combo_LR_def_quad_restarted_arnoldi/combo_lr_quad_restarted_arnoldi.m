function Y = combo_lr_quad_restarted_arnoldi(A, x, m, k, max_iter, tol, min_decay)
    % Input:
    %      A - n x n matrix
    %      x - n x 1 vector
    %      m - No. of critical values for which defation has to be undergone.
    %      k - no. of iterations for the kryl ov's subspace
    %      max_iter - Maximum no.of iterations for the restart of the Arnoldi decomposition
    %      tol - Set tolerance for stopping criteria
    %      min_decay - the decay rate of error after each iteration.

    addpath("C:\Users\jkjbt\Documents\GitHub\Master_thesis_Accelerating_the_computation_of_the_matrix_sign_function\3_implementations\1_Matlab_implementations\2_sign_fn\combo_LR_def_quad_restarted_arnoldi");

    %% Step 1: Compute left and right eigenvectors

    [Rm, Lm, Dm] = compute_eigenvectors(A, m);

    %% Step 2: Compute triangular matrix for critical eigenvalues

    f_Tm = compute_sign_function_diag(Dm);

    %% Step 3: Compute x_ominus = (1 − Rm * Lm' ) * x 

    % ' -> transpose conjugate
    % .' -> transpose
    x_ominus = compute_x_ominus(Rm, Lm, x);

    %% Step 4: Construct an orthonormal basis for the Krylov subspace using the arnoldi process

    fA_x_ominus = Quadrature_based_restarted_arnoldi(A, x_ominus, k, max_iter, tol, min_decay);

    %% Step 6: Compute the approximation to f(A)x

    Y = Rm * (f_Tm * (Lm' * x)) + fA_x_ominus;
end
