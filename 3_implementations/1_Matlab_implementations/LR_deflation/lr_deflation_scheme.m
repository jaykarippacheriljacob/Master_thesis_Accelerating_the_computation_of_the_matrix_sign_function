function result = lr_deflation_scheme(A, x, m, k)
    % Step 1: Compute left and right eigenvectors
    [Lm, Rm, D] = compute_eigenvectors(A, m);

    % Step 2: Compute triangular matrix for critical eigenvalues
    f_Tm = compute_sign_function(D);

    % Step 3: Compute x_ominus = (1 − Rm * conj(Lm') ) * x
    x_ominus = compute_x_ominus(Rm, Lm, x);

    % Step 4: Construct an orthonormal basis for the Krylov subspace using GMRES
    [~, x, Vk] = gmres_method(A, x_ominus, x_ominus, k, 1e-8);
    
    % Compute Hessenberg matrix
    %Hk = Vk' * A * Vk;
    Hk = Vk(1:size(A, 1), 1:size(A,2));

    % Step 5: Compute (the first column of) f(Hk) using Roberts’ iterative method
    f_Hk_column = compute_sign_function(Hk);

    % Step 6: Compute the approximation to f(A)x
    result = Rm * f_Tm * (Lm' \ x) + x_ominus * f_Hk_column(1);
end
