function Y = lr_deflation_scheme(A, x, m, k)
    % Step 1: Compute left and right eigenvectors
    [Lm, Rm, Dm] = compute_eigenvectors(A, m);

    % Step 2: Compute triangular matrix for critical eigenvalues
    f_Tm = compute_sign_function_diag(Dm);

    % Step 3: Compute x_ominus = (1 − Rm * Lm' ) * x 
    % ' -> transpose conjugate
    % .' -> transpose
    x_ominus = compute_x_ominus(Rm, Lm, x);

    % Step 4: Construct an orthonormal basis for the Krylov subspace using GMRES
    [~, ~, H, V] = gmres_method(A, x_ominus, x_ominus, k, 1e-8);
    % Task check the arnoldi relation is fulfilled.
    
    % Compute Hessenberg matrix
    %Hk = Vk' * A * Vk;
    Hk = H(1:size(A, 1), 1:size(A,2));
    Vk = V(1:size(A, 1), 1:size(A,2));

    % Step 5: Compute (the first column of) f(Hk) using Roberts’ iterative method
    f_Hk = compute_sign_function(Hk);

    % Step 6: Compute the approximation to f(A)x
    beta = norm(x_ominus);
    e1 = [1; zeros(size(A, 1) - 1, 1)];
    Y0 = x;
    Y = Rm * f_Tm * (conj(Lm') \ Y0) + beta .* Vk * f_Hk * e1;   
end
