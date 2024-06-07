function Y = lr_deflation_scheme(A, x, m, k)
    % Step 1: Compute left and right eigenvectors
    [Rm, Lm, Dm] = compute_eigenvectors(A, m);

    % Step 2: Compute triangular matrix for critical eigenvalues
    f_Tm = compute_sign_function_diag(Dm);

    % Step 3: Compute x_ominus = (1 − Rm * Lm' ) * x 
    % ' -> transpose conjugate
    % .' -> transpose
    x_ominus = compute_x_ominus(Rm, Lm, x);

    % Step 4: Construct an orthonormal basis for the Krylov subspace using
    % the arnoldi process
    [~, H, V] = Arnoldi_method(A, x_ominus, zeros(size(x_ominus, 1),1), k);
    
    % Compute Hessenberg matrix
    %Hk = Vk' * A * Vk;
    k = size(V, 2) - 1;
    Hk = H(1:k, 1:k);
    Vk = V(1:size(A, 1), 1:k);

    % Task check the arnoldi relation is fulfilled.
    % A*Vm - Vm*Hm ~= h(m+1,m)*q(m+1)*e(m).' *******
    verify_Arnoldi(A, Vk, Hk, H(k+1, k), V(:, k+1));
    
    % Step 5: Compute (the first column of) f(Hk) using Roberts’ iterative method
    f_Hk = compute_sign_function(Hk);
    %disp(f_Hk);

    % Step 6: Compute the approximation to f(A)x
    beta = norm(x_ominus);
    e1 = [1; zeros(k-1, 1)];
    Y = Rm * (f_Tm * (Lm' * x)) + beta * (Vk * (f_Hk * e1));
end
