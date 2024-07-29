function fA_b = right_precondi_poly_fom(A, b, m, m1)
    % Input:
    %      A - N x N matrix
    %      b - N x 1 vector
    %      m - no. of iterations for the krylov's subspace
    %      m1 - No. of iterations for the krylov's subspace to be used in pre-conditioning polynomial
    % Output: 
    %      fA_b - f(A)b

    %% Finding the Preconditioned polynomial using Newton's interpolation using divided differences

    p_A = pre_condi_poly(A, b, m1);
    
    %% implementation of right polynomially preconditioned arnoldi process for A^{-1/2}
    
    N = size(A, 1);
    V = zeros(N, m+1); % m+1 arnoldi vectors
    e1 = zeros(m, 1);
    e1(1) = 1;
    b = A * b;
    
    % Generate basis Vm of Km(A, b)
    beta = norm(b);
    V(:, 1) = b / beta;
    
    H = zeros(N, m+1);
    
    for j = 1:m
        % Apply matrix A to the last basis vector
        y = p_A * V(:, j);
        u = p_A * y;
        w = A * (A * u);
    
        % Arnoldi process: Orthogonalization
        for i = 1:j
            % Compute the coefficient for orthogonalization
            H(i, j) = V(:, i)' * w;
    
            % Subtract the projection of v onto the already computed basis vectors
            w = w - H(i, j) * V(:, i);
        end
    
        % Store the new basis vector
        H(j+1,j) = norm(w);
        V(:, j+1) = w / H(j+1, j);
    end
    
    % Compute the approximation to f(A)x
    fA_b = p_A * (V(:, 1:m) * (inv(sqrtm(H(1:m, 1:m))) * e1 * beta));
end