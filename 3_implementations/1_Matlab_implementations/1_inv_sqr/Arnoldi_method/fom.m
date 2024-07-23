function result = fom(A, b, m)
    % Input: 
    %      A - N x N matrix
    %      b - N x 1 vector
    %      m - no. of iterations for the krylov's subspace
    % Output: 
    %      result - f(A)b

    % Step 1: Initializing the necessary variables
    N = size(A, 1);
    V = zeros(N, m+1); % m+1 arnoldi vectors

    % Step 2: Generate basis Vm of Km(A, b)
    r0 = b; % Setting the residual to be the vector b
    beta = norm(r0);
    V(:, 1) = r0 / beta;

    H = zeros(N, m+1);

    for j = 1:m
        % Step 2.1: Apply matrix A to the last basis vector
        v = A * V(:, j);

        % Step 2.2: Arnoldi process: Orthogonalization
        for i = 1:j
            % Compute the coefficient for orthogonalization
            H(i, j) = V(:, i)' * v;

            % Subtract the projection of v onto the already computed basis vectors
            v = v - H(i, j) * V(:, i);
        end

        % Step 2.3: Store the new basis vector
        H(j+1,j) = norm(v);
        V(:, j+1) = v / H(j+1, j);
    end
    
    % Step 3: Verify arnoldi iteration
    verify_Arnoldi(A, V(:, 1:m), H(1:m, 1:m), H(m+1, m), V(:, m+1));
    V = V(:, 1:m);

    % Step 5: Compute thin QR decomposition Vm = Qm * Rm
    [Qm, Rm] = qr(V(:,1:j), 0);

    % Step 4: Compute the approximation to f(A)x
    f = compute_sign_function(Qm' * ((A * V(:, 1:m)) * inv(Rm)));
    result = V * (inv(Rm) * (f * (Qm' * b)));
end
