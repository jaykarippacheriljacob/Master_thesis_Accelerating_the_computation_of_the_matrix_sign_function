function [x, resvec,H] = gmres_method(A, b, x0, k, tol)
    n = length(b);
    m = min(n, k + 1);
    Q = zeros(n, m + 1);
    H = zeros(m + 1, m);

    % Initialize residual vector
    r0 = b - A * x0;
    beta = norm(r0);
    Q(:, 1) = r0 / beta;

    resvec = zeros(m, 1);

    for j = 1:m
        % Apply matrix A to the last basis vector
        v = A * Q(:, j);

        % Arnoldi process: Orthogonalization
        for i = 1:j
            % Compute the coefficient for orthogonalization
            H(i, j) = Q(:, i)' * v;

            % Subtract the projection of v onto the already computed basis vectors
            v = v - H(i, j) * Q(:, i);
        end

        % Store the new basis vector
        H(j + 1, j) = norm(v);
        Q(:, j + 1) = v / H(j + 1, j);

        % Construct the (j+1)-dimensional Krylov subspace matrix
        G = H(1:j+1, 1:j) \ [beta; zeros(j, 1)];

        % Compute the solution minimizing the residual norm within the Krylov subspace
        x = x0 + Q(:, 1:j) * G;

        % Compute the new residual norm
        resvec(j) = norm(b - A * x);

        % Check convergence
        if resvec(j) < tol
            break;
        end
    end
end
