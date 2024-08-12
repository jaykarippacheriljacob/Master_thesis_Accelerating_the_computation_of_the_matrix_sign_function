function [H, V] = Arnoldi_method(A, b, k, trunc)
    % Inputs:
    %       A: The coefficient matrix (n x n)
    %       b: The right-hand side vector (n x 1)
    %       k: Maximum number of iterations
    %       trunc - Truncate orthogonalization to the last 'trunc' vector
    % Outputs:
    %       H: Upper Hessenberg matrix from Arnoldi process
    %       V: Arnoldi vector matrix
    
    n = size(b, 1);
    m = min(n, k + 1);
    V = zeros(n, m + 1);
    H = zeros(m + 1, m);

    % Initial residual vector
    beta = norm(b);
    V(:, 1) = b / beta;

    for j = 1:m
        % Apply matrix A to the last basis vector
        v = A * V(:, j);

        % Arnoldi process: Orthogonalization
        for i = max(1, j - trunc + 1):j
            % Compute the coefficient for orthogonalization
            H(i, j) = V(:, i)' * v;
            % H(i, j) = v' * V(:, i); % check which direction is correct?

            % Subtract the projection of v onto the already computed basis vectors
            v = v - H(i, j) * V(:, i);
        end
        % Store the new basis vector
        H(j + 1, j) = norm(v);

        % add a check to see if there is a zero division or not.
        V(:, j + 1) = v / H(j + 1, j);
        
    end
end
