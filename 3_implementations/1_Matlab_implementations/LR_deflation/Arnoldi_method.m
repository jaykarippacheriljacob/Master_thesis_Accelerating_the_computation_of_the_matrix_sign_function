function [resvec, H, V] = Arnoldi_method(A, b, x0, k)
    % GMRES_METHOD Solves a linear system using the GMRES method.
    % Inputs:
    %       A: The coefficient matrix (n x n)
    %       b: The right-hand side vector (n x 1)
    %       x0: The initial guess for the solution (n x 1)
    %       k: Maximum number of iterations
    %       tol: Tolerance for convergence xxxxxxxx
    % Outputs:
    %       x: The approximate solution vector (n x 1) xxxxxxxx
    %       resvec: Residual norms at each iteration xxxxxxxxx
    %       H: Upper Hessenberg matrix from Arnoldi process
    %       V: Arnoldi vector matrix
    
    n = size(b, 1);
    m = min(n, k + 1);
    V = zeros(n, m + 1);
    H = zeros(m + 1, m);

    % Initialize residual vector
    % the initial vector should always be zero, pass zero vector for x0
    r0 = b - A * x0;
    beta = norm(r0);
    V(:, 1) = r0 / beta;

    resvec = zeros(m, 1);

    for j = 1:m
        % Apply matrix A to the last basis vector
        v = A * V(:, j);

        % Arnoldi process: Orthogonalization
        for i = 1:j
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
