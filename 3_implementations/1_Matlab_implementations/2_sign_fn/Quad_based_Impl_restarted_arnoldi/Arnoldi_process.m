function [v, H, V_big, h] = Arnoldi_process(A, m, s, V_big, H)
    % Inputs:
    %       A: The coefficient matrix (n x n)
    %       b: The right-hand side vector (n x 1)
    %       m: Maximum number of iterations
    %       s: Starting iteration number based on thick restart
    % Outputs:
    %       H: Upper Hessenberg matrix from Arnoldi process
    %       V: Arnoldi vector matrix
    %       v: End Krylov subspace vector

    for j = s:m
        % Apply matrix A to the last basis vector
        v = A * V_big(:, j);

        % Arnoldi process: Orthogonalization
        for i = 1:j
            % Compute the coefficient for orthogonalization
            H(i, j) = V_big(:, i)' * v;
            % H(i, j) = v' * V(:, i); % check which direction is correct?

            % Subtract the projection of v onto the already computed basis vectors
            v = v - H(i, j) * V_big(:, i);
        end
        % Store the new basis vector
        H(j + 1, j) = norm(v);

        % add a check to see if there is a zero division or not.
        v = v / H(j + 1, j);
        if j < m
            V_big(:, j + 1) = v;
        end
        
    end
    h = H(m+1, m);
    H = H(1:m, 1:m);
end
