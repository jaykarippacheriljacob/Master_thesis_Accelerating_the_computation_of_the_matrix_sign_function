function [v, H, V_big, h] = Arnoldi_proces(A, m, s, V_big, H)
    %% Computes the Krylov subspace basis vectors and its
    %  corresponding hessenberg matrix for a given dimension of subspace.
    % Inputs:
    %       A: The coefficient matrix (n x n)
    %       b: The right-hand side vector (n x 1)
    %       m: Maximum number of iterations
    %       s: Starting iteration number based on thick restart
    %       H: Upper Hessenberg matrix from Arnoldi process
    %       V: Arnoldi vector matrix
    % Outputs:
    %       H: Upper Hessenberg matrix from Arnoldi process
    %       V: Arnoldi vector matrix
    %       v: End Krylov subspace vector

    for j = s:m
        % Apply matrix A to the last basis vector
        v = A * (A * V_big(:, j));

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

    %% Verification on whether the arnoldi relationship is fulfilled or not
    % verify_Arnoldi_4_sqr(A, V_big(:,1:m+1), H);
end
