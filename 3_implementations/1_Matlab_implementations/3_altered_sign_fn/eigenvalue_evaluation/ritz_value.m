function theta = ritz_value(A, m)
    %% Returns the sorted ritz values in vector form
    % Input: 
    %      A - N x N matrix
    %      m - no. of iterations for the krylov's subspace
    % Output: 
    %      theta - The sorted ritz values in vector form
    
    rng(2130); % setting random seed generator for reproducibility
    %% Arnoldi process to obtain ritz values
    
    N = size(A, 1);
    V = zeros(N, m+1); % m+1 arnoldi vectors
    
    % Generate basis Vm of Km(A, b)
    %% Correction 1
    % beta = norm(b);
    % V(:, 1) = b / beta;
    r0 = randn(N, 1);
    beta = norm(r0);
    V(:, 1) = r0 / beta;
    
    H = zeros(N, m+1);
    
    for j = 1:m
        % Apply matrix A to the last basis vector
        v = A * (V(:, j));
    
        % Arnoldi process: Orthogonalization
        for i = 1:j
            % Compute the coefficient for orthogonalization
            H(i, j) = V(:, i)' * v;
    
            % Subtract the projection of v onto the already computed basis vectors
            v = v - H(i, j) * V(:, i);
        end
    
        % Store the new basis vector
        H(j+1,j) = norm(v);
        V(:, j+1) = v / H(j+1, j);
    end
        
    % Verify arnoldi iteration
    % verify_Arnoldi_4_sqr(A, V, H);
    
    %% implementation of Ritz value based Newton iteration for A^{-1/2}
    %% Correction 2
    % theta = sort(eig(H(1:m,1:m))); % Ritz values <- eigen values of the Hessenberg matrix Hm.
    theta = leja_sort(eig(H(1:m,1:m)));
end