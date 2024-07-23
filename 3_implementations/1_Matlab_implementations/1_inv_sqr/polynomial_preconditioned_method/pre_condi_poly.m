function p_A = pre_condi_poly(A, b, m)
    % Input: 
    %      A - N x N matrix
    %      m - no. of iterations for the krylov's subspace
    %      b - N x 1 vector
    % Output: 
    %      p_A - preconditioned polynomial based on Newton interpolation of
    %            divided differences.
    
    %% Arnoldi process to obtain ritz values
    
    N = size(A, 1);
    V = zeros(N, m+1); % m+1 arnoldi vectors
    
    % Generate basis Vm of Km(A, b)
    beta = norm(b);
    V(:, 1) = b / beta;
    
    H = zeros(N, m+1);
    
    for j = 1:m
        % Apply matrix A to the last basis vector
        v = A * V(:, j);
    
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
    verify_Arnoldi(A, V(:, 1:m), H(1:m, 1:m), H(m+1, m), V(:, m+1));
    
    %% implementation of Ritz value based Newton iteration for A^{-1/2}

    theta = sort(eig(H(1:m,1:m))); % Ritz values <- eigen values of the Hessenberg matrix Hm.
    p_A = zeros(N); % Polynomial for preconditioning of matrix A
    I = eye(N); % Identity matrix of size N
    f = @(t) t.^(-1/2); % Define the function
    dd = zeros(m, m); % Initialize the divided difference table
    
    % Compute the first column of divided differences (function values at points theta)
    for i = 1:m
        dd(i, 1) = f(theta(i));
    end
    
    % Compute the divided differences
    for j = 2:m
        for i = 1:(m-j+1)
            dd(i, j) = (dd(i+1, j-1) - dd(i, j-1)) / (theta(i+j-1) - theta(i));
        end
    end
    
    % Extract the coefficients
    coefficients = dd(1, :);
    
    % Expressing the polynomial using newton basis
    for i=1:m
        if i > 1
                p_term = p_term * (A - theta(i-1) * I);
        else
            p_term = I;
        end
        p_A = p_A + coefficients(i) * p_term;
    end
end