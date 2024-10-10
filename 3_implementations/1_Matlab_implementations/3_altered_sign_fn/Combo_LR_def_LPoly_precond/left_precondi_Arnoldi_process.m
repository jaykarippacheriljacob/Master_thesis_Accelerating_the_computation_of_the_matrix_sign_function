function [V ,H, beta] = left_precondi_Arnoldi_process(A, b, m, m1)
    %% Computes the Left preconditioned Krylov subspace basis vectors and its
    %  corresponding hessenberg matrix for a given dimension of subspace.
    % Input: 
    %      p_A - preconditioned polynomial based on Newton interpolation of
    %            divided differences.
    %      A - N x N matrix.
    %      b - N x 1 vector.
    %      m - dimension of the (preconditioned) krylov subspace.
    %      m1 - no of interpolaiton points for the prec. polynomial 
    %           which has degree m1-1 .
    % Output: 
    %      V - matrix with m+1 columns holding the Arnoldi vectors for the
    %      Krylov subspace with matrix A^2p^2(A^2) and vector p^(A)*A*b.
    %      H - magtrix of orth. coefficients.
    %      beta: norm of the vector p^(A)*A*b.

    %% Initializing variables
    N = size(A, 1);
    V = zeros(N, m+1); % m+1 arnoldi vectors 
    H = zeros(N, m+1);
    
    %% determine interpolation nodes for the polynomials
    theta = ritz_value(A, m1);

    %% left polynomially preconditioned Arnoldi process for A^2    
    % Generate basis Vm of Km(A, b)
    c = eval_pre_condi_poly(A, A*b, theta, m1);
    beta = norm(c);
    V(:, 1) = c / beta;
    
    for j = 1:m
        % Apply matrix A to the last basis vector
        u = A * (A * V(:, j));
        y = eval_pre_condi_poly(A, u, theta, m1);
        w = eval_pre_condi_poly(A, y, theta, m1);
    
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
    
end