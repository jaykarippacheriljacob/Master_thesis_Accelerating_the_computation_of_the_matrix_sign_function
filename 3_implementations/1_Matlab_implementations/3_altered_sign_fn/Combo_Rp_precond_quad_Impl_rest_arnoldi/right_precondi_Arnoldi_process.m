function [w, H, V_big, h] = right_precondi_Arnoldi_process(A, m1, m, s, V_big, H, theta)
    % Input: 
    %      p_A  - preconditioned polynomial based on Newton interpolation of
    %             divided differences.
    %      A    - N x N matrix
    %      b    - N x 1 vector
    %      m    - dimension of the (preconditioned) krylov subspace
    %      m1   - no of interpolaiton points for the prec. polynomial 
    %             which has degree m1-1 
    % Output: 
    %      V    - matrix with m+1 columns holding the Arnoldi vectors for the
    %             Krylov subspace with matrix A^2p^2(A^2) and vector p^(A)*A*b
    %      H    - matrix of orth. coefficients
    %      beta - norm of the vector p^(A)*A*b

    %% right polynomially preconditioned Arnoldi process for A^2
    for j = s:m
        % Apply A^2*p^2(A^2) to the last basis vector v = V(:,j)
        y = eval_pre_condi_poly(A, V_big(:, j), theta, m1);  %p(A^2)*v
        u = eval_pre_condi_poly(A, y, theta, m1);        %p^2(A^2)*v 
        w = A * (A * u);                                 % A^2*p^2(A^2)*v
    
        % Arnoldi process: Orthogonalization
        for i = 1:j
            % Compute the coefficient for orthogonalization
            H(i, j) = V_big(:, i)' * w;
    
            % Subtract the projection of v onto the already computed basis vectors
            w = w - H(i, j) * V_big(:, i);
        end
    
        % Store the new basis vector
        H(j+1,j) = norm(w);

        % add a check to see if there is a zero division or not.
        w = w / H(j + 1, j);
        if j < m
            V_big(:, j + 1) = w;
        end

    end
    h = H(m+1, m);
    H = H(1:m, 1:m);
end