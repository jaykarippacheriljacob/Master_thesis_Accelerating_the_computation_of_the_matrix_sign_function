function fm_k = Quadrature_based_restarted_arnoldi(A, b, m, max_iter)
    % Quadrature-based restarted Arnoldi approximation for f(A)b.
    % Input: 
    %      A - N x N matrix
    %      b - N x 1 vector
    %      m - no. of iterations for the kryl ov's subspace
    %      max_iter - Maximum no.of iterations for the restart of the Arnoldi decomposition
    % Output: 
    %      fm_k - f(A)b

    % Set tolerance
    tol = 1e-10;

    % Step 1: Compute arnoldi decomposition wrt A and b.
    [Hm, Vm] = Arnoldi_method(A, b, m);

    % Step 2: Set fm_1 := norm(b)*Vm_1*f(Hm_1)*e1
    b_norm = norm(b);
    e1 = zeros(m, 1);
    e1(1) = 1;
    fm_1 = b_norm * Vm(:, 1:m) * compute_sign_function(Hm(1:m, 1:m)) * e1;
    
    % Step 3: Restart cycles untill convergence
    for k = 2:max_iter
        
        % Step 4: Compute Arnoldi decomposition wrt A and vm+1_k-1
        [~, Vm] = Arnoldi_method(A, Vm(:, m+1), m);

        % Step 5: Compute h1m_k = em_k-1(Hm_k)e1 by quadrature of order l1
        %         and ompute h2m_k = em_k-1(Hm_k)e1 by quadrature of order
        %         l2.
        AVm = A * Vm(:, 1:m);
        [~, ~, ~, hm] = Quadrature_rule_invsqrt(Vm(:, 1:m), AVm, Vm(:, m+1), tol);

        % Step 6: Compute fm_k = fm_k-1 + norm(b) * Vm_k * hm_k
        fm_k = fm_1 + b_norm * Vm(:, 1:m) * hm;

        % Step 7: Check the relative error to check for the necesity of further
        %         restart.
        rel_err = norm(fm_1 - fm_k);
        if rel_err < tol
            % disp(rel_err);
            disp(['Number of restarts done: ', num2str(k - 1)]);
            break;
        else
            fm_1 = fm_k;
        end
    end
end