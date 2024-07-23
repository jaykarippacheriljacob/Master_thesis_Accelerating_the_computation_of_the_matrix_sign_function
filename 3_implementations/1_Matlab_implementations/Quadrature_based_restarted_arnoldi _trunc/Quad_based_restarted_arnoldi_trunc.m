function [fm_k, iter, f] = Quad_based_restarted_arnoldi_trunc(A, b, m, max_iter, tol, min_decay, trunc)
    % Quadrature-based restarted Arnoldi approximation for f(A)b.
    % Input: 
    %      A - N x N matrix
    %      b - N x 1 vector
    %      m - no. of iterations for the kryl ov's subspace
    %      max_iter - Maximum no.of iterations for the restart of the Arnoldi decomposition
    %      tol - Set tolerance for stopping criteria
    %      min_decay - the decay rate of error after each iteration.
    %      trunc - Truncate orthogonalization to the last 'trunc' vector
    % Output: 
    %      fm_k - f(A)b
    %      No.of restarts done
    %      f - The f(A)*b calculated between each restarts

    subdiag = []; % subdiagonal entries of hessenberg matrices (for computing
                  % the norm of v_{m+1})
    f = [];
    active_nodes = []; % Interpolation nodes (Ritz values) from each
                         % restart cycle

    % Step 1: Compute arnoldi decomposition wrt A and b.
    [Hm, Vm] = Arnoldi_method(A, b, m, trunc);
    active_nodes = [active_nodes; sort(eig(Hm(1:m, 1:m)))];
    subdiag = [subdiag; diag(Hm(1:m, 1:m), -1); Hm(m+1, m)];

    % Step 2: Set fm_1 := norm(b)*Vm_1*f(Hm_1)*e1
    b_norm = norm(b);
    e1 = zeros(m, 1);
    e1(1) = 1;
    f_Hm = inv(sqrtm(Hm(1:m, 1:m)));
    fm_1 = b_norm * Vm(:, 1:m) * f_Hm * e1;
    f = [f, fm_1];
    
    % Step 3: Restart cycles untill convergence
    for k = 2:max_iter
        
        % Step 4: Compute Arnoldi decomposition wrt A and vm+1_k-1
        [Hm, Vm] = Arnoldi_method(A, Vm(:, m+1), m, trunc);

        % Step 5: Compute h1m_k = em_k-1(Hm_k)e1 by quadrature of order l1
        %         and ompute h2m_k = em_k-1(Hm_k)e1 by quadrature of order
        %         l2
        subdiag = [subdiag; diag(Hm(1:m, 1:m), -1); Hm(m+1, m)];
        h2 = Quadrature_rule_invsqrt(A, active_nodes, subdiag, Hm(1:m, 1:m), tol);
        active_nodes = [active_nodes; sort(eig(Hm(1:m, 1:m)))];

        % Step 6: Compute fm_k = fm_k-1 + norm(b) * Vm_k * hm_k
        h_new = b_norm * h2;
        fm_k = fm_1 + Vm(:, 1:m) * h_new;

        % Step 7: Check the stopping criteria necessary for further restart.
        %         { decay for norm of updates are checked to see if it's
        %           below the stopping accuracy }
        % norm_update = [norm_update, b_norm * h_new];
        err_new = b_norm * h_new;
        f = [ f, fm_k];
        % if (norm_update(k-1) / norm(f)) < tol
        if k ~= 2
            if err_new / err_old > min_decay
                disp(['Number of restarts done: ', num2str(k - 1)]);
                break;
            end
        end
        if (err_new / norm(f)) < tol
            % disp(rel_err);
            disp(['Number of restarts done: ', num2str(k - 1)]);
            break;
        elseif k == max_iter
            disp(['Number of restarts done: ', num2str(k - 1)]);
        else
            fm_1 = fm_k;
            err_old = err_new;
        end
    end
    iter = k;
end