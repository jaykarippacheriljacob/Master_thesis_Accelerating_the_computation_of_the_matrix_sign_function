function [fm_k, iter, f] = Quadrature_based_restarted_arnoldi(A, b, m, max_iter, thick_num, tol, min_decay)
    % Quadrature-based restarted Arnoldi approximation for f(A)b.
    % Input: 
    %      A - N x N matrix
    %      b - N x 1 vector
    %      m - no. of iterations for the kryl ov's subspace
    %      max_iter - Maximum no.of iterations for the restart of the Arnoldi decomposition
    %      thick_num - Number of target eigenvalues for implicit deflation
    %      tol - Set tolerance for stopping criteria
    %      min_decay - the decay rate of error after each iteration.
    % Output: 
    %      fm_k - f(A)b
    %      No.of restarts done
    %      f - The f(A)*b calculated between each restarts

    subdiag = []; % subdiagonal entries of hessenberg matrices (for computing
                  % the norm of v_{m+1})
    f = [];

    ell = 0; % thick restart param

    % Vm = zeros(length(b), m);
    % Hm = zeros(m, m);
    Vm = [];
    Hm = [];

    % Step 1: Compute arnoldi decomposition wrt A and b.
    [Hm, Vm, h, v] = Arnoldi_method(A, b, m, Hm, Vm, ell);
    subdiag = [subdiag; diag(Hm(end-m+1:end, end-m+1:end), -1); h];

    % Step 2: Set fm_1 := norm(b)*Vm_1*f(Hm_1)*e1
    b_norm = norm(b);
    e1 = zeros(m, 1);
    e1(1) = 1;
    f_Hm = inv(sqrtm(Hm));
    fm_1 = b_norm * Vm * f_Hm * e1;
    f = [f, fm_1];
    l = 8; % numbers of quadrature nodes

    % Step 3: Setting up thick interpolation for 1^{st} restart
    thick_interpol{1} = eig(Hm);

    % Step 4: Computation of the schur decomposition of the Hessenberg
    %         matrix
    [U, T] = schur(Hm, 'real');
    % D = ordeig(T);
    
    %% Need to be completely reordered
    % Step 5: Restart cycles untill convergence
    for k = 2:max_iter
        
        % Step : Reordering of Schur form and eigenvalues such that the ell
        %        "wanted" ones occur first
        ell_prev = ell;
        ell = thick_num;
        [~, jj] = sort(abs(real(D)), 1, 'ascend');
        ind = 0 * D;
        ind(jj(1:ell)) = 1;
        [U, T] = ordschur(U, T, ind);
        D = ordeig(T);
        % do not break conjugate eigenvalues in the real Schur form
        if(length(D) > ell && abs(conj(D(ell+1)) - D(ell)) < 100 * eps * abs(D(ell)))
            ell = ell + 1;
        end

        % Step : Compute A-invariant subspace of prev. cycle (thick restart)
        thick_replaced{k-1} = D(1:ell);
        thick_interpol{k-1} = D(ell+1:end);
        if ell
            U = U(:,1:ell);
            Vm(:,1:ell) = Vm(:,1:m+ell_prev) * U;
            H_hat = U' * Hm * U;
            Hm = [H_hat; h * U(end,:)];
        % else
            % Hm = [];
        end
        Vm(:, ell+1) = v;

        % Step : Compute Arnoldi decomposition wrt A and vm+1_k-1
        [Hm, Vm, h, v] = Arnoldi_method(A, v, m, Hm, Vm, ell);

        % Step : Setting up thick interpolation for k^{th} restart
        thick_interpol{k} = eig(Hm);

        % Step : Computation of the schur decomposition of the Hessenberg
        %         matrix
        [U, T] = schur(Hm, 'real');
        D = ordeig(T);

        % Step : Interpolation nodes currently active in f:
        %    k = 1 -> []
        %    k = 2 -> [ out1. {1} ; out1.thick_replaced{1} ]
        %    k = 3 -> [ out1.thick_interpol{1} ; out1.thick_interpol{2} ; out1.thick_replaced{2} ]
        %    ...
        % (note that the "replaced nodes" are going to be replaced in the next
        % sweep, respectively)
        active_nodes = [];
        for kk = 1:k-1
            active_nodes = [ active_nodes ; thick_interpol{kk} ];
        end
        active_nodes = [ active_nodes ; thick_replaced{k-1} ];
        

        % Step : Compute h1m_k = em_k-1(Hm_k)e1 by quadrature of order l1
        %         and ompute h2m_k = em_k-1(Hm_k)e1 by quadrature of order
        %         l2
        [h2, l] = Quadrature_rule_invsqrt(A, active_nodes, subdiag, thick_replaced, Hm, m, tol, ell, k, l);
    
        % Step : keep track of subdiagonal entries of Hessenberg matrix
        subdiag = [subdiag; diag(Hm(end-m+1:end, end-m+1:end), -1); h];

        % Step : Compute fm_k = fm_k-1 + norm(b) * Vm_k * hm_k
        h_new = b_norm * h2(1:m+ell, 1);
        if(size(Vm, 2) > length(h_new))
            h_new(size(Vm, 2), 1) = 0;
        end
        fm_k = fm_1 + Vm * h_new;

        % Step : Check the stopping criteria necessary for further restart.
        %         { decay for norm of updates are checked to see if it's
        %           below the stopping accuracy };
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