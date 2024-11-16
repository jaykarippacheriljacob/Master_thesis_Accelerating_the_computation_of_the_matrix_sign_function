function [f, iter, fm, cost] = Quad_based_imp_rest_arnoldi(A, b, m, max_iter, thick_num, tol, min_decay)
    %% Quadrature-based Implicit restarted Arnoldi approximation for f(A)b.
    % Input: 
    %      A         - N x N matrix
    %      b         - N x 1 vector
    %      m         - each restart cycle consists of m Arnoldi iterations
    %      max_iter  - Maximum no.of restart cycles
    %      thick_num - Number of target eigenvalues for implicit deflation
    %      tol       - Set tolerance for stopping criteria
    %      min_decay - the decay rate of error after each iteration.
    % Output: 
    %      f         - f(A)b
    %      iter      - No.of restarts done
    %      fm        - The f(A)*b calculated between each restarts
    %      cost      - No.of matrix{A} vector multiplications

    addpath(fullfile(pwd, 'Quad_based_Impl_restarted_arnoldi'));

    ell = 0;    % thick restart param
    subdiag = []; % subdiagonal entries of hessenberg matrices (for computing the norm of v_{m+1})
    f = zeros(size(b));
    b = A * b;
    b_norm = norm(b);
    fm = [];
    l = 8; % initial number of quadrature points

    % allocate memory [ >= m+max(ell) ]
    alloc = m + 20;
    V_big = zeros(length(b),alloc);
    v = b / b_norm;

    % restart loop starts here
    for k = 1:max_iter

        % compute A-invariant subspace of prev. cycle (thick restart)
        if k > 1
            ell_prev = ell;

            %% Rewrite it to a function as shown in the funm_quad.m
            % reordered Schur form and eigenvalues such that the ell "wanted" ones occur first
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
            %%
     
            thick_replaced{k-1} = D(1:ell);
            thick_interpol{k-1} = D(ell+1:end);
            if ell
                U = U(:,1:ell);
                V_big(:,1:ell) = V_big(:,1:m+ell_prev)*U;
                H_hat = U'*H*U;
                H = [H_hat; eta*U(end,:)];
            else
                H = [];
            end
        else
            H = [];
        end
     
        V_big(:,ell+1) = v;
     
        %% Rewrite the arnoldi process.
        % compute/extend Krylov decomposition
        % [ v,H,eta,breakdown, accuracy_flag ] = arnoldi( A,m+ell,H,ell+1,param );
        [v, H, V_big,eta] = Arnoldi_process(A, m+ell, ell+1, V_big, H);
        %%
    
        thick_interpol{k} = eig(H);
     
        % Schur form of H (for thick restarts)
        if isreal(H)
            [U,T] = schur(H,'real');
        else
            [U,T] = schur(H);
        end
        D = ordeig(T);
     
        % interpolation nodes currently active in f:
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
    
        %% Missing earlier
        if k>1
            active_nodes = [active_nodes; thick_replaced{k-1}];
        end
        %%
     
        % approximate the restart function:
        if k == 1 % in the first iteration g_1 = f and we use the "closed" form
            h2 = inv(sqrtm(H));
        else
            % Step : Compute h1m_k = em_k-1(Hm_k)e1 by quadrature of order l1
            %         and ompute h2m_k = em_k-1(Hm_k)e1 by quadrature of order
            %         l2
            [h2, l] = Quadrature_rule_invsqrt(A, active_nodes, subdiag, thick_replaced, H, m, tol, ell, k, l);
            fprintf('numbers of quadrature nodes = %d\n', l);
        end
     
        h_big = b_norm*h2(1:m+ell,1);
        if size(V_big,2) > length(h_big)
            h_big(size(V_big,2),1) = 0;
        end
     
        % update Krylov approximation
        f = V_big*h_big + f;
     
        fm = [ fm, f];
     
     
        % keep track of subdiagonal entries of Hessenberg matrix
        if m ~= 1
            subdiag = [ subdiag ; diag(H(end-m+1:end,end-m+1:end),-1) ; eta ];
        else
            subdiag = [ subdiag ; eta ];
        end
        %% Need to check if necessary or not.
        % s = eta;
        %%
        err_new = b_norm * h_big;
        % if (norm_update(k-1) / norm(f)) < tol
        if k > 2
            if err_new / err_old > min_decay
                disp(['Number of restarts done: ', num2str(k - 1)]);
                break;
            end
        end
        if (err_new / norm(fm)) < tol
            % disp(rel_err);
            disp(['Number of restarts done: ', num2str(k - 1)]);
            break;
        elseif k == max_iter
            disp(['Number of restarts done: ', num2str(k - 1)]);
        else
            err_old = err_new;
        end
    end
    iter = k;
    cost = 1 + 2 * m * iter;
end