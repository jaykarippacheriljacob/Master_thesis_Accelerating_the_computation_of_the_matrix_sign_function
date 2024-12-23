function result = Quadrature_based_sketched_FOM(A, b, m, s)
    % Input: 
    %      A - N x N matrix
    %      b - N x 1 vector
    %      m - no. of iterations for the krylov's subspace, m < min(N, s)
    %      s - sketch matrix row dimension
    % Output: 
    %      result - f(A)b

    %% Step 1: Draw sketching matrix S € C(s*N)
    N = size(b,1);
    D = speye(N);
    D = D(randperm(N, s), :);
    e = rand(N, 1);
    e = (e >= 0.5) - (e < 0.5);
    E = spdiags(e, 0, N, N);
    S = @(M) sqrt(N/s)*(D*(dct(E*M)));

    % check if m < min(N, s)
    if m > min(N, s)
        disp('Error: "The condition m < min(N, s) is not satisfied."');
    end

    V = zeros(N, m+1);
    SV = zeros(s, m);
    SAV = zeros(s, m);
    
    r0 = A * b;
    % Sb = S(A * b);
    Sb = S(r0);

    tol = 1e-10;

    % Step 2: Generate basis Vm of Km(A, b), as well as SVm
    % and SAVm.
    % r0 = A * b;
    beta = norm(r0);
    V(:, 1) = r0 / beta;

    H = zeros(N, m+1);

    for j = 1:m
        % Apply matrix A to the last basis vector
        v = A * (A * V(:, j));

        SV(:, j) = S(V(:, j));
        SAV(:,j) = S(v);

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
    V = V(:, 1:m);

    % Step 3: Compute thin QR decomposition SVm = Qm * Rm
    [Qm, Rm] = qr(SV(:,1:j), 0);

    % Step 4: Compute SVm = Qm, SAVm = (SAVm)*inv(Rm), Vm = Vm*inv(Rm)
    SV = Qm;
    SAV = SAV(:, 1:j) / Rm;
    V = V / Rm;

    % Step 5: Compute the quadrature rules ql1(S, A, Vm, b) and ql2(S, A, Vm, b)
    [~, ~, ~, h] = Quad_rule_invsqrt(SV, SAV, Sb, tol);

    % Step 6: Compute the approximation to f(A)x = Vm*ql2(S, A, Vm, b)
    result = V * h;
end