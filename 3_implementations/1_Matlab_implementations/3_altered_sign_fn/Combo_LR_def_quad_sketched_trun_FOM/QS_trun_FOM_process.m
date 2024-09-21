function [H, V, beta, SV, SAV, Sb] = QS_trun_FOM_process(A, b, m, s, trunc)
    % Input: 
    %      A     - N x N matrix
    %      b     - N x 1 vector
    %      m     - Maximum number of iterations for the krylov's subspace, m < min(N, s)
    %      s     - sketch matrix row dimension
    %      trunc - Truncate orthogonalization to the last 'trunc' vector
    % Output: 
    %      H    - Upper Hessenberg matrix from Arnoldi process
    %      V    - Arnoldi vector matrix
    %      sV   - Sketched arnoldi vector matrix
    %      SAV  - Sketched Arnoldi matrix-vector products
    %      Sb   - Sketched vector of b

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
    r0 = A*b;

    Sb = S(r0);

    %% Step 2: Generate basis Vm of Km(A, b), as well as SVm
    %          and SAVm.
    beta = norm(r0);
    V(:, 1) = r0 / beta;

    H = zeros(N, m+1);

    for j = 1:m
        % Apply matrix A to the last basis vector
        v = A * (A * V(:, j));

        SV(:, j) = S(V(:, j));
        SAV(:,j) = S(v);

        % Truncated arnoldi process: Orthogonalization
        for i = max(1, j - trunc + 1):j
            % Compute the coefficient for orthogonalization
            H(i, j) = V(:, i)' * v;

            % Subtract the projection of v onto the already computed basis vectors
            v = v - H(i, j) * V(:, i);
        end

        % Store the new basis vector
        H(j+1,j) = norm(v);
        V(:, j+1) = v / H(j+1, j);
    end
end