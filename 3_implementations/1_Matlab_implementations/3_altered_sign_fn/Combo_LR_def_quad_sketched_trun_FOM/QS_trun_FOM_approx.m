function [fA_b] = QS_trun_FOM_approx(V, H, SV, SAV, Sb, m, tol)
    %% returns the k-th Arnoldi matrix function approximation to the inverse square
    %  root applied to the vector beta*V(:,1)
    % Input:
    %      V    - Arnoldi basis matrix (size Nx(maxIter+1))
    %      H    - Upper Hessenberg matrix from Arnoldi process (size maxIter x maxIter)
    %      SV   - Sketched Arnoldi basis vectors (size sketchDim x maxIter)
    %      SAV  - Sketched matrix-vector products (size sketchDim x maxIter)
    %      Sb   - Sketched vector of b (size sketchDim x 1)
    %      m    - Number of iterations used for the Krylov subspace
    %      tol  - Tolerance for quadrature rule
    %
    % Output:
    %      fA_b - Approximation of f(A)b

    n = size(H,2);  
    if m > n
        disp('Error: subspace too small');
    else
        V = V(:, 1:m);
    
        %% Step 1: Compute thin QR decomposition SVm = Qm * Rm
        [Qm, Rm] = qr(SV(:,1:m), 0);
    
        %% Step 2: Compute SVm = Qm, SAVm = (SAVm)*inv(Rm), Vm = Vm*inv(Rm)
        SV = Qm;
        SAV = SAV(:, 1:m) / Rm;
        V = V / Rm;
    
        %% Step 3: Compute the quadrature rules ql1(S, A, Vm, b) and ql2(S, A, Vm, b)
        [~, ~, ~, h] = Quadrature_rule_invsqrt(SV, SAV, Sb, tol);
    
        %% Step 4: Compute the approximation to f(A)x = Vm*ql2(S, A, Vm, b)
        fA_b = V * h;
    end
end