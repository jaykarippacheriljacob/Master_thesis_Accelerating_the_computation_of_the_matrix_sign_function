function [fA_b, cost] = lr_deflation(A, x, m_values, k_values)
    % Input:
    %      A        - n x n matrix
    %      x        - n x 1 vector
    %      m_values - vector of increasing no. of critical values for which defation has to be undergone.
    %      k_values - vector of increasing Krylov subspace dimensions. 
    % Output: 
    %      fA_b     - Approximation of f(A)b for k_values
    %      cost     - No.of matrix{A} vector multiplications

    addpath(fullfile(pwd, 'LR_deflation'));
    
    %% Initializing variables
    no_k = length(k_values);
    no_m = length(m_values);
    n = size(A,2);
    fA_x_ominus = zeros(n,no_k);
    kmax = max(k_values);
    cost = zeros(no_k*no_m, 1);
    fA_b = zeros(n, no_k*no_m);
    j = 1;
    
    %% Step 1: Compute left and right eigenvectors
    mmax = max(m_values);
    [Rmmax, Lmmax, Dmmax] = compute_eigenvectors(A, mmax);

    for i = 1:length(m_values)
        m = m_values(i);
        Rm = Rmmax(:, 1:m);
        Lm = Lmmax(:, 1:m);
        Dm = Dmmax(1:m, 1:m);

        %% Step 2: Compute f(A) for critical eigenvalues
        f_Tm = compute_sign_function_diag(Dm);
    
        %% Step 3: Compute x_ominus = (1 − Rm * Lm' ) * x 
        % ' -> transpose conjugate
        % .' -> transpose
        x_ominus = compute_x_ominus(Rm, Lm, x);
    
        %% Step 4: Construct an orthonormal basis for the Krylov subspace using the arnoldi process
        %          We do the Arnoldi process only once and then extract the Arnoldi approximations 
        %          for the different dimensions in k_values 
        [H, V, beta] = Arnoldi_process(A, x_ominus, kmax);
    
        % Checking whether the arnoldi relation is fulfilled.
        % A*Vm - Vm*Hm ~= h(m+1,m)*q(m+1)*e(m).' *******
        verify_Arnoldi(A, V, H);
        
        %% Step 5: Compute f(A)x_ominus
        for l=1:no_k  %l-th column holds approx. for subspace dimension k(l)
            fA_x_ominus(:,l) = Arnoldi_approx(V,H,beta,k_values(l));
            cost(j+l-1) = k_values(l);
        end
        
        %% Step 6: Compute the approximation to f(A)x
        fA_b(:,j:(j+no_k-1)) = Rm * (f_Tm * (Lm' * x))*ones(1,no_k) + fA_x_ominus;
        j = j+no_k;
    end
end
