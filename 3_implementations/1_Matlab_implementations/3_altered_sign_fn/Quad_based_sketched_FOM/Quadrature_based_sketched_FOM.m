function [fA_b, cost] = Quadrature_based_sketched_FOM(A, b, k_values, s, tol)
    %% Quadrature based sketched FOM approximation for f(A)b.
    % Input: 
    %      A        - N x N matrix
    %      b        - N x 1 vector
    %      k_values - vector of increasing subspace dimensions. We return the deflated 
    %                 sketched Arnoldi matrix approximation for all values in k, k < min(N, s)
    %      s        - sketch matrix row dimension
    %      tol      - parameter for setting the tolerance of the quadrature
    %                 rule
    % Output: 
    %      fA_b     - Approximation of f(A)b for k_values
    %      cost     - No.of matrix{A} vector multiplications
    
    addpath(fullfile(pwd, 'Quad_based_sketched_FOM'));

    n = size(A,2);
    no_k = length(k_values);
    fA_b = zeros(n,no_k);
    kmax = max(k_values);
    cost = zeros(no_k, 1);

    % Setting the tolerance of the quadrature rule if not provided.
    if nargin < 5
        tol = 1e-10;
    end
    
    %% Construct a sketched basis for the Krylov subspace using the arnoldi process
    %  We do the sketeched Arnoldi process only once and then extract the sketched Arnoldi approximations 
    %  for the different dimensions in k_values
    [H, V, ~, SV, SAV, Sb] = QS_FOM_process(A, b, kmax, s);
    
    % Checking whether the arnoldi relation is fulfilled.
    % A*Vm - Vm*Hm ~= h(m+1,m)*q(m+1)*e(m).' *******
    verify_Arnoldi_4_sqr(A, V, H);

    for l=1:no_k  %l-th column holds approx. for subspace dimension k(l)
        fA_b(:,l) = QS_FOM_approx(V, H, SV, SAV, Sb, k_values(l), tol);
        cost(l) = 1 + 2 * k_values(l);
    end 
    
end