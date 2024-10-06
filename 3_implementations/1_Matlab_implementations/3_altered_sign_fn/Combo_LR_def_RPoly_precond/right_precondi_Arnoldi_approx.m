function fA_b = right_precondi_Arnoldi_approx(V,H,beta,k)
    %% returns the k-th Arnoldi matrix function approximation to the inverse square
    %  root applied to the vector beta*V(:,1)
    % Input: 
    %         V    - matrix of Arnoldi vectors, size is n x (m+1)
    %         H    - matrix of orthogonalization coefficients, size is (m+1) x m
    %         k    - number of approximation, k <= m
    %         beta - see above
    % Output: 
    %         fA_b - the Arnoldi approximation on the Krylov subspace of dimension k 
 
 m = size(H,2);  
 if k > m
     disp('Error: subspace too small');
 else
    e1 = zeros(k,1); 
    e1(1) = 1;
    fA_b = V(:,1:k) * inv(sqrtm(H(1:k, 1:k))) * e1 * beta;
 end


end

