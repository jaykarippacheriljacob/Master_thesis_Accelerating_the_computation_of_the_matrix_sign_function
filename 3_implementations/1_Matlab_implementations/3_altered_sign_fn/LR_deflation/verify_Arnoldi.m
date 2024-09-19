function verify_Arnoldi(A, V, H)
% This function checks whether the Arnoldi relation is fulfilled for
% matrices A^2, Vm, and Hm using the Arnoldi relation equation, 
% A^2*V(:,1:m) - V(:,1:m+1)*H(1_m+1,1:m) should be zero in exact arithmetic

%   Inputs:
%       A: The matrix A.
%       V: matrix of Arnoldi vectors
%       H: matrix of orthogonalizaiion coeffs
%
%   Outputs:
%       None. Prints a message indicating whether the Arnoldi relation is
%       verified or not.

    m = size(V,2)-1;
    em = zeros(m, 1);
    em(m) = 1;
    lhs = (A*V(:,1:m))-V(:,1:m)*H(1:m, 1:m);
    rhs = H(m+1, m) * V(:, m+1) * em.';
    rel_error = norm(lhs - rhs, 'fro')/norm(lhs,'fro'); % since V is orthogonal, we can use norm(H) instead of norm(V*H)
    
    % Display the error
    disp(['Relative error in Arnoldi relation= ', num2str(rel_error)]);
    
end