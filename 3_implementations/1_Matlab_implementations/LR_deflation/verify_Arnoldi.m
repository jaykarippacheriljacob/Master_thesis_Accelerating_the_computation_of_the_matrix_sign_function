function verify_Arnoldi(A, Vm, Hm, hm1m, qm1)
% This function checks whether the Arnoldi relation is fulfilled for
% matrices A, Vm, and Hm using the Arnoldi relation equation, 
% A*Vm - Vm*Hm ~= h(m+1,m)*q(m+1)*e(m).'
% A*Vm = V(n, m+1) * H(m+1, m)

%   Inputs:
%       A: The matrix A.
%       Vm: The matrix Vm.
%       Hm: The matrix Hm.
%       hm1m: The (m+1,m)-th entry of Hm.
%       qm1: The (m+1)-th column of Vm.
%
%   Outputs:
%       None. Prints a message indicating whether the Arnoldi relation is
%       verified or not.

    % Compute the left-hand side of the Arnoldi relation
    lhs = A * Vm - Vm * Hm;
    %disp(lhs);
    
    % Compute the right-hand side of the Arnoldi relation using the additional equation
    m = size(Vm, 2);
    em = zeros(m, 1);
    em(m) = 1;
    rhs = hm1m * qm1 * em.';
    %disp(rhs);
    
    % Compute the residual
    residual = norm(lhs - rhs, 'fro')/norm(lhs,'fro');
    disp(residual)
    
    % Set tolerance level
    tol = 1e-10;
    
    % Check if the residual is sufficiently small
    if residual < tol
        disp('Arnoldi relation verified.');
    else
        disp('Arnoldi relation not verified.');
    end
end
