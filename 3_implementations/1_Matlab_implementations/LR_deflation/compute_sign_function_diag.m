
function f_Dm = compute_sign_function_diag(Dm)
    % Inputs: 
    %       Dm: Input matrix
    % Outputs:
    %       f_Dm: Resultant sign function applied to the diagonal elements 
    
    % sign(A) = A(A^{2})^{-1/2}
    f_Dm = speye(size(Dm, 1));
    for i = 1:size(Dm, 1)
        z = Dm(i,i);
        f_Dm(i, i) = z * ((sqrt(z * z))^(-1));
    end
end