function f_Dm = compute_sign_function_diag(Dm)
    %% Computes the sign function applied to the diagonal elements 
    % Inputs: 
    %       Dm   - Input matrix
    % Outputs:
    %       f_Dm - Resultant of the sign function applied to the diagonal elements 
    
    % sign(A) = A(A^{2})^{-1/2}

    % Conversion of the diagonal elements of the matrix Dm to a vector
    d = diag(Dm);
    % disp(d);

    % Computation of the sign function for the matrix Dm
    f_Dm = diag(d./sqrt(d.*d));
    % disp(f_Dm)

end