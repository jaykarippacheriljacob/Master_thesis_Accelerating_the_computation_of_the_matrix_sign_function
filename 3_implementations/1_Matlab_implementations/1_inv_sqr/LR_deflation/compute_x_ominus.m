function result = compute_x_ominus(Rm, Lm, x)
    % Inputs:
    %       Rm: Matrix Rm
    %       Lm: Matrix Lm
    %       x: Vector x
    % Outputs:
    %       result: Resulting vector
    
    % Compute x_ominus = (1 − Rm * Lm' ) * x
    result = x - Rm * (Lm' * x);
end
