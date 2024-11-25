function result = compute_x_ominus(Rm, Lm, x)
    %% Computes x_ominus = (1 − Rm * Lm' ) * x
    % Inputs:
    %       Rm     - Matrix Rm.
    %       Lm     - Matrix Lm.
    %       x      - Vector x.
    % Outputs:
    %       result - Resulting vector.
     
    result = x - Rm * (Lm' * x);
end
