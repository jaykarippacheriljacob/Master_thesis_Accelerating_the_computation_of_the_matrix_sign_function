function f_Tm_column = compute_f_Tm_column(Tm)
    % Initialize S0 to identity matrix of the same size as Tm
    S0 = eye(size(Tm));

    % Iterate Roberts' iterative method
    S = S0;
    for i = 1:size(Tm, 1)
        S = 0.5 * (S + inv(S));
    end

    % Compute the first column of f(Tm)
    f_Tm_column = S(:, 1);
end
