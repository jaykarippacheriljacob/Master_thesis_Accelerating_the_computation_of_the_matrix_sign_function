function result = compute_x_ominus(Rm, Lm, x)
    % Compute x_ominus = (1 − Rm * conj(Lm') ) * x

    % Compute conjugate transpose of Lm
    %Lm_conj_transpose = Lm';

    % Compute Rm * conjugate transpose of Lm
    %RmLm_conj_transpose = Rm * Lm';

    % Compute (1 - RmL†m) * x
    %result = (ones(size(x)) - RmLm_conj_transpose) * x;
    %result = (ones(size(x)) - Rm * Lm') * x;
    result = x - Rm * (Lm' * x);
end
