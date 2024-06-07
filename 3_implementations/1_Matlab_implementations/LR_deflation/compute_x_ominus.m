function result = compute_x_ominus(Rm, Lm, x)
    % Compute x_ominus = (1 − Rm * conj(Lm') ) * x

    % Compute conjugate transpose of Lm
    Lm_conj_transpose = conj(Lm');

    % Compute Rm * conjugate transpose of Lm
    RmLm_conj_transpose = Rm * Lm_conj_transpose;

    % Compute (1 - RmL†m) * x
    result = (eye(size(x)) - RmLm_conj_transpose) * x;
end
