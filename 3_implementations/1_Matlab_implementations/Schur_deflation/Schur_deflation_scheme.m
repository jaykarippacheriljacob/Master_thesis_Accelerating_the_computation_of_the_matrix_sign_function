% wrong......... Need to resolve it.
function result = Schur_deflation_scheme(A, x, m, k)
    
    % Step 1: Compute the corresponding Schur matrix Sm and uppertrianglar
    % matrix, Tm = Sm' * A * Sm
    [Sm, ~] = eigs(A, m, 'smallestabs');
    Tm = Sm' * (A * Sm);

    % Step 2: Compute the triangular matrix f(Tm)
    f_Tm = compute_sign_function(Tm);

    % Step 3: Compute x_perp = (1 - Sm*Sm')*x
    x_perp = x - (Sm * (Sm' * x));

    % Step 4: Construction of an orthogonal basis for hte projected Krylov
    % subspace Kk_perp(A, x_perp).
    [Hk, Vk] = Arnoldi_method(A, x_perp, k);

    % Step 5: Computing f(H) and solving the Sylvester equation for Y
    f_Hk = compute_sign_function(Hk(1:k, 1:k));
    Y = sylvester_equation(A, Sm, Tm, Hk(1:k, 1:k), Vk(:, 1:k), f_Tm, f_Hk);

    % Compute f(A) * x
    result = Sm * (f_Tm * (Sm' * x)) + [Sm, Vk(:, 1:k)] * ([Y; f_Hk] * (Vk(:, 1:k)' * x));
end
