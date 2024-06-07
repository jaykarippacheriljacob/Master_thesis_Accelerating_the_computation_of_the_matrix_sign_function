% wrong......... Need to resolve it.
function [Sm, Tm, D] = computeSchur(A, m)
    % Compute all eigenvalues and eigenvectors of A
    [V, D] = eig(A);

    % Sort eigenvalues and eigenvectors based on magnitude
    [D_sorted, idx] = sort(diag(D), 'descend');
    V_sorted = V(:, idx);

    % Select the m critical eigenvalues and eigenvectors
    D = diag(D_sorted(1:m));
    V = V_sorted(:, 1:m);

    % Construct the Schur matrix Sm
    Sm = V;

    % Construct the upper triangular matrix Tm = S†mASm
    Tm = Sm' * A * Sm;
end
