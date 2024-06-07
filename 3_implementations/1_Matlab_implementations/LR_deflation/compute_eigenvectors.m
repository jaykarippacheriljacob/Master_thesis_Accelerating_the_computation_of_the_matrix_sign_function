% function [Lm, Rm, D] = compute_eigenvectors(A, m)
%     % Compute right eigenvectors using 'LM' (Largest Magnitude) option
%     [Rm, D] = eigs(A, m, 'LM'); % Compute right eigenvectors and eigenvalues
% 
%     % Compute left eigenvectors
%     if isreal(A)  % If A is real symmetric
%         Lm = inv(Rm); % Compute left eigenvectors
%     else  % If A is complex Hermitian
%         Lm = pinv(Rm); % Compute left eigenvectors
%     end
% end
function [Lm, Rm, D] = compute_eigenvectors(A, m)
    % Compute right eigenvectors using 'LM' (Largest Magnitude) option
    [V, D] = eig(A);

    % Extract eigenvalues and sort them in descending order
    [D, idx] = sort(diag(D), 'descend');

    % Rearrange eigenvectors according to the sorted eigenvalues
    Rm = V(:, idx);

    % Compute left eigenvectors
    if isreal(A)  % If A is real symmetric
        Lm = inv(Rm); % Compute left eigenvectors
    else  % If A is complex Hermitian
        Lm = pinv(Rm); % Compute left eigenvectors
    end
    Rm = Rm(1:m, 1:m);
    Lm = Lm(1:m, 1:m);

    disp(Lm);
    disp(Rm);
end