function [Vm, Wm, D1] = compute_eigenvectors(A, m)
    % Inputs:
    %       A: The matrix for which eigenvalues and eigenvectors are computed.
    %       m: The number of eigenvalues and eigenvectors to compute.
    % Outputs:
    %       Vm: Right eigenvectors of A.
    %       Wm: Left eigenvectors of A.
    %       D1: Eigenvalues corresponding to the right eigenvectors.

    A = A / normest(A);

    % Compute Right eigenvalues and eigenvectors for A
    [Vm, D1] = eigs(A, m, 'smallestabs');
    
    % Compute Left eigenvalues and eigenvectors for A
    if min(size(A,1), m+1) == m+1
        [Wm, D2] = eigs(A', m+1, 'smallestabs');
    else
        [Wm, D2] = eigs(A', m, 'smallestabs');
    end

    % disp(D1);
    % disp(D2);

    % Evaluation of the eigen values to see if the left and right
    % eigenvectors are ordered in pairs
    tol = 1e-8; % Tolerance for comparing floating point values
    for i = 1:size(D1, 1)
        
        % Check if D2m(i,i) is the conjugate of D1m(i,i)
        if abs(conj(D1(i,i)) - D2(i,i)) > tol
            %disp(['Conjugate eigenvalue correction needed for eigenvalue ' num2str(i)]);
            
            % Find the conjugate of D1m(i,i) from the elements inside D2m
            for j = i+1:size(D2, 1)
                if abs(conj(D1(i, i)) - D2(j,j)) < tol
                    
                    % Swap places
                    %D2(i, i) = D2(i, i) + D2(j, j);
                    %D2(j, j) = D2(i, i) - D2(j, j);
                    %D2(i, i) = D2(i, i) - D2(j, j);
                    temp = D2(i, i);
                    D2(i, i) = D2(j, j);
                    D2(j, j) = temp;

                    
                    % Swap places for corresponding eigenvectors
                    temp = Wm(:, j);
                    Wm(:, j) = Wm(:, i);
                    Wm(:, i) = temp;
                    %Wm(:, [i, j]) = Wm(:, [j, i]);
                    break;
                end
            end
        end
    end

    % Remove the last column of Wm
    Wm = Wm(:, 1:m);
    %Wm = Wm(:, 1:end-1);
    % disp(D1);
    % disp(D2);
    % disp(Wm);
end