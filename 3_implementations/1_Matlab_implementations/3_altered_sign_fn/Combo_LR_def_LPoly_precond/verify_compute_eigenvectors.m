function verify_compute_eigenvectors(B, m, symm_type)
    %% Verify and compute eigenvectors of a matrix.
    % Inputs:
    %   B         - Input matrix
    %   m         - Number of critical eigenvalues to compute
    %   symm_type - Type of symmetric transformation (1 for (B * B'), 0 for (B + B'))
    
    % Check if symm_type is provided, otherwise set default value
    if nargin < 3
        symm_type = 1;
    end

    % make it symmetric
    if symm_type == 1
        A = (B * B');
    else
        A = (B + B');
    end

    % Compute eigenvectors using compute_eigenvectors
    %m = size(A, 1);
    [Rm, Lm, ~] = compute_eigenvectors(A, m);
    
    % Check if left and right eigenvectors are identical
    rel_error = norm(abs(Rm) - abs(Lm), 'fro')/norm(Rm,'fro'); 
    
    % Set tolerance level
    tol = 1e-10;
    
    % Check if the relative error is sufficiently small
    if rel_error < tol
        disp(['Left and right eigen vectors are properly computed for ', num2str(m), ' critical eigen values.']);
    else
        disp(['Left and right eigen vectors are not properly computed for ', num2str(m), ' dimension square matrix.']);
    end
end
